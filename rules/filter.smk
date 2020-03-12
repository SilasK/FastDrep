
localrules: get_orig_filenames

rule get_orig_filenames:
    input:
        dir=input_genome_folder
    output:
        filenames="filter/filenames_all.tsv",
        dir=directory(temp('filter/all_bins'))
    run:
        import pandas as pd
        from glob import glob
        filenames= pd.DataFrame()
        filenames['FilenameOriginal']=glob(os.path.join(input.dir,"*.f*"))



        filenames.index= filenames.FilenameOriginal.apply(io.simplify_path)
        filenames.index.name='Bin'
        filenames['Filename']= filenames.index.map(lambda g: os.path.join(output.dir,g+'.fasta'))
        filenames.to_csv(output.filenames,sep='\t')

        os.makedirs(output.dir)
        for g in filenames.index:
            os.symlink(filenames.loc[g,'FilenameOriginal'],
                       filenames.loc[g,'Filename'])






rule calculate_stats:
    input:
        rules.get_orig_filenames.output.filenames,
    output:
        "filter/genome_stats.tsv"
    threads:
        config['threads']
    run:
        from common.genome_stats import get_many_genome_stats
        import pandas as pd
        filenames= pd.read_csv(input[0],sep='\t',index_col=0,squeeze=True)
        get_many_genome_stats(filenames.FilenameOriginal,output[0],threads)


if config.get('skip_filter',False):
    filtered_filenames= "filter/filenames_all.tsv"
else:
    filtered_filenames="filter/filenames_filtered_quality.tsv"


    localrules: filter_genomes_by_size, filter_genomes_by_quality
    rule filter_genomes_by_size:
        input:
            dir=rules.get_orig_filenames.output.dir,
            genome_stats= rules.calculate_stats.output[0],
            filenames= rules.get_orig_filenames.output.filenames
        output:
            dir=directory(temp('filter/bins_filtered')),
            filenames=temp("filter/filenames_filtered_size.tsv")
        params:
            genome_filter=config['filter_criteria']
        run:
            import pandas as pd

            genome_stats= pd.read_csv(input.genome_stats,index_col=0,sep='\t')
            filenames= pd.read_csv(input.filenames,index_col=0,sep='\t',squeeze=True)



            filtered=  filenames.loc[genome_stats.query(params.genome_filter).index]
            filtered.to_csv(output.filenames,sep='\t')

            if filtered.shape[0]<2:
                raise Exception("Less than 2 genomes pass fragementation filter")

            os.makedirs(output.dir)
            symlink_relative(filtered.index+'.fasta',input.dir,output.dir)




    localrules: get_predefined_quality, combine_checkm_quality


    if 'genome_qualities' in config:
        ruleorder: get_predefined_quality> merge_checkm
        rule get_predefined_quality:
            input:
                config['genome_qualities']
            output:
                "filter/Genome_quality.tsv"
            shell:
                "cp {input} {output}"





    rule filter_genomes_by_quality:
        input:
            filenames=rules.filter_genomes_by_size.output.filenames,
            filtered_dir= rules.filter_genomes_by_size.output.dir,
            quality="filter/Genome_quality.tsv",
            stats= "filter/genome_stats.tsv"
        output:
            filenames=temp(filtered_filenames)
        params:
            quality_filter=config['qualityfilter_criteria'],
        run:
            import pandas as pd
            from glob import glob
            from numpy import log


            Q= gd.load_quality(input.quality)
            assert not Q.index.duplicated().any(), f"duplicated indexes in {input.quality}"

            stats= pd.read_csv(input.stats,index_col=0,sep='\t')
            stats['logN50']=log(stats.N50)

            Q=Q.join(stats.loc[Q.index,stats.columns.difference(Q.columns)])

            Q['quality_score']= Q.eval(config['quality_score'])

            Filenames= pd.read_csv(input.filenames,sep='\t',index_col=0)


            intersection = Q.index.intersection(Filenames.index)

            missing_quality= Filenames.index.difference(intersection)
            if len(missing_quality) >0:
                raise Exception(f"missing quality information for following files: {missing_quality}")

            #missing_fasta= Q.index.difference(intersection)
            #if len(missing_fasta) >0:
            #    raise Exception(f"missing fasta file for following genomes: {missing_fasta}")



            Q= Q.loc[intersection]

            filtered=  Filenames.loc[Q.query(params.quality_filter).index]
            filtered.to_csv(output.filenames,sep='\t')


            if filtered.shape[0]<2:
                raise Exception("Less than 2 genomes pass quality filter")





def gen_names_for_range(N,prefix='',start=1):
    """generates a range of IDS with leading zeros so sorting will be ok"""
    n_leading_zeros= len(str(N))
    format_int=prefix+'{:0'+str(n_leading_zeros)+'d}'
    return [format_int.format(i) for i in range(start,N+start)]





localrules: rename_genomes, decompress_genomes, rename_quality
checkpoint rename_genomes:
    input:
        filenames= filtered_filenames,
        stats=rules.calculate_stats.output[0]
    output:
        genome_folder= directory(genome_folder),
        mapping= "tables/renamed_filenames.tsv",
        stats="tables/genome_stats.tsv"
    params:
        prefix= config.get('mag_prefix','MAG'),
        method= config.get('rename_method','prefix')

    run:

        import pandas as pd
        import shutil
        from glob import glob

        Mapping= pd.read_csv(input.filenames,sep='\t',index_col=0)

        if (params.method is None) or (params.method=='None'):
            Mapping['Genome']=Mapping.index
        elif params.method=='prefix':
            Mapping['Genome']= gen_names_for_range(Mapping.shape[0],params.prefix)
        else:
            raise NotImplementedError("config['rename_method'] should be one of 'None', 'prefix'")


        # new filenames
        Mapping['Filename']= Mapping.Genome.apply(lambda g: os.path.join(output.genome_folder,g+'.fasta'))
        Mapping.to_csv(output.mapping,sep='\t')
        #newFilenames= Mapping.set_index('Genome').Filename


        os.makedirs(output.genome_folder)

        for _,row in Mapping.iterrows():

            shutil.copy(row.FilenameOriginal,
                        row.Filename
                        )
        #Rename stats
        Stats= pd.read_csv(input.stats, sep='\t',index_col=0)
        Stats= Stats.rename(index=Mapping.Genome).loc[Mapping.Genome]
        Stats.to_csv(output.stats,sep='\t')

#TODO: joun logN50 to quality for cluster_species
rule rename_quality:
    input:
        mapping= rules.rename_genomes.output.mapping,
        quality="filter/Genome_quality.tsv",
    output:
        quality="tables/Genome_quality.tsv",


    run:
        import pandas as pd

        Mapping= pd.read_csv(input.mapping,sep='\t',index_col=0)

        Q= gd.load_quality(input.quality)
        Q= Q.rename(index=Mapping.Genome).loc[Mapping.Genome]
        Q.to_csv(output.quality,sep='\t')






# rule decompress_genomes:
#     input:
#         "genomes.tar.gz"
#     output:
#         directory("genomes")
#     shell:
#         "tar -xzf {input}"
# ruleorder: decompress_genomes>rename_genomes
