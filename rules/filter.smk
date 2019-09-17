
rule calculate_stats:
    input:
        input_genome_folder,
    output:
        "filter/genome_stats.tsv"
    threads:
        10
    run:
        from common.genome_stats import get_genome_stats
        from common import genome_pdist as gd
        from multiprocessing import Pool

        import pandas as pd


        pool = Pool(threads)

        fasta_files= glob(f"{input[0]}/*.fasta")

        results= pool.map(get_genome_stats,fasta_files)
        Stats= pd.DataFrame(results,columns=['Length','N50'],index=fasta_files)
        Stats.index= gd.simplify_index(Stats.index)

        Stats.to_csv(output[0],sep='\t')

if config.get('skip_filter',False):
    filter_genome_folder= input_genome_folder
else:
    filter_genome_folder='filter/bins_filtered_quality'




    localrules: filter_genomes_by_size, filter_genomes_by_quality,get_orig_filenames

    rule get_orig_filenames:
        input:
            dir=input_genome_folder
        output:
            "filter/bin2filename.tsv"
        run:
            import pandas as pd
            filenames= pd.Series(os.listdir(input.dir),name='Filename')
            filenames.index= filenames.apply(lambda f: os.path.splitext(f)[0])

            filenames.index.name='Bin'

            filenames.to_csv(output[0],sep='\t',header=True)


    rule filter_genomes_by_size:
        input:
            dir=os.path.abspath(input_genome_folder),
            genome_stats= rules.calculate_stats.output[0],
            filenames= rules.get_orig_filenames.output[0]
        output:
            directory(temp('filter/bins_filtered_size'))
        params:
            genome_filter=config['filter_criteria']
        run:
            import pandas as pd

            genome_stats= pd.read_csv(input.genome_stats,index_col=0,sep='\t')
            filenames= pd.read_csv(input.filenames,index_col=0,sep='\t',squeeze=True)



            filtered=  filenames.loc[genome_stats.query(params.genome_filter).index]

            if filtered.shape[0]<2:
                raise Exception("Less than 2 genomes pass fragementation filter")

            os.makedirs(output[0])
            for f in filtered.values:
                os.symlink(os.path.join(input.dir,f),
                           os.path.join(output[0],f))


    localrules: get_predifined_quality
    ruleorder: get_predifined_quality> merge_checkm
    rule get_predifined_quality:
        input:
            config['genome_qualities']
        output:
            "filter/Genome_quality.tsv"
        shell:
            "cp {input} {output}"





    rule filter_genomes_by_quality:
        input:
            filtered_dir= rules.filter_genomes_by_size.output[0],
            dir= os.path.abspath(input_genome_folder),
            quality=config['genome_qualities'],
        output:
            directory(filter_genome_folder)
        params:
            quality_filter=config['qualityfilter_criteria'],
        run:
            import pandas as pd


            Q= gd.load_quality(input.quality)
            assert not Q.index.duplicated().any(), f"duplicated indexes in {input.quality}"

            Q['quality_score']= Q.eval(config['quality_score'])



            files_in_folder= pd.Series(os.listdir(input.filtered_dir))

            files_in_folder.index= files_in_folder.apply(lambda f: os.path.splitext(f)[0])

            intersection = Q.index.intersection(files_in_folder.index)

            missing_quality= files_in_folder.index.difference(intersection)
            if len(missing_quality) >0:
                logger.error(f"missing quality information for following files: {missing_quality}")

            missing_fasta= Q.index.difference(intersection)
            if len(missing_fasta) >0:
                logger.error(f"missing fasta file for following genomes: {missing_fasta}")


            Q= Q.loc[intersection]

            filtered=  files_in_folder.loc[Q.query(params.quality_filter).index]


            if filtered.shape[0]<2:
                raise Exception("Less than 2 genomes pass quality filter")

            os.makedirs(output[0])
            for f in filtered.values:
                os.symlink(os.path.join(input.dir,f),
                           os.path.join(output[0],f))



def gen_names_for_range(N,prefix='',start=1):
    """generates a range of IDS with leading zeros so sorting will be ok"""
    n_leading_zeros= len(str(N))
    format_int=prefix+'{:0'+str(n_leading_zeros)+'d}'
    return [format_int.format(i) for i in range(start,N+start)]


genome_folder='mags'

localrules: rename_genomes, decompress_genomes, rename_quality
checkpoint rename_genomes:
    input:
        genome_folder= filter_genome_folder,
        stats= rules.calculate_stats.output[0],

    output:
        genome_folder= directory(genome_folder),
        mapping= "tables/renamed_genomes.tsv",
        stats="tables/genome_stats.tsv",
    params:
        prefix= config.get('mag_prefix','MAG')

    run:

        import pandas as pd
        import shutil
        from glob import glob

        Mapping= pd.DataFrame()

        Mapping['Original_fasta'] = glob(os.path.join(input.genome_folder,"*.f*"))
        Mapping.index = gd.simplify_index(Mapping.Original_fasta)
        Mapping.index.name='Original'
        Mapping.sort_index(inplace=True)

        Mapping['Genome']= gen_names_for_range(Mapping.shape[0],params.prefix)
        Mapping.to_csv(output.mapping,sep='\t')

        #Rename stats
        Stats= pd.read_csv(input.stats, sep='\t',index_col=0)
        Stats= Stats.rename(index=Mapping.Genome).loc[Mapping.Genome]
        Stats.to_csv(output.stats,sep='\t')




        os.makedirs(output.genome_folder)

        for _,row in Mapping.iterrows():

            shutil.copy(row.Original_fasta,
                        os.path.join(output.genome_folder,row.Genome+'.fasta')
                        )

rule rename_quality:
    input:
        mapping= rules.rename_genomes.output.mapping,
        quality=config['genome_qualities']
    output:
        quality="tables/Genome_quality.tsv"

    run:
        import pandas as pd

        Mapping= pd.read_csv(input.mapping,sep='\t',index_col=0)

        Q= gd.load_quality(input.quality)
        Q= Q.rename(index=Mapping.Genome).loc[Mapping.Genome]
        Q.to_csv(output.quality,sep='\t')



rule decompress_genomes:
    input:
        "genomes.tar.gz"
    output:
        directory("genomes")
    shell:
        "tar -xzf {input}"
ruleorder: decompress_genomes>rename_genomes
