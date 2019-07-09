

if config.get('skip_filter',False):
    filter_genome_folder= genome_folder
else:
    filter_genome_folder='filtered_bins'

    localrules: filter_genomes
    rule filter_genomes:
        input:
            dir=os.path.abspath(genome_folder),
            quality=config['genome_qualities']
        output:
            directory(filter_genome_folder)
        params:
            filter=config['filter_criteria']
        run:
            import pandas as pd

            Q= pd.read_csv(input.quality,sep='\t',index_col=0)
            assert not Q.index.duplicated().any(), f"duplicated indexes in {input.quality}"

            Q.index= gd.simplify_index(Q.index)

            files_in_folder= pd.Series(os.listdir(input.dir))

            files_in_folder.index= files_in_folder.apply(lambda f: os.path.splitext(f)[0])

            intersection = Q.index.intersection(files_in_folder.index)

            missing_quality= files_in_folder.index.difference(intersection)
            if len(missing_quality) >0:
                logger.error(f"missing quality information for following files: {missing_quality}")

            missing_fasta= Q.index.difference(intersection)
            if len(missing_fasta) >0:
                logger.error(f"missing fasta file for following genomes: {missing_fasta}")


            Q= Q.loc[intersection]


            filtered=  files_in_folder.loc[Q.query(params.filter).index]


            if filtered.shape[0]<2:
                raise Exception("Less than 2 genomes pass quality filter")

            os.makedirs(output[0])
            for f in filtered.values:
                os.symlink(os.path.join(input.dir,f),
                           os.path.join(output[0],f))



rule calculate_stats:
    input:
        filter_genome_folder,
        genome_folder
    output:
        "genome_stats.tsv"
    run:
        from common.genome_stats import get_genome_stats
        import pandas as pd

        Stats= pd.DataFrame(columns=['Total_length','N50'])
        for fasta_file in glob(f"{input[0]}/*.fasta"):
            sample = os.path.splitext(fasta_file)[0].replace(f'{input[0]}/','')
            Stats.loc[sample]= get_genome_stats(fasta_file)
        Stats.to_csv(output[0],sep='\t')
