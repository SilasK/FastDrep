

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
            assert not Q.index.duplicated().any()

            files_in_folder= set(os.listdir(input.dir))
            intersection = Q.index.intersection(files_in_folder)

            missing_quality= files_in_folder.difference(intersection)
            if len(missing_quality) >0:
                logger.info("missing quality information for following files:"
                            "\n".join(missing_quality))


            missing_fasta= Q.index.difference(intersection)
            if len(missing_fasta) >0:
                logger.info("missing fasta file for following genomes:"
                            "\n".join(missing_fasta))

            Q= Q.loc[intersection]


            filtered= Q.query(params.filter).index



            os.makedirs(output[0])
            for f in filtered:
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
