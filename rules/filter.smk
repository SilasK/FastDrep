

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
            filtered= Q.query(params.filter).index



            os.makedirs(output[0])
            for f in filtered:
                os.symlink(os.path.join(input.dir,f),
                           os.path.join(output[0],f))
