#TODO: expand precluster for species cluaterin

# duplicated rule in case not imported
rule calculate_stats:
    input:
        "tables/Genome_quality.tsv",
    output:
        "tables/genome_stats.tsv"
    threads:
        config['threads']
    run:


        from common.genome_stats import get_many_genome_stats
        import pandas as pd
        d= pd.read_csv(input[0],sep='\t',index_col=0,squeeze=True,usecols=[0])
        filenames = d.index.map(lambda s: f"genomes/{s}{config['fasta_extension']}")
        del d
        get_many_genome_stats(filenames,output[0],threads)

checkpoint get_alignment_subsets:
    input:
        rules.precluster.output.edgelist
    output:
        dir=temp(directory('alignment/subsets'))
    params:
        N=config['subset_size_alignments'],

    run:

        os.makedirs(output.dir)


        fout=None
        with open(input[0],'r') as edgelist:
            for i,line in enumerate(edgelist):
                if (i % int(params.N)) ==0:
                    n_file= int(i // params.N )+1
                    if fout is not None: fout.close()
                    fout= open(f"{output.dir}/subset_{n_file}.txt",'w')

                pair=sorted(line.strip().split())

                fout.write("\t".join(pair)+'\n')


            fout.close()
        assert i>=1,"No connection correspond the criteria for preclustering"

rule many_minimap:
    input:
        genome_folder=genome_folder,
        alignment_list="alignment/subsets/{subset}.txt",
        genome_stats= "tables/genome_stats.tsv"
    output:
        alignments_stats="alignment/alignments_stats/{subset}.tsv",
    log:
        "logs/minimap/{subset}.txt"
    benchmark:
        "logs/benchmarks/minimap/{subset}.txt"
    threads:
        config['threads']
    resources:
        time=config['runtime']['minimap']
    conda:
        "../envs/minimap2.yaml"
    params:
        minimap_extra= config['minimap_extra'], #asm5/asm10/asm20: asm-to-ref mapping, for ~0.1/1/5% sequence divergence
        paf_folder="alignment/paf",
        extension=config['fasta_extension']
    script:
        "../scripts/many_minimap.py"


def get_combine_paf_input(wildcards):


    subset_dir= checkpoints.get_alignment_subsets.get().output[0]
    logger.info(f"subset dir : {subset_dir}")

    subsets= glob_wildcards(os.path.join(subset_dir,'{subset}.txt')).subset

    return expand("alignment/alignments_stats/{subset}.tsv",subset=subsets)


localrules: combine_paf
rule combine_paf:
    input:
        get_combine_paf_input
    output:
        "tables/minimap_dists.tsv"
    run:

        import pandas as pd

        sep='\t'

        D = pd.read_csv(input[0],sep=sep)
        n_cols= D.shape[1]
        D.to_csv(output[0],sep=sep)

        for file in input[1:]:
            D = pd.read_csv(file,sep=sep)
            assert n_cols== D.shape[1], f"File {file} doen't have the same number of columns as the one before. Cannot concatenate."
            D.to_csv(output[0],sep=sep,header=False,mode='a')


        in_dir= os.path.dirname(input[0])
        import shutil
        shutil.rmtree(in_dir)



## Mummer


def estimate_time_mummer(N,threads):
    "retur time in minutes"

    time_per_mummer_call = 1/60 # 1min in h

    return int(N*time_per_mummer_call)//threads + 5



rule run_mummer:
    input:
        comparison_list="alignment/subsets/{subset}.txt",
        genome_folder= genome_folder,
        genome_stats="tables/genome_stats.tsv",
    output:
        temp("alignment/ANI/{subset}.tsv")
    threads:
        config['threads']
    conda:
        "../envs/mummer.yaml"
    resources:
        time= lambda wc, input, threads: estimate_time_mummer(config['subset_size_alignments'],threads),
        mem=config['mem'].get('mummer',20)
    log:
        "logs/mummer/workflows/{subset}.txt"
    benchmark:
        "logs/benchmarks/mummer/{subset}.txt"
    benchmark:
        "logs/benchmarks/mummer/{subset}.txt"
    params:
        path= os.path.dirname(workflow.snakefile)
    shell:
        "snakemake -s {params.path}/rules/mummer.smk "
        "--config comparison_list='{input.comparison_list}' "
        " genome_folder='{input.genome_folder}' "
        " subset={wildcards.subset} "
        " genome_stats={input.genome_stats} "
        " --rerun-incomplete "
        "-j {threads} --nolock 2> {log}"



def get_merge_mummer_ani_input(wildcards):

    if '.gz' in config['fasta_extension']:
        raise Exception("Mummer doesn't handle gziped files, "
                        "Select in the config file minimap instead"
                        )
    subset_dir= checkpoints.get_alignment_subsets.get().output[0]
    subsets= glob_wildcards(os.path.join(subset_dir,'{subset}.txt')).subset


    return expand("alignment/ANI/{subset}.tsv",subset=subsets)

localrules: merge_mummer_ani
rule merge_mummer_ani:
    input:
        get_merge_mummer_ani_input
    output:
        "tables/mummer_dists.tsv"
    run:
        import pandas as pd
        import shutil

        Mummer={}
        for file in input:
            Mummer[io.simplify_path(file)]= pd.read_csv(file,index_col=[0,1],sep='\t')

        M= pd.concat(Mummer,axis=0)
        M.index= M.index.droplevel(0)
        M.to_csv(output[0],sep='\t')


        ani_dir= os.path.dirname(input[0])
        shutil.rmtree(ani_dir)
        #sns.jointplot('ANI','Coverage',data=M.query('ANI>0.98'),kind='hex',gridsize=100,vmax=200)
