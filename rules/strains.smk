
rule Dstrain:
    input:
        ANI="tables/dist_strains.tsv",
        ani_dir='mummer/ANI',
        delta_dir="mummer/delta"
    shell:
        " tar -czf {input.delta_dir}.tar.gz {input.delta_dir} ;"
        "rm -rf {input.delta_dir} {input.ani_dir}"









def estimate_time_mummer(N,threads):
    "retur time in minutes"

    time_per_mummer_call = 10/60

    return int(N*time_per_mummer_call)//threads + 5

localrules: get_deltadir,decompress_delta,Dstrain
rule get_deltadir:
    output:
        directory("mummer/delta")
    run:
        os.makedirs(output[0])

rule decompress_delta:
    input:
        ancient("mummer/delta.tar.gz")
    output:
        directory("mummer/delta")
    shell:
        "tar -xzf {input}"
ruleorder: decompress_delta>get_deltadir

def get_merge_mummer_ani_input(wildcards):

    subsets=get_mummer_subsets(wildcards)

    return expand("mummer/ANI/{subset}.tsv",subset=subsets)

rule merge_mummer_ani:
    input:
        get_merge_mummer_ani_input
    output:
        "tables/dist_mummer.tsv"
    run:
        import pandas as pd
        Mummer={}
        for file in input:
            Mummer[io.simplify_path(file)]= pd.read_csv(file,index_col=[0,1],sep='\t')

        M= pd.concat(Mummer,axis=0)
        M.index= M.index.droplevel(0)
        M.to_csv(output[0],sep='\t')

        #sns.jointplot('ANI','Coverage',data=M.query('ANI>0.98'),kind='hex',gridsize=100,vmax=200)




rule run_mummer:
    input:
        comparison_list="mummer/subsets/{subset}.txt",
        genome_folder= genome_folder,
        genome_stats="tables/genome_stats.tsv",
        delta_dir="mummer/delta"
    output:
        temp("mummer/ANI/{subset}.tsv")
    threads:
        config['threads']
    conda:
        "../envs/mummer.yaml"
    resources:
        time= lambda wc, input, threads: estimate_time_mummer(config['mummer']['subset_size'],threads),
        mem= 2
    log:
        "logs/mummer/workflows/{subset}.txt"
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


localrules: cluster_strains
checkpoint cluster_strains:
    input:
        dists=rules.merge_mummer_ani.output[0],
        quality ="tables/Genome_quality.tsv",
        mag2species= "tables/mag2species.tsv"
    output:
        mag2strain="tables/mag2strains.tsv"
    params:
        #treshold=config['species_treshold'],
        linkage_method=config.get('linkage_method','average'),
    script:
        "../scripts/group_strains.py"


checkpoint get_ref_bbsplit:
    input:
        dir= genome_folder,
        cluster_file= "tables/mag2strains.tsv"
    output:
        dir= directory("representatives/merged_strains"),
    run:

        import pandas as pd
        df= pd.read_csv(input.cluster_file,sep='\t')


        os.makedirs(output.dir)

        for species_name,sp in df.groupby('Species'):
            with open(os.path.join(output.dir,species_name+'.fasta'),'w') as fasta_out:
                for genome in sp.Representative_Strains.unique():

                    with open(os.path.join(input.dir,genome+'.fasta')) as fasta_in:
                        fasta_out.write(fasta_in.read())
