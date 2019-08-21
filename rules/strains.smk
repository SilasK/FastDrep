

localrules: species_subsets
checkpoint species_subsets:
    input:
        cluster_file="mash/clusters.tsv",
    output:
        subsets_dir= directory("mummer/subsets")
    run:
        import pandas as pd
        labels= pd.read_csv(input[0],sep='\t',header=None,squeeze=True,index_col=0)


        os.makedirs(output.subsets_dir)

        for species in labels.unique():

            genomes_of_cluster= labels.index[labels==species].values
            if len(genomes_of_cluster)>1:

                with open(f"{output.subsets_dir}/{species}.txt","w") as f:

                    f.write(''.join([g+'.fasta\n' for g in genomes_of_cluster ]))


def get_species(wildcards):

    dir=checkpoints.species_subsets.get().output.subsets_dir

    return glob_wildcards(f"{dir}/{{species}}.txt").species


def estimate_time_mummer(input,threads):
    "retur time in minutes"

    N= len(open(input.genome_list).read().split())

    time_per_mummer_call = 10/60

    return int(N**2/2*time_per_mummer_call + N/2)//threads + 5

localrules: get_deltadir,decompress_delta,Dstrain
rule get_deltadir:
    output:
        directory("mummer/delta")
    run:
        os.makedirs(output[0])

rule decompress_delta:
    input:
        "mummer/delta.tar.gz"
    output:
        directory("mummer/delta")
    shell:
        "tar -xzf {input}"
ruleorder: decompress_delta>get_deltadir


rule Dstrain:
    input:
        ANI=lambda wc: expand("mummer/ANI/{species}.tsv",species=get_species(wc)),
        subsets_dir="mummer/subsets",
        delta_dir="mummer/delta"
    output:
        "mummer/delta.tar.gz"
    shell:
        " tar -czf {input.delta_dir}.tar.gz {input.delta_dir} ;"
        "rm -rf {input.subsets_dir} {input.delta_dir}"

rule run_mummer:
    input:
        genome_list="mummer/subsets/{species}.txt",
        genome_folder= genome_folder,
        genome_stats="tables/genome_stats.tsv",
        delta_dir="mummer/delta"
    output:
        "mummer/ANI/{species}.tsv"
    threads:
        config['threads']
    conda:
        "../envs/mummer.yaml"
    resources:
        time= lambda wc, input, threads: estimate_time_mummer(input,threads),
        mem=1
    log:
        "logs/mummer/workflows/{species}.txt"
    params:
        path= os.path.dirname(workflow.snakefile)
    shell:
        "snakemake -s {params.path}/rules/mummer.smk "
        "--config genome_list='{input.genome_list}' "
        " genome_folder='{input.genome_folder}' "
        " species={wildcards.species} "
        " genome_stats={input.genome_stats} "
        " --rerun-incomplete "
        "-j {threads} --nolock 2> {log}"
