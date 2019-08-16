

def estimate_time_mummer(genome_list):

    N= len(open(genome_list).read().split())

    (N**2/2)*0.33 + N +10

rule Dstrain:
    input:
        lambda wc: expand("mummer/alignements/species_{i}.tsv",i=get_species_numbers(wc))

rule run_mummer:
    input:
        genome_list="mash/clusters/{species}.txt",
        genome_folder= genome_folder,
        genome_stats="genome_stats.tsv"
    output:
        "mummer/alignements/{species}.tsv"
    threads:
        8
    conda:
        "../envs/mummer.yaml"

    params:
        path= os.path.dirname(workflow.snakefile)
    shell:
        "snakemake -s {params.path}/rules/mummer.smk "
        "--config genome_list='{input.genome_list}' "
        " genome_folder='{input.genome_folder}' "
        " species={wildcards.species} "
        " genome_stats={input.genome_stats} "
        "-j {threads} "
