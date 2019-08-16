

def estimate_time_mummer(input,threads):
    "retur time in minutes"

    N= len(open(input.genome_list).read().split())

    time_per_mummer_call = 10/60

    return int((N**2/2*time_per_mummer_call + N)/threads +10)

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
        "-j {threads} --nolock 2> {log}"
