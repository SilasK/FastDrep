

rule run_bbsketch:
    input:
        genome_folder= genome_folder,
    output:
        sketch="bbsketch/mags_{NTorAA}.sketch.gz",
        dists= "tables/bbsketch_{NTorAA}.tsv"
    wildcard_constraints:
        NTorAA="(nt|aa)"
    threads:
        config['threads']
    conda:
        "../envs/bbmap.yaml"
    resources:
        mem= 50
    log:
        "logs/bbsketch/workflow_{NTorAA}.log"
    params:
        path= os.path.dirname(workflow.snakefile),
        amino=lambda wildcards: wildcards.NTorAA=='aa'
    shell:
        "snakemake -s {params.path}/rules/bbsketch.smk "
        " --reason "
        "--config  "
        " genome_folder='{input.genome_folder}' "
        " sketch={output.sketch} "
        " amino={params.amino} "
        " dists={output.dists} "
        " --rerun-incomplete "
        " --quiet "
        "-j {threads} --nolock 2> {log}"



rule sendsketch:
    input:
        "bbsketch/mags_aa.sketch.gz"
    output:
        "tables/refseq_mapping.tsv"
    params:
        minid=0.9
    log:
        "logs/bbsketch/sendsketch.log"
    conda:
        "../envs/bbmap.yaml"
    threads:
        1
    resources:
        mem= 1,
        time=10
    shell:
        "sendsketch.sh in={input} out={output} protein format=3 minid={params.minid} 2> {log}"
