
rule sendsketch:
    input:
        os.path.join(genome_folder,"{genome}.fasta")
    output:
        "taxonomy/sendsketch/{genome}.tsv"
    conda:
        "../envs/bbmap.yaml"
    threads:
        1
    resources:
        mem= 1,
        time=10
    shell:
        "sendsketch.sh in={input} out={output} protein -Xmx{resources.mem}g"

rule combine_tax_sketch:
    input:
        lambda wildcards: expand(rules.sendsketch.output[0], genome= get_representative_species(wildcards))
