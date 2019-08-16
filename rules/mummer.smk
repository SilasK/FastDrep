from os import path
from itertools import combinations

genomes= open(config['genome_list']).read().strip().replace('.fasta','').split()
genome_folder= config['genome_folder']
species= config['species']

assert (path.isdir(genome_folder) & path.exists(genome_folder))


rule all:
    input:
        f"mummer/alignements/{species}.tsv"


rule run_mummer:
    input:
        ref=path.join(genome_folder,"{ref}.fasta"),
        query=path.join(genome_folder,"{query}.fasta")
    output:
        "mummer/delta/{ref}-{query}.delta"
    params:
        out_prefix= "mummer/delta/{ref}-{query}",
        options="--mincluster 65 --maxgap 90 ",
        method= "mum"
    log:
        "logs/mummer/{ref}-{query}.txt"
    threads:
        1
    shell:
        "nucmer --{params.method} --prefix {params.out_prefix} {params.options} {input.ref} {input.query} 2> {log}"

rule delta_filter:
    input:
        rules.run_mummer.output[0]
    output:
        temp("mummer/delta/{ref}-{query}.delta.filtered")
    params:
        options="-r -q"
    shell:
        "delta-filter {params.options} {input} > {output}"


def parse_delta(filename):
    """
    """
    aln_length, sim_errors = 0, 0
    for line in open(filename, 'rU'):
        fields = line.strip().split()


        if (fields[0] != 'NUCMER') and  (not fields[0].startswith('>') ) and (len(fields) == 7):
        # We only process lines with seven columns:

            aln_length += abs(int(fields[1]) - int(fields[0]))
            sim_errors += int(fields[4])
    return aln_length, sim_errors

rule parse_delta:
    input:
        rules.delta_filter.output[0]
    output:
        temp("mummer/delta/{ref}-{query}.txt")
    run:

        aln_length, sim_errors = parse_delta(input[0])
        with open(output[0],'w') as f:
            f.write(f"{wildcards.ref}\t{wildcards.query}\t{aln_length}\t{sim_errors}\n")



rule combine:
    input:
        ["mummer/delta/{}-{}.txt".format(*pair) for pair in combinations(genomes,2)]
    output:
        f"mummer/alignements/{species}.tsv"
    shell:
        "cat {input} > {output}"
