from os import path
from itertools import combinations

genomes= open(config['genome_list']).read().strip().replace('.fasta','').split()
genome_folder= config['genome_folder']
species= config['species']
genome_stats= config['genome_stats']

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
        "logs/mummer/mummer/{ref}-{query}.txt"
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
        temp(f"mummer/alignements/{species}.txt")
    shell:
        "cat {input} > {output}"


rule calculate_ANI:
    input:
        alignments=rules.combine.output[0],
        genome_stats= genome_stats
    output:
        f"mummer/alignements/{species}.tsv"
    run:
        import pandas as pd

        M= pd.read_csv(input.alignments,index_col=[0,1],header=None,sep='\t')
        genome_stats= pd.read_csv(input.genome_stats,sep='\t',index_col=0)

        M.columns=['Aln_length','Sim_errors']
        M.index.names=['Genome1','Genome2']
        M['ANI']= 1- M.Sim_errors / M.Aln_length
        M['GenomeSize1']= genome_stats.loc[M.index.get_level_values(0),'Length'].values
        M['GenomeSize2']= genome_stats.loc[M.index.get_level_values(1),'Length'].values
        M['Coverage']= M.Aln_length / M[['GenomeSize1','GenomeSize2']].max(1)

        M.to_csv(output[0],sep='\t')
