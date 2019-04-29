import os
from glob import glob

#genomes=glob("/Users/silas/Atlas/bins/DASTool3/genomes/*.fasta") + \
#glob("/Users/silas/PhD/06_Projects/Warm_claire_microbiota/Metagenome/WD/genomes/genomes/*.fasta")



localrules: gen_genome_list
rule gen_genome_list:
    input:
        genomes=genomes
    output:
        genome_list=temp('fastANI_genome_list.txt')
    run:
        with open(output.genome_list,'w') as fl:
            for g in input.genomes:
                fl.write(f'{g}\n')



rule fastANI:
    input:
        genome_list= rules.gen_genome_list.output[0],
        genomes= genomes
    output:
        "ANIs.txt",
    resources:
        mem=30
    threads:
        config['threads']
    conda:
        "envs/fastANI.yaml"
    log:
        "logs/fastANI.log"
    benchmark:
        "benchmarks/fastANI.txt"
    shell:
        " fastANI "
        " --threads {threads} "
        " --queryList {input.genome_list} --refList {input.genome_list} "
        " -k {config[fastani][k]} "
        " --fragLen {config[fastani][fragLen]} "
        " --minFrag {config[fastani][minFrag]} "
        " -o {output} "
        " 2> {log}"
