from os import path

genome_folder= config['genome_folder']
assert (path.isdir(genome_folder) & path.exists(genome_folder))
genomes= glob_wildcards(os.path.join(genome_folder,'{genome}.fasta')).genome

sketch=config['sketch']
dists= config['dists']
amino= config['amino']

if amino:
    k="9,12"
else:
    k=31

rule all:
    input:
        sketch,
        dists

rule bbsketch:
    input:
        input=os.path.join(genome_folder,"{genome}.fasta")
    output:
        out=temp("bbsketch/sketches/{genome}.sketch")
    params:
        k=k,
        translate=amino,
        overwrite=True,
        command="bbsketch.sh",
        name0="{genome}"
    threads:
        1
    script:
        "../scripts/runBB.py"

rule mergesketch:
    input:
        expand("bbsketch/sketches/{genome}.sketch",
                                        genome= genomes)
    output:
        out=sketch
    threads:
        1
    shell:
        "cat {input} | gzip > {output}"



rule allvall:
    input:
        ref=sketch
    output:
        out=dists
    params:
        amino=amino,
        overwrite=True,
        command="comparesketch.sh alltoall",
        format=3
    threads:
        16
    script:
        "../scripts/runBB.py"
