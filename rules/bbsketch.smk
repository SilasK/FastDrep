from os import path

genome_folder= config['genome_folder']
assert (path.isdir(genome_folder) & path.exists(genome_folder))
genomes= glob_wildcards(os.path.join(genome_folder,'{genome}.fasta')).genome

sketch=config['sketch']
dists= config['dists']
amino= config['amino']

if amino:
    k="9,12"
    postfix='_aa'
else:
    k="31,24"
    postfix='_nt'

rule all:
    input:
        sketch,
        dists

rule bbsketch:
    input:
        input=os.path.join(genome_folder,"{genome}.fasta")
    output:
        out=f"bbsketch/sketches{postfix}/{{genome}}.sketch"
    params:
        k=k,
        translate=amino,
        overwrite=True,
        command="bbsketch.sh",
        name0="{genome}"
    resources:
        mem= 1
    log:
        f"logs/bbsketch/sketch{postfix}/{{genome}}.log"
    threads:
        1
    script:
        "../scripts/runBB.py"

rule mergesketch:
    input:
        expand(rules.bbsketch.output[0],
                                        genome= genomes)
    output:
        out=sketch
    threads:
        1
    run:
        import gzip
        import shutil
        with gzip.open(output[0], 'wb') as f_out:
            for f in input:
                with open(f, 'rb') as f_in:
                    shutil.copyfileobj(f_in, f_out)




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
    resources:
        mem= 100
    log:
        f"logs/bbsketch/alltoall.log"
    threads:
        16
    script:
        "../scripts/runBB.py"
