localrules: mash_list
rule mash_list:
    input:
        genome_folder
    output:
        temp("sketches/genome_list.txt")
    run:
        from glob import glob
        with open(output[0],'w') as f:
            f.write('\n'.join(glob(f'{input[0]}/*.fasta'))+'\n' )


rule mash_sketch_genome:
    input:
        rules.mash_list.output
    output:
        "sketches/genomes.msh"
    params:
        k= config['mash']['k'],
        s= config['mash']['sketch_size'],
        out_name = lambda wc, output: os.path.splitext(output[0])[0]
    threads:
        config['threads']
    conda:
        "../envs/mash.yaml"
    log:
        "logs/mash/sketch.log"
    shell:
        "mash sketch -o {params.out_name} -p {threads} -s {params.s} -k {params.k} -l {input[0]} 2> {log}"


rule mash_calculate_dist:
    input:
        genomes=rules.mash_sketch_genome.output
    output:
        "tables/mash_dists.txt"
    params:
        d= config['mash']['max_dist']
    threads:
        config['threads']
    conda:
        "../envs/mash.yaml"
    log:
        "logs/mash/dist.log"
    shell:
        "mash dist -p {threads} -d {params.d} {input.genomes} {input.genomes} > {output[0]} 2> {log}"


rule bindash_sketch_genome:
    input:
        rules.mash_list.output
    output:
        "sketches/genomes.bdsh",
        "sketches/genomes.bdsh.dat",
        "sketches/genomes.bdsh.txt"
    params:
        k= config['bindash']['k'],
        sketchsize64= config['bindash']['sketch_size'],
        extra=config['bindash'].get('extra',"")
    threads:
        config['threads']
    conda:
        "../envs/bindash.yaml"
    log:
        "logs/bindash/sketch.log"
    shell:
        "bindash sketch "
        "--outfname={output[0]} "
        "--nthreads={threads} "
        "--sketchsize64={params.sketchsize64} "
        "--kmerlen={params.k} "
        "{params.extra} "
        "--listfname={input[0]} 2> {log}"

rule bindash_dist:
    input:
        rules.bindash_sketch_genome.output[0]
    output:
        "tables/bindash_dists.tsv"
    params:
        d= config['bindash']['max_dist']
    threads:
        config['threads']
    conda:
        "../envs/mash.yaml"
    log:
        "logs/bindash/dist.log"
    shell:
        "bindash dist "
        "--nthreads={threads} "
        "--mthres={params.d} "
        "--outfname={output} {input[0]} 2> {log}"


def get_minhash()

localrules: filter_minhash
checkpoint filter_minhash:
    input:
        rules.bindash_calculate_dist.output[0] if (config['sketcher']=='bindash') else rules.mash_calculate_dist.output[0]
    output:
        temp(directory('mummer/subsets'))
    params:
        treshold=config['mummer']['max_dist'],
        N=config['mummer']['subset_size']
    run:

        F= gd.load_mash(input[0])
        G= gd.to_graph(F.query(f"Distance<={params.treshold}"))
        G.remove_edges_from(G.selfloop_edges())

        os.makedirs(output[0])

        fout=None
        for i,e in enumerate(G.edges()):
            if (i % params.N) ==0:
                n= int(i // params.N )+1
                if fout is not None: fout.close()
                fout= open(f"{output[0]}/subset_{n}.txt",'w')

            fout.write("\t".join(sorted(e))+'\n')

def get_mummer_subsets(wildcards):
    subset_dir= checkpoints.filter_minhash.get().output[0]

    subsets= glob_wildcards(os.path.join(subset_dir,'{subset}.txt')).subset

    return subsets
