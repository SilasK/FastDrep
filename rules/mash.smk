

rule mash_sketch_genome:
    input:
        genome_folder
    output:
        "mash/genomes.msh"
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
        "mash sketch -o {params.out_name} -p {threads} -s {params.s} -k {params.k} {input[0]}/* 2> {log}"


rule mash_calculate_dist:
    input:
        genomes=rules.mash_sketch_genome.output
    output:
        "mash/dists.txt"
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



localrules: cluster_mash
rule cluster_mash:
    input:
        mash_dists=rules.mash_calculate_dist.output[0],
    output:
        cluster_file="mash/clusters.tsv",
    params:
        threshold = 1- config['mash']['dist_treshold'],
        fillna=0.8,
        linkage_method='average',
        square=False
    run:
        M= gd.load_mash(input.mash_dists)
        labels= gd.group_species_linkage(M,**params)

        if min(labels)==0: labels+=1
        n_leading_zeros= len(str(max(labels)))
        format_int='Species{:0'+str(n_leading_zeros)+'d}'

        labels=labels.apply(format_int.format)

        labels.to_csv(output[0],sep='\t',header=False)
