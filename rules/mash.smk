

rule mash_sketch_genome:
    input:
        filter_genome_folder,
        genome_folder
    output:
        "genomes.msh"
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
        "mash_dists.txt"
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



localrules: group_species
checkpoint group_species:
    input:
        mash_dists=rules.mash_calculate_dist.output[0],
    output:
        cluster_file="mash/clusters.tsv",
        subsets_dir= directory("mash/clusters")
    params:
        threshold = 1- config['mash']['dist_treshold'],
        fillna=0.8,
        linkage_method='average',
        square=False
    run:
        M= gd.load_mash(input.mash_dists)
        labels= gd.group_species_linkage(M,**params)
        labels.to_csv(output[0],sep='\t',header=False)

        os.makedirs(output.subsets_dir)

        for i in labels.unique():

            genomes_of_cluster= labels.index[labels==i].values

            with open(f"{output.subsets_dir}/species_{i}.txt","w") as f:

                f.write(''.join([g+'.fasta\n' for g in genomes_of_cluster ]))


def get_species_numbers(wildcards):

    dir=checkpoint.group_species.get().output.subsets_dir

    return glob_wildcards(f"{dir}/species_{{i}}.txt").i
