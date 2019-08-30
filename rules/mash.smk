

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



localrules: cluster_species,get_representatives
checkpoint cluster_species:
    input:
        dists=rules.mash_calculate_dist.output[0],
        quality ="tables/Genome_quality.tsv",
    output:
        cluster_file="tables/mag2species.tsv",
        scores="tables/evaluation_species_clustering.tsv"
    params:
        treshold=config['species_treshold'],
        linkage_method=config.get('linkage_method','average'),
    script:
        "../scripts/group_species.py"


def get_species(wildcards):
    import pandas as pd
    cluster_file=checkpoints.cluster_species.get().output.cluster_file

    df= pd.read_csv(cluster_file,sep='\t',index_col=0)
    return list(df.Species.unique())

def get_representative_mapping(cluster_file,resolution_level):
    import pandas as pd
    df= pd.read_csv(cluster_file,sep='\t')

    if resolution_level=='species':
        rep= "Representative_Species"
    elif resolution_level=='strains':
        rep= "Representative_Strain"
    else: raise Exception(f"taxrank should be strains or species got {resolution_level} ")

    return df[rep]


def get_representatives(wildcards):
    import pandas as pd

    resolution_level= wildcards.resolution_level

    if resolution_level=='species':

        cluster_file=checkpoints.cluster_species.get().output.cluster_file
    elif resolution_level=='strains':
        cluster_file=checkpoints.cluster_strains.get().output.cluster_file

    mapping= get_representative_mapping(cluster_file, resolution_level)

    return list(mapping.unique())


checkpoint get_representatives:
    input:
        dir= genome_folder,
        cluster_file= "tables/mag2{resolution_level}.tsv"
    output:
        dir= directory("representatives/{resolution_level}"),
    run:
        mapping = get_representative_mapping(input.cluster_file,wildcards.resolution_level)

        output_dir = output.dir
        os.makedirs(output_dir)
        input_dir= os.path.relpath(input.dir,start=output_dir)

        for genome in mapping.unique():
            os.symlink(os.path.join(input_dir,genome+'.fasta'),
                       os.path.join(output_dir,genome+'.fasta')
                   )
