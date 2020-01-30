



ANI_file= f"tables/{config['sketcher']}_dists.tsv"

localrules: cluster_species,get_representatives
checkpoint cluster_species:
    input:
        dists=ANI_file,
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




localrules: cluster_strains
checkpoint cluster_strains:
    input:
        dists=ANI_file,
        quality ="tables/Genome_quality.tsv",
        mag2species= "tables/mag2species.tsv"
    output:
        mag2strain="tables/mag2strains.tsv"
    params:
        #treshold=config['species_treshold'],
        linkage_method=config.get('linkage_method','average'),
    script:
        "../scripts/group_strains.py"


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



checkpoint get_ref_bbsplit:
    input:
        dir= genome_folder,
        cluster_file= "tables/mag2strains.tsv"
    output:
        dir= directory("representatives/merged_strains"),
    run:

        import pandas as pd
        df= pd.read_csv(input.cluster_file,sep='\t')


        os.makedirs(output.dir)

        for species_name,sp in df.groupby('Species'):
            with open(os.path.join(output.dir,species_name+'.fasta'),'w') as fasta_out:
                for genome in sp.Representative_Strains.unique():

                    with open(os.path.join(input.dir,genome+'.fasta')) as fasta_in:
                        fasta_out.write(fasta_in.read())
