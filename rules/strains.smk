


def estimate_time_mummer(N,threads):
    "retur time in minutes"

    time_per_mummer_call = 1 # min

    return int(N*time_per_mummer_call)//threads + 5



rule run_mummer:
    input:
        comparison_list="mummer/subsets/{subset}.txt",
        genome_folder= genome_folder,
        genome_stats="tables/genome_stats.tsv",
    output:
        temp("mummer/ANI/{subset}.tsv")
    threads:
        config['threads']
    conda:
        "../envs/mummer.yaml"
    resources:
        time= lambda wc, input, threads: estimate_time_mummer(config['mummer']['subset_size'],threads),
        mem= 2
    log:
        "logs/mummer/workflows/{subset}.txt"
    params:
        path= os.path.dirname(workflow.snakefile)
    shell:
        "snakemake -s {params.path}/rules/mummer.smk "
        "--config comparison_list='{input.comparison_list}' "
        " genome_folder='{input.genome_folder}' "
        " subset={wildcards.subset} "
        " genome_stats={input.genome_stats} "
        " --rerun-incomplete "
        "-j {threads} --nolock 2> {log}"



def get_merge_mummer_ani_input(wildcards):

    subsets=get_mummer_subsets(wildcards)

    return expand("mummer/ANI/{subset}.tsv",subset=subsets)

localrules: merge_mummer_ani
rule merge_mummer_ani:
    input:
        get_merge_mummer_ani_input
    output:
        "tables/mummer_dist.tsv"
    run:
        import pandas as pd
        import shutil

        Mummer={}
        for file in input:
            Mummer[io.simplify_path(file)]= pd.read_csv(file,index_col=[0,1],sep='\t')

        M= pd.concat(Mummer,axis=0)
        M.index= M.index.droplevel(0)
        M.to_csv(output[0],sep='\t')


        ani_dir= os.path.dirname(input[0])
        shutil.rmtree(ani_dir)
        #sns.jointplot('ANI','Coverage',data=M.query('ANI>0.98'),kind='hex',gridsize=100,vmax=200)



localrules: cluster_species,get_representatives
checkpoint cluster_species:
    input:
        dists=rules.merge_mummer_ani.output[0],
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
        dists=rules.merge_mummer_ani.output[0],
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
