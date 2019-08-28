
rule Dstrain:
    input:
        ANI="tables/dist_strains.tsv",
        ani_dir='mummer/ANI',
        delta_dir="mummer/delta"
    output:
        "mummer/delta.tar.gz"
    shell:
        " tar -czf {input.delta_dir}.tar.gz {input.delta_dir} ;"
        "rm -rf {input.delta_dir} {input.ani_dir}"


localrules: species_subsets
checkpoint species_subsets:
    input:
        cluster_file=rules.cluster_mash.output.cluster_file,
    output:
        subsets_dir= temp(directory("mummer/subsets"))
    run:
        import pandas as pd
        labels= pd.read_csv(input[0],sep='\t',index_col=0).Species


        os.makedirs(output.subsets_dir)

        for species in labels.unique():

            genomes_of_cluster= labels.index[labels==species].values
            if len(genomes_of_cluster)>1:

                with open(f"{output.subsets_dir}/{species}.txt","w") as f:

                    f.write(''.join([g+'.fasta\n' for g in genomes_of_cluster ]))


def get_species_for_sub_clustering(wildcards):
    import pandas as pd
    subset_dir= checkpoints.species_subsets.get() # subsetdir must be present
    cluster_file=checkpoints.cluster_mash.get().output.cluster_file

    df= pd.read_csv(cluster_file,sep='\t',index_col=0)

    Nspecies= df.groupby('Species').size()


    return list(Nspecies.index[Nspecies>1])



def estimate_time_mummer(input,threads):
    "retur time in minutes"

    N= len(open(input.genome_list).read().split())

    time_per_mummer_call = 10/60

    return int(N**2/2*time_per_mummer_call + N/2)//threads + 5

localrules: get_deltadir,decompress_delta,Dstrain
rule get_deltadir:
    output:
        directory("mummer/delta")
    run:
        os.makedirs(output[0])

rule decompress_delta:
    input:
        "mummer/delta.tar.gz"
    output:
        directory("mummer/delta")
    shell:
        "tar -xzf {input}"
ruleorder: decompress_delta>get_deltadir


rule merge_mummer_ani:
    input:
        lambda wc: expand("mummer/ANI/{species}.tsv",species=get_species_for_sub_clustering(wc))
    output:
        "tables/dist_strains.tsv"
    run:
        import pandas as pd
        Mummer={}
        for file in input:
            Mummer[io.simplify_path(file)]= pd.read_csv(file,index_col=[0,1],sep='\t')

        M= pd.concat(Mummer,axis=0)
        M['Species']=M.index.get_level_values(0)
        M.index= M.index.droplevel(0)
        M.to_csv(output[0],sep='\t')

        #sns.jointplot('ANI','Coverage',data=M.query('ANI>0.98'),kind='hex',gridsize=100,vmax=200)




rule run_mummer:
    input:
        genome_list="mummer/subsets/{species}.txt",
        genome_folder= genome_folder,
        genome_stats="tables/genome_stats.tsv",
        delta_dir="mummer/delta"
    output:
        temp("mummer/ANI/{species}.tsv")
    threads:
        config['threads']
    conda:
        "../envs/mummer.yaml"
    resources:
        time= lambda wc, input, threads: estimate_time_mummer(input,threads),
        mem= 2
    log:
        "logs/mummer/workflows/{species}.txt"
    params:
        path= os.path.dirname(workflow.snakefile)
    shell:
        "snakemake -s {params.path}/rules/mummer.smk "
        "--config genome_list='{input.genome_list}' "
        " genome_folder='{input.genome_folder}' "
        " species={wildcards.species} "
        " genome_stats={input.genome_stats} "
        " --rerun-incomplete "
        "-j {threads} --nolock 2> {log}"


localrules: cluster_strains
checkpoint cluster_strains:
    input:
        dists=rules.merge_mummer_ani.output[0],
        quality ="tables/Genome_quality.tsv",
        mag2species= "tables/mag2species.tsv"
    output:
        cluster_file="tables/mag2strains.tsv"
    params:
        #treshold=config['species_treshold'],
        linkage_method=config.get('linkage_method','average'),
    script:
        "../scripts/group_strains.py"


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
