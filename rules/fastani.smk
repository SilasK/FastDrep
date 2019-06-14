import os


localrules: gen_genome_lists,combine_fastANI
checkpoint gen_genome_lists:
    input:
        filter_genome_folder
    output:
        clusters=directory("clusters/fastani")
    params:
        cluster_size=config['fastani']['subset_size']
    run:
        from glob import glob
        genomes= glob(os.path.join(input[0],"*.fasta"))

        os.makedirs(output[0])

        c, file=0, None
        for i,g in enumerate(genomes):
            if (i % params.cluster_size) == 0:
                c+=1

                if file is not None: file.close()
                file = open(os.path.join(output.clusters,f"Cluster{c}.txt"),'w')

            file.write(f"{g}\n")


        file.close()


rule fastANI:
    input:
        list1= "clusters/fastani/Cluster{i}.txt",
        list2= "clusters/fastani/Cluster{j}.txt",
        genomes= filter_genome_folder,
        all_Genomes= genome_folder
    output:
        "clusters/fastani/ANIcluster{i}-{j}.txt",
    resources:
        mem=30
    threads:
        config['threads']
    conda:
        "../envs/fastANI.yaml"
    log:
        "logs/fastANI/cluster{i}-{j}.log"
    benchmark:
        "logs/benchmarks/fastANI/cluster{i}-{j}.log"
    shell:
        " fastANI "
        " --threads {threads} "
        " --queryList {input.list1} --refList {input.list2} "
        " -k {config[fastani][k]} "
        " --fragLen {config[fastani][fragLen]} "
        " --minFrag {config[fastani][minFrag]} "
        " -o {output} "
        " 2> {log}"


def get_allFastANI(wildcards):
    ALL_I=  glob_wildcards(os.path.join(checkpoints.gen_genome_lists.get(**wildcards).output.clusters,
                                        "Cluster{i}.txt")).i

    return expand("clusters/fastani/ANIcluster{i}-{j}.txt", i=ALL_I,j=ALL_I )


rule combine_fastANI:
    input:
        get_allFastANI
    output:
        "ANI.tsv"
    shell:
        "cat {input} > {output};"
        "rm -r clusters "
