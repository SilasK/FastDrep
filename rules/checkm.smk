
if not 'database_dir' in config:
    raise Exception("Expect to find a path to the 'database_dir' in the config."
                    "If you haven't yet downloaded the checkm databse run this snakefile with 'download' as target")

DBDIR = os.path.realpath(config["database_dir"])
CHECKMDIR = os.path.join(DBDIR, "checkm")
CHECKM_ARCHIVE="checkm_data_2015_01_16.tar.gz"
CHECKM_ARCHIVE_HASH= "631012fa598c43fdeb88c619ad282c4d"
CHECKM_ADRESS= "https://data.ace.uq.edu.au/public/CheckM_databases"


import os
import hashlib
def md5(fname):
    # https://stackoverflow.com/questions/3431825/generating-an-md5-checksum-of-a-file
    hash_md5 = hashlib.md5()
    if not os.path.exists(fname):
        return None
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()

genome_folder=config['genome_folder']
checkm_dir= config.get('checkm_dir','checkm')
CHECKM_init_flag= os.path.join(checkm_dir,'init.txt')

rule all:
    input:
        expand("{checkm_dir}/checkm/{file}",checkm_dir=checkm_dir,
               file=['taxonomy.tsv',"completeness.tsv","storage/tree/concatenated.fasta"])



checkpoint get_subsets_for_checkm:
    input:
        dir=genome_folder
    output:
        dir=directory(temp("{checkm_dir}/subsets"))
    params:
        subset_size=500
    run:

        for i,f in enumerate():

        files= os.listdir(input.dir)
        n=1



            n= i % params.subset_size

            output_dir= os.path.join(output.dir,f'subset{n:d}')
            os.path.makedirs(output_dir)
            def symlink(files,input_dir,output_dir):
                """create symlink with and adjust for relative path"""


            input_dir_rel= os.path.relpath(input.dir, output_dir)

            for f in files:
                os.symlink(os.path.join(input_dir_rel,f),
                           os.path.join(output_dir,f))


rule run_checkm:
    shell:


rule combine_checkm_quality:

    output:
        "filter/Genome_quality.tsv"







rule download:
    output:
        touch(f'{DBDIR}/checkm_downloaded'),
        tarbal=temp(f"{DBDIR}/{CHECKM_ARCHIVE}"),
        dir = CHECKMDIR
    threads:
        1
    run:
        shell("wget -O {output.tarbal} '{CHECKM_ADRESS}/{CHECKM_ARCHIVE}' ")
        if not CHECKM_ARCHIVE_HASH == md5(output.tarbal):
            raise OSError(2, "Invalid checksum", output.tarbal)
        os.makedirs(output.dir)
        shell("tar -zxf {output.tarbal} --directory {output.dir}")





rule run_checkm_lineage_wf:
    input:
        touched_output = CHECKM_init_flag,
        bin_dir = genome_folder
    output:
        "{checkm_dir}/completeness.tsv",
        "{checkm_dir}/storage/tree/concatenated.fasta"
    conda:
        "../envs/checkm.yaml"
    threads:
        config.get("threads",8)
    log:
        "{checkm_dir}/logs/lineage_wf.txt"
    shell:
        """
        rm -r {wildcards.checkm_dir} ;
        checkm lineage_wf \
            --file {output[0]} \
            --tab_table \
            --quiet \
            --extension fasta \
            --threads {threads} \
            {input.bin_dir} \
            {wildcards.checkm_dir} 2> {log}
        """



rule run_checkm_tree_qa:
    input:
        tree="{checkm_dir}/completeness.tsv"
    output:
        "{checkm_dir}/taxonomy.tsv",
    conda:
        "../envs/checkm.yaml"
    threads:
        1
    log:
        "{checkm_dir}/logs/tree_qa.txt"
    shell:
        """
            checkm tree_qa \
               {wildcards.checkm_dir} \
               --out_format 4 \
               --file {output[0]} 2> {log}

        """





localrules: initialize_checkm
rule initialize_checkm:
    input:
        ancient(rules.download.output[0])
    output:
        touched_output = touch(CHECKM_init_flag)
    params:
        database_dir = CHECKMDIR,
    conda:
        "../envs/checkm.yaml"
    log:
        "logs/initialize_checkm.log"
    shell:
        """
            checkm data setRoot {params.database_dir} 2> {log}

        """



#
# rule checkm_tetra:
#     input:
#         contigs=BINNING_CONTIGS,
#     output:
#         "checkm/tetranucleotides.txt"
#     log:
#         "checkm/tetra.txt"
#     conda:
#         "../envs/checkm.yaml"
#     threads:
#         config.get("threads", 8)
#     shell:
#         """
#             checkm tetra \
#             --threads {threads} \
#             {input.contigs} {output} 2> {log}
#         """
#
#
# rule checkm_outliers:
#     input:
#         tetra= "checkm/tetranucleotides.txt",
#         bin_folder= bin_folder,
#         checkm = "checkm/completeness.tsv"
#     params:
#         checkm_folder = lambda wc, input: os.path.dirname(input.checkm),
#         report_type = 'any',
#         treshold = 95 #reference distribution used to identify outliers; integer between 0 and 100 (default: 95)
#     output:
#         "checkm/outliers.txt"
#     log:
#         "logs/checkm/outliers.txt"
#     conda:
#         "../envs/checkm.yaml"
#     threads:
#         config.get("threads", 8)
#     shell:
#         """
#             checkm outliers \
#             --extension fasta \
#             --distributions {params.treshold} \
#             --report_type {params.report_type} \
#             {params.checkm_folder} \
#             {input.bin_folder} \
#             {input.tetra} \
#             {output} 2> {log}
#         """
#
#
# rule refine_bins:
#     input:
#         expand("checkm/outliers.txt",sample=SAMPLES,binner=config['binner'])
#
# rule find_16S:
#     input:
#         contigs=BINNING_CONTIGS,
#         bin_dir= bin_folder
#     output:
#         summary="{sample}/binning/{binner}/SSU/ssu_summary.tsv",
#         fasta="{sample}/binning/{binner}/SSU/ssu.fna",
#     params:
#         output_dir = lambda wc, output: os.path.dirname(output[0]),
#         evalue = 1e-05,
#         concatenate = 200 #concatenate hits that are within the specified number of base pairs
#     conda:
#         "../envs/checkm.yaml"
#     threads:
#         1
#     shell:
#         """
#         rm -r {params.output_dir} && \
#            checkm ssu_finder \
#                --extension fasta \
#                --threads {threads} \
#                --evalue {params.evalue} \
#                --concatenate {params.concatenate} \
#                {input.contigs} \
#                {input.bin_dir} \
#                {params.output_dir}
#         """
#
#
# rule get_all_16S:
#     input:
#         summaries= expand(rules.find_16S.output.summary,sample=SAMPLES,binner=config['final_binner']),
#         fastas= expand(rules.find_16S.output.fasta,sample=SAMPLES,binner=config['final_binner'])
#     output:
#         fasta="genomes/SSU/ssu.fasta",
#         summary ="genomes/SSU/ssu_summary.tsv"
#     run:
#         shell("cat {input.fastas} > {output.fasta}")
#
#         import pandas as pd
#         summary= pd.DataFrame()
#
#         for file in input.summaries:
#             try:
#                 d = pd.read_csv(file,index_col=0,sep='\t')
#                 summary=summary.append(d)
#             except:
#                 pd.errors.EmptyDataError
#
#         summary.to_csv(output.summary,sep='\t')
