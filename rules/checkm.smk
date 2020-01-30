
if not 'database_dir' in config:
    raise Exception("Expect to find a path to the 'database_dir' in the config."
                    "If you haven't yet downloaded the checkm databse run this snakefile with 'download' as target")

input_genome_folder=config['genome_folder']
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

def symlink_relative(files,input_dir,output_dir):
    """create symlink with and adjust for relative path"""


    input_dir_rel= os.path.relpath(input_dir, output_dir)

    for f in files:
        os.symlink(os.path.join(input_dir_rel,f),
                   os.path.join(output_dir,f))

CHECKM_init_flag= 'checkm/init.txt'
#
# rule all:
#      input:
#          checkm="filter/Genome_quality.tsv",
#          markers= "filter/checkm_markers.fasta"


localrules: get_subsets_for_checkm, merge_checkm

checkpoint get_subsets_for_checkm:
    input:
        dir=input_genome_folder
    output:
        dir=directory(temp("checkm/subsets"))
    params:
        subset_size=100
    run:
        files= os.listdir(input.dir)

        N_subsets= len(files)//params.subset_size + 1
        for n in range(N_subsets):
            output_dir= os.path.join(output.dir,f'subset{n:d}')
            os.makedirs(output_dir)
            subset_files=files[n*params.subset_size: min((n+1)*params.subset_size,len(files))]
            symlink_relative(subset_files,input.dir,output_dir)



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





rule run_checkm:
    input:
        touched_output = CHECKM_init_flag,
        bin_dir = "checkm/subsets/{subset}"
    output:
        directory("checkm/results/{subset}"),
    conda:
        "../envs/checkm.yaml"
    threads:
        config.get("threads",8)
    log:
        "checkm/logs/{subset}.txt"
    shell:
        """
        checkm lineage_wf \
            --file {output[0]}/completeness.tsv \
            --tab_table \
            --extension fasta \
            --threads {threads} \
            {input.bin_dir} \
            {output[0]} 2> {log}

            checkm tree_qa \
               {output[0]} \
               --out_format 4 \
               --file {output[0]}/taxonomy.tsv 2>> {log}

        """


def get_subsets_results_dir(wildcards):

    subset_folder= checkpoints.get_subsets_for_checkm.get(**wildcards).output[0]
    subsets= [d for d in os.listdir(subset_folder) if d.startswith('subset')]

    checkmdirs= expand(rules.run_checkm.output[0],
           subset= subsets)
    return checkmdirs


rule merge_checkm:
    input:
        checkmdirs= get_subsets_results_dir
    output:
        checkm="filter/Genome_quality.tsv",
        markers= "filter/checkm_markers.fasta"
    run:

        import pandas as pd
        import shutil
        D=[]
        fout= open(output.markers,'wb')
        for i in range(len(input)):

            completeness= os.path.join(input.checkmdirs[i],'completeness.tsv')
            taxonomy= os.path.join(input.checkmdirs[i],'taxonomy.tsv')
            fasta   = os.path.join(input.checkmdirs[i],"storage/tree/concatenated.fasta")


            df= pd.read_csv(completeness,index_col=0,sep='\t')
            df= df.join(pd.read_csv(taxonomy,index_col=0,sep='\t'))
            D.append(df)

            shutil.copyfileobj(open(fasta,'rb'),fout)

        D= pd.concat(D,axis=0)
        D.to_csv(output.checkm,sep='\t')
        fout.close()




localrules: initialize_checkm
rule initialize_checkm:
    input:
        ancient(f'{DBDIR}/checkm_downloaded')
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
            checkm data setRoot {params.database_dir} &> {log}

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
