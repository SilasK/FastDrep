DBDIR = os.path.realpath(config["database_dir"])
CHECKMDIR = os.path.join(DBDIR, "checkm")
CHECKM_ARCHIVE = "checkm_data_v1.0.9.tar.gz"
CONDAENV='envs'
CHECKM_init_flag= os.path.join(CHECKM_ARCHIVE,'init.txt')

import os
from glob import glob



genome_folder=config['genome_folder']
genome_fasta= glob(os.path.join(genome_folder,"*"))
genomes = [os.path.splitext(g)[0] for g in genome_fasta]

genome_fasta= dict(zip(genomes, genome_fasta))
def get_genome_fasta(wildcards):
    return genome_fasta[wildcards.genome]

rule all:
    input:
        expand("{genome}/checkm/{file}",genome=genomes,
               file=['taxonomy.tsv',"completeness.tsv","storage/tree/concatenated.fasta"])

rule init_checkm:
    input:
        get_genome_fasta
    output:
        "{genome}/{genome}.fasta"
    shell:
        "cp {input} {output}"


rule run_checkm_lineage_wf:
    input:
        touched_output = CHECKM_init_flag,
        bins = ["{genome}/{genome}.fasta"] # actualy path to fastas
    output:
        "{genome}/checkm/completeness.tsv",
        "{genome}/checkm/storage/tree/concatenated.fasta"
    params:
        output_dir = lambda wc, output: os.path.dirname(output[0]),
        bin_dir= lambda wc, input: os.path.dirname(input.bins[0]),
    conda:
        "%s/checkm.yaml" % CONDAENV
    threads:
        config.get("threads", 1)
    shell:
        """
        rm -r {params.output_dir}
        checkm lineage_wf \
            --file {params.output_dir}/completeness.tsv \
            --tab_table \
            --quiet \
            --extension fasta \
            --threads {threads} \
            {prams.bin_dir} \
            {params.output_dir}
        """



rule run_checkm_tree_qa:
    input:
        tree="{checkmfolder}/completeness.tsv"
    output:
        summary="{checkmfolder}/taxonomy.tsv",
    params:
        tree_dir = lambda wc, input: os.path.dirname(input.tree),
    conda:
        "%s/checkm.yaml"  % CONDAENV
    threads:
        1
    shell:
        """
            checkm tree_qa \
               {params.tree_dir} \
               --out_format 4 \
               --file {output.netwick}

        """


CHECKMFILES=[   "%s/taxon_marker_sets.tsv" % CHECKMDIR,
        "%s/selected_marker_sets.tsv" % CHECKMDIR,
        "%s/pfam/tigrfam2pfam.tsv" % CHECKMDIR,
        "%s/pfam/Pfam-A.hmm.dat" % CHECKMDIR,
        "%s/img/img_metadata.tsv" % CHECKMDIR,
        "%s/hmms_ssu/SSU_euk.hmm" % CHECKMDIR,
        "%s/hmms_ssu/SSU_bacteria.hmm" % CHECKMDIR,
        "%s/hmms_ssu/SSU_archaea.hmm" % CHECKMDIR,
        "%s/hmms_ssu/createHMMs.py" % CHECKMDIR,
        "%s/hmms/phylo.hmm.ssi" % CHECKMDIR,
        "%s/hmms/phylo.hmm" % CHECKMDIR,
        "%s/hmms/checkm.hmm.ssi" % CHECKMDIR,
        "%s/hmms/checkm.hmm" % CHECKMDIR,
        "%s/genome_tree/missing_duplicate_genes_97.tsv" % CHECKMDIR,
        "%s/genome_tree/missing_duplicate_genes_50.tsv" % CHECKMDIR,
        "%s/genome_tree/genome_tree.taxonomy.tsv" % CHECKMDIR,
        "%s/genome_tree/genome_tree_reduced.refpkg/phylo_modelJqWx6_.json" % CHECKMDIR,
        "%s/genome_tree/genome_tree_reduced.refpkg/genome_tree.tre" % CHECKMDIR,
        "%s/genome_tree/genome_tree_reduced.refpkg/genome_tree.log" % CHECKMDIR,
        "%s/genome_tree/genome_tree_reduced.refpkg/genome_tree.fasta" % CHECKMDIR,
        "%s/genome_tree/genome_tree_reduced.refpkg/CONTENTS.json" % CHECKMDIR,
        "%s/genome_tree/genome_tree.metadata.tsv" % CHECKMDIR,
        "%s/genome_tree/genome_tree_full.refpkg/phylo_modelEcOyPk.json" % CHECKMDIR,
        "%s/genome_tree/genome_tree_full.refpkg/genome_tree.tre" % CHECKMDIR,
        "%s/genome_tree/genome_tree_full.refpkg/genome_tree.log" % CHECKMDIR,
        "%s/genome_tree/genome_tree_full.refpkg/genome_tree.fasta" % CHECKMDIR,
        "%s/genome_tree/genome_tree_full.refpkg/CONTENTS.json" % CHECKMDIR,
        "%s/genome_tree/genome_tree.derep.txt" % CHECKMDIR,
        "%s/.dmanifest" % CHECKMDIR,
        "%s/distributions/td_dist.txt" % CHECKMDIR,
        "%s/distributions/gc_dist.txt" % CHECKMDIR,
        "%s/distributions/cd_dist.txt" % CHECKMDIR

localrules: initialize_checkm
rule initialize_checkm:
    input:
        ancient(CHECKMFILES)
    output:
        touched_output = CHECKM_init_flag
    params:
        database_dir = CHECKMDIR,
        script_dir = os.path.dirname(os.path.abspath(workflow.snakefile))
    conda:
        "%s/checkm.yaml" % CONDAENV
    log:
        "logs/initialize_checkm.log"
    shell:
        """
        python {params.script_dir}/rules/initialize_checkm.py \
            {params.database_dir} \
            {output.touched_output} \
            {log}
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
#         "%s/checkm.yaml" % CONDAENV
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
#         "%s/checkm.yaml" % CONDAENV
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
#         "%s/checkm.yaml" % CONDAENV
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
