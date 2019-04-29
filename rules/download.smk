# The main entry point of your workflow.
# After configuring, running snakemake -n in a clone of this repository should successfully execute a dry-run of the workflow.


#configfile: "config.yaml"


rule all:
    input:
        "ftp_list.txt"



import pandas as pd

acc_table= pd.read_csv(config['acc_table'],sep='\t',index_col=0)
accessions= acc_table.index[:3]

from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
FTP = FTPRemoteProvider()

rule ftp_summary:
    input:
        # only keeping the file so we can move it out to the cwd
        FTP.remote("ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt", keep_local=True)
    output:
        "assembly_summary.txt"
    run:
        shell("mv {input} ./")

rule get_fasta:
    input:
        summary="assembly_summary.txt"
    output:
        "ftp_list.txt"
    run:
        import pandas as pd
        D=pd.read_csv(input.summary,sep='\t',skiprows=0)
        D.loc[accessions].to_csv(output[0],sep='\t')


download_genome:
    input:
        "ftp_list.txt"
    output:
        expand("genomes",accession=accessions)


include: "rules/other.smk"
