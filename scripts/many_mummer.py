import sys, os
import pandas as pd

sys.stdout = open(snakemake.log[0], "w")
sys.stderr = open(snakemake.log[0], "a")

from snakemake.shell import shell



def run_mummer_and_filter(
    ref,
    query,
    delta_file,
    tmp_folder=snakemake.params.tmpfolder,
    options=snakemake.params.mummer_options,
    threads=snakemake.threads,
    log=snakemake.log[0],
    method= "mum",
    filter_options="-r -q"
):

    tmp_delta_prefix= f"{tmp_folder}/{delta_file.replace('/','_').replace('.delta','')}"

    os.makedirs(os.path.dirname(delta_file), exist_ok=True)

    shell(
        "nucmer --{method} --prefix {tmp_delta_prefix} {options} {ref} {query} 2> {log} ;"
        "delta-filter {filter_options} {tmp_delta_prefix}.delta > {delta_file} 2> {log}"

    )


def many_mummer(alignment_list,
                genome_folder,
                out_folder,
                genome_stats_file,
                extension=".fasta"
                ):

    stats = pd.read_csv(genome_stats_file, index_col=0, sep="\t")

    os.makedirs(out_folder, exist_ok=True)

    delta_files = []

    with open(alignment_list) as f:
        for line in f:

            # get larger genome as ref and smaler as query
            genome_pair = line.strip().split()
            genome_query, genome_ref = (
                stats.loc[genome_pair, "Length"].sort_values().index
            )

            delta_file = os.path.join(out_folder, genome_ref, genome_query + ".delta")

            if not os.path.exists(delta_file):

                # run minimap
                fasta_ref = os.path.join(genome_folder, genome_ref + extension)
                fasta_query = os.path.join(genome_folder, genome_query + extension)
                run_mummer_and_filter(fasta_ref, fasta_query, delta_file)

            delta_files.append(delta_file)

    return delta_files





def parse_delta(filename):
    """
    """
    aln_length, sim_errors = 0, 0
    for line in open(filename, 'rU'):
        fields = line.strip().split()


        if (fields[0] != 'NUCMER') and  (not fields[0].startswith('>') ) and (len(fields) == 7):
        # We only process lines with seven columns:

            aln_length += abs(int(fields[1]) - int(fields[0]))
            sim_errors += int(fields[4])
    return aln_length, sim_errors

def parse_many_delta_files(delta_files,output_table,genome_stats):
    from io import StringIO
    # Create table in memory as string

    table= ""
    for deltaf in delta_files:
        aln_length, sim_errors= parse_delta(deltaf)
        _,_,genome1,genome2=deltaf.replace('.delta','').split('/')

        table+=f"{genome1}\t{genome2}\t{aln_length}\t{sim_errors}\n"

    # create pandas DataFrame
    M= pd.read_csv(StringIO(table),index_col=[0,1],header=None,sep='\t')
    M.columns=['Aln_length','Sim_errors']
    M.index.names=['Genome1','Genome2']
    del table

    #load genome stats
    genome_stats= pd.read_csv(genome_stats,sep='\t',index_col=0)


    # calculate ANI, Coverage add genome Sizes
    M['ANI']= 1- M.Sim_errors / M.Aln_length
    M['GenomeSize1']= genome_stats.loc[M.index.get_level_values(0),'Length'].values
    M['GenomeSize2']= genome_stats.loc[M.index.get_level_values(1),'Length'].values
    M['Coverage']= M.Aln_length / M[['GenomeSize1','GenomeSize2']].min(1)

    #coverage as min genome


    #save table
    M.to_csv(output_table,sep='\t')





delta_files = many_mummer(
    snakemake.input.alignment_list,
    snakemake.input.genome_folder,
    snakemake.params.out_folder,
    snakemake.input.genome_stats,
)

parse_many_delta_files(delta_files, snakemake.output.output_table,snakemake.input.genome_stats)
