
import sys,os


sys.stdout = open(snakemake.log[0], 'w')

from snakemake.shell import shell



def run_mummer(ref,query,out_prefix,options="-c 65 -g 90 ",#snakemake.params,
                threads=snakemake.threads,log=snakemake.log[0]):
    shell("nucmer --{method} -p {out_prefix} {options} {ref} {query} 2>> {log}")

    shell("delta-filter -r -q {out_prefix}.delta > {out_prefix}.delta.filtered 2>> {log}")



def parse_delta(filename):
    """Returns (alignment length, similarity errors) tuple from passed .delta.
    - filename - path to the input .delta file
    Extracts the aligned length and number of similarity errors for each
    aligned uniquely-matched region, and returns the cumulative total for
    each as a tuple.
    """
    aln_length, sim_errors = 0, 0
    for line in open(filename, 'rU'):
        fields = line.strip().split()


        if (fields[0] != 'NUCMER') and  (not fields[0].startswith('>') ) and (len(fields) == 7):
        # We only process lines with seven columns:

            aln_length += abs(int(fields[1]) - int(fields[0]))
            sim_errors += int(fields[4])
    return aln_length, sim_errors

def process_deltafiles(deltafiles, org_lengths, logger=None, **kwargs):

    Table = {'querry':[],'reference':[],'alignment_length':[],'similarity_errors':[],
            'ref_coverage':[],'querry_coverage':[],'ani':[], 'reference_length':[],
            'querry_length':[],'alignment_coverage':[]}

    # Process .delta files assuming that the filename format holds:
    # org1_vs_org2.delta
    coverage_method = kwargs.get('coverage_method')
    #logging.debug('coverage_method is {0}'.format(coverage_method))

    for deltafile in deltafiles:
        qname, sname =  os.path.splitext(os.path.split(paf_file)[-1])[0].split('-')


        tot_length, tot_sim_error = parse_delta(deltafile)
        if tot_length == 0 and logger is not None:
            logging.info("Total alignment length reported in " +
                               "%s is zero!" % deltafile)

        if not (tot_length==0):
            perc_id = 1 - float(tot_sim_error) / tot_length
            query_cover = float(tot_length) / org_lengths[qname]
            sbjct_cover = float(tot_length) / org_lengths[sname]

            Table['querry'].append(qname)
            Table['querry_length'].append(org_lengths[qname])
            Table['reference'].append(sname)
            Table['reference_length'].append(org_lengths[sname])
            Table['alignment_length'].append(tot_length)
            Table['similarity_errors'].append(tot_sim_error)
            Table['ani'].append(perc_id)
            Table['ref_coverage'].append(sbjct_cover)
            Table['querry_coverage'].append(query_cover)

            Table['alignment_coverage']= tot_length/max(org_lengths[qname],org_lengths[sname])

        # if coverage_method == 'total':
        #     Table['alignment_coverage'].append((tot_length * 2)/(org_lengths[qname]\
        #                                                      + org_lengths[sname]))
        # elif coverage_method == 'larger':
        #     Table['alignment_coverage'].append(max((tot_length/org_lengths[qname]),\
        #         (tot_length/org_lengths[sname])))

    df = pd.DataFrame(Table)
    return df






def many_mummer(genome_list,genome_folder,delta_folder,extension='.fasta'):

    os.makedirs(delta_folder,exist_ok=True)

    delta_files=[]

    genomes = open(genome_list).read().strip().split()

    for genome1 in genomes:
        for genome2 in genomes:
            deltaf= os.path.join(delta_folder,f"{genome1}-{genome2}.delta")
            ref= os.path.join(genome_folder,genome1+extension)
            query =os.path.join(genome_folder,genome2+extension)

            run_minimap(ref,query,deltaf)
            delta_files.append(deltaf)

    return delta_files



paf_files = many_minimap(snakemake.input.alignment_list,
                       snakemake.input.genome_folder,
                       snakemake.params.paf_folder,
                       snakemake.params.extension
                       )




parse_paf_files(paf_files,snakemake.input.genome_stats,snakemake.output.alignments_stats)
