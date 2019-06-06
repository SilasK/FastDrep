
genome_folder= 'all_bins'

genomes_of_species = dict(ERR675604_metabat_01=['ERR675524_metabat_01', 'ERR675522_maxbin_09', 'FR4_maxbin_13',
       'ERR675523_metabat_01', 'ERR675501_metabat_17', 'F39_maxbin_01_sub',
       'C13_maxbin_22_sub', 'ERR675514_metabat_04', 'SRR4116659_maxbin_23',
       'ERR675581_metabat_01', 'ERR675622_metabat_01', 'ERR675520_maxbin_04',
       'ERR675512_metabat_01', 'ERR675655_metabat_02', 'CR1_maxbin_07',
       'ERR675527_metabat_02', 'ERR675626_maxbin_08', 'C15_metabat_02',
       'ERR675528_maxbin_17', 'ERR675508_metabat_01', 'C4_metabat_01',
       'ERR675531_maxbin_06', 'CR8_metabat_21', 'C7_maxbin_33',
       'ERR675506_metabat_01_sub', 'C11_metabat_06', 'ERR675631_maxbin_32',
       'ERR675502_metabat_01', 'ERR675584_metabat_32', 'F40_metabat_03',
       'ERR675519_maxbin_01', 'C14_maxbin_20_sub', 'FR3_metabat_04',
       'F37_metabat_03', 'F34_metabat_02', 'ERR675510_maxbin_05',
       'ERR675634_metabat_01', 'ERR675617_metabat_05', 'ERR675515_metabat_01',
       'C16_maxbin_06', 'ERR675605_metabat_01', 'CR6_maxbin_47',
       'ERR675607_maxbin_5', 'ERR675620_maxbin_7', 'C6_metabat_05',
       'ERR675521_maxbin_13', 'ERR675629_maxbin_25', 'ERR675627_maxbin_21',
       'C1_metabat_22', 'ERR675632_metabat_10', 'ERR675633_maxbin_09',
       'F27_metabat_04', 'ERR675621_maxbin_1', 'ERR675529_maxbin_08',
       'ERR675625_metabat_1', 'C10_maxbin_37_sub', 'ERR675526_maxbin_25',
       'ERR675635_maxbin_08', 'ERR675604_metabat_01', 'C9_metabat_01',
       'ERR675525_metabat_01', 'ERR675516_metabat_01', 'ERR675624_metabat_1',
       'C2_maxbin_13', 'ERR675500_maxbin_06', 'ERR675530_maxbin_11',
       'ERR675505_metabat_25', 'F33_metabat_10', 'F36_maxbin_20_sub',
       'ERR675517_metabat_08', 'F35_metabat_04', 'C3_maxbin_17',
       'ERR675507_metabat_01', 'C8_metabat_01', 'ERR675582_metabat_02',
       'ERR675583_metabat_01', 'ERR675511_metabat_01', 'ERR675504_metabat_01',
       'CR7_maxbin_51', 'C5_metabat_24', 'F32_maxbin_14_sub',
       'ERR675628_maxbin_24', 'ERR675518_maxbin_24', 'CR2_maxbin_21',
       'ERR675509_metabat_01', 'ERR675665_maxbin_02', 'ERR675636_maxbin_18',
       'ERR675503_maxbin_21', 'ERR675606_maxbin_2', 'F30_metabat_03',
       'FR6_metabat_27', 'ERR675630_maxbin_16', 'FR2_metabat_004',
       'FR1_maxbin_03', 'ERR675498_maxbin_05', 'ERR675619_maxbin_02',
       'C12_maxbin_33_sub', 'ERR675618_metabat_01', 'F38_metabat_08'])

rule all:
    input:
        "minimap/{species}/alignments_stats.tsv".format(species='ERR675604_metabat_01')

rule minimap:
    input:
        ref= f"{genome_folder}/{{genome2}}.fasta",
        querry=f"{genome_folder}/{{genome1}}.fasta",
    output:
        "minimap/paf/{genome1}-{genome2}.paf"
    threads:
        3
    params:
        preset= "asm10" #asm5/asm10/asm20: asm-to-ref mapping, for ~0.1/1/5% sequence divergence
    shell:
        "minimap2 -x {params.preset} -t {threads} {input.querry} {input.ref}   > {output}"

import pandas as pd
def load_paf(paf_file):

    M= pd.read_csv(paf_file,header=None,sep='\t',usecols=range(12))
    M.columns=['Contig2','Length2','Start2','End2',
               'Strand',
               'Contig1','Length1','Start1','End1',
               'Nmatches','Allength','Quality']
    M['Identity']= M.Nmatches/M.Allength

    return M

from itertools import combinations
def parse_paf_input(wildcards):

    genomes = genomes_of_species[wildcards.species]
    return [f"minimap/paf/{g1}-{g2}.paf" for g1,g2 in combinations(genomes,2)]


localrules: parse_paf
rule parse_paf:
    input:
        parse_paf_input
    output:
        "minimap/{species}/alignments_stats.tsv"
    run:
        Out=[]
        for paf_file in input:

            M= load_paf(paf_file)

            parsed= {}
            M.sort_values('Identity',inplace=True,ascending=False)

            parsed['average_id'] = (M.Allength*M.Identity).sum() / M.Allength.sum()
            parsed['Id90'] = M.loc[0.9*M.Allength.sum() <= M.Allength.cumsum()].iloc[0].Identity
            parsed['allength99'] = M.query('Identity>=0.99').Allength.sum()
            parsed['allength95'] = M.query('Identity>=0.95').Allength.sum()

            parsed= pd.Series(parsed,name=(wildcards.genome1,wildcards.genome2))

            Out.append(parsed)


        pd.concat(Out).to_csv(output[0],sep='\t')
