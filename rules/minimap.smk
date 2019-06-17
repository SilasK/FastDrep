


Mapping= pd.Series()


rule minimap:
    input:
        ref= f"{genome_folder}/{{genome2}}.fasta",
        querry=f"{genome_folder}/{{genome1}}.fasta",
    output:
        "minimap/paf/{genome1}-{genome2}.paf"
    conda:
        "../envs/minimap2.yaml"
    group:
        "minimap2"
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


def parse_paf_input(wildcards):

    alignment_list = checkpoints.filter_mash.get().output[0]
    paf_files=[]

    with open(alignment_list) as f:
        for line in f:
            g1,g2= f.strip().split()
            paf_files.append(f"minimap/paf/{g1}-{g2}.paf")
    return paf_files




localrules: parse_paf
rule parse_paf:
    input:
        paf= parse_paf_input,
        genome_stats= "genome_stats.tsv"
    output:
        "minimap/{species}_alignments.tsv"
    group:
        "minimap2"
    run:

        stats= pd.read_table(input.genome_stats,index_col=0)

        with open(output[0],'w') as out:

            out.write("\t".join(['genome1','genome2','Identity','Length',
                                 'Length_at99id','Length_at95id','Id_at90length','Id_at50length',
                                 'AlignedFraction','AlignedFraction95','AlignedFraction99'])+'\n')

            for paf_file in input.paf:

                M= load_paf(paf_file)

                M.sort_values('Identity',inplace=True,ascending=False)

                genome1,genome2 = os.path.splitext(os.path.split(paf_file)[-1])[0].split('-')
                Identity = M.Nmatches.sum() / M.Allength.sum()
                Length= M.Allength.sum()

                min_genome_size= stats.loc[[genome1,genome2],'Total_length'].min()


                Length_at99id = M.query('Identity>=0.99').Allength.sum()
                Length_at95id = M.query('Identity>=0.95').Allength.sum()

                AlignedFraction =Length / min_genome_size
                AlignedFraction95= Length_at95id / min_genome_size
                AlignedFraction99= Length_at99id / min_genome_size


                Id_at90length = M.loc[0.9*M.Allength.sum() <= M.Allength.cumsum()].iloc[0].Identity
                Id_at50length = M.loc[0.5*M.Allength.sum() <= M.Allength.cumsum()].iloc[0].Identity

                out.write(f"{genome1}\t{genome2}\t{Identity}\t{Length}\t"
                          f"{Length_at99id}\t{Length_at95id}\t{Id_at90length}\t{Id_at50length}\t"
                          f"{AlignedFraction}\t{AlignedFraction95}\t{AlignedFraction99}\n")
