import pandas as pd
import os



MINIMAP_HEADERS=['Contig2','Length2','Start2','End2',
               'Strand',
               'Contig1','Length1','Start1','End1',
               'Nmatches','Allength','Quality']
MINIMAP_DATATYPES= [str,int,int,int,str,str,int,int,int,int,int,int]

assert len(MINIMAP_HEADERS)==len(MINIMAP_DATATYPES)



minimap_dtypes_map= {'i':int,'f':float}


def parse_minimap_tag(tag,out={}):

    name,dtype,value = tag.split(':')

    dtype=minimap_dtypes_map.get(dtype,str)

    out[name]= dtype(value)

def parse_minimap_line(line):
    """parses a minmap paf line, return a dict.
    reads tags and converts datatyes"""
    elements= line.strip().split()
    print([e[:3] for e in elements])
    out={}

    if not len(elements)==0:

        try:




            for i,h in enumerate(MINIMAP_HEADERS):


                dtype = MINIMAP_DATATYPES[i]
                out[h]= dtype(elements[i])



            for i in range(len(MINIMAP_HEADERS),len(elements)):

                parse_minimap_tag(elements[i],out)


        except Exception as e:
            raise IOError(f'Error during parsing paf line : {elements}') from e

    return out


def load_paf(paf_file):

    try:

        parsed=[]
        with open(paf_file) as f:
            for line in f:
                line= f.readline()
                parsed.append(parse_minimap_line(line))

        M=pd.DataFrame(parsed).dropna(how='all')

        # some values are 0.0000, some are negative
        M['Identity']= 1- M.de #.replace('0.0000',0.00005).astype(float).apply(lambda d: max(d,0))

        headers=MINIMAP_HEADERS+['Identity']
        #rearange headers
        M=M.loc[:, headers+list(M.columns.drop(headers))]


    except Exception as e:
        raise IOError(f'Error during parsing paf file: {paf_file}') from e

    return M

def parse_paf_files(paf_files,genome_stats_file,output_file):

    stats= pd.read_csv(genome_stats_file,index_col=0,sep='\t')

    with open(output_file,'w') as out:

        out.write("\t".join(['genome1','genome2','Identity','Length',
                             'Length_at99id','Length_at95id','Id_at90length','Id_at50length',
                             'AlignedFraction','AlignedFraction95','AlignedFraction99'])+'\n')

        for paf_file in paf_files:

            M= load_paf(paf_file)


            M.sort_values('Identity',inplace=True,ascending=False)

            genome1,genome2 = os.path.splitext(os.path.split(paf_file)[-1])[0].split('-')
            Identity = (M.Identity*M.Allength).sum() / M.Allength.sum()
            Length= M.Allength.sum()

            min_genome_size= stats.loc[[genome1,genome2],'Total_length'].min()
            AlignedFraction =Length / min_genome_size



            Length_at99id = M.query('Identity>=0.99').Allength.sum()
            Length_at95id = M.query('Identity>=0.95').Allength.sum()


            AlignedFraction95= Length_at95id / min_genome_size
            AlignedFraction99= Length_at99id / min_genome_size


            Id_at90length = M.loc[0.9*M.Allength.sum() <= M.Allength.cumsum()].iloc[0].Identity
            Id_at50length = M.loc[0.5*M.Allength.sum() <= M.Allength.cumsum()].iloc[0].Identity

            out.write(f"{genome1}\t{genome2}\t{Identity}\t{Length}\t"
                      f"{Length_at99id}\t{Length_at95id}\t{Id_at90length}\t{Id_at50length}\t"
                      f"{AlignedFraction}\t{AlignedFraction95}\t{AlignedFraction99}\n")
