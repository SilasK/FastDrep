import pandas as pd
def load_paf(paf_file):

    M= pd.read_csv(paf_file,header=None,sep='\t',usecols=range(12))
    M.columns=['Contig2','Length2','Start2','End2',
               'Strand',
               'Contig1','Length1','Start1','End1',
               'Nmatches','Allength','Quality']
    M['Identity']= M.Nmatches/M.Allength

    return M
