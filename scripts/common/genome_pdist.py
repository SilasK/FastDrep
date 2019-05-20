import pandas as pd

import seaborn as sns
import networkx as nx
import warnings
import os




def load_ani_table_(dist_file,header=None,simplify_indexes=False):

    F = pd.read_csv(dist_file,sep='\t',header=None,index_col=[0,1])

    F.columns= header
    F.index.names=['Genome1','Genome2']

    if simplify_indexes:

        path= F.index[0][0]
        extension= os.path.splitext(path)[-1]
        dir_name= os.path.dirname(path)+'/'

        def simplify(index):
            return index.replace(extension,'').replace(dir_name,'')

        F= F.rename(index=simplify)




    return F








def load_fastani(dist_file,simplify_indexes=True):
    """Loads fastANI output calculates overlap.
    Outputs a table with ['Genome1','Genome2','Identity','Nmapped','Ntotal','Overlap' ] in header"""

    F = load_ani_table_(dist_file,['Identity','Nmapped','Ntotal'],simplify_indexes=simplify_indexes)
    F.loc[:,'Overlap']= F.Nmapped.values/F.Ntotal.values
    F['Identity']/=100.

    return F

def load_mash(dist_file,simplify_indexes=True):
    """Loads fastANI output calculates overlap.
    Outputs a table with ['Genome1','Genome2','Distance','Pvalue','Fraction','Identity'] in header"""
    F = load_ani_table_(dist_file,['Distance','Pvalue','Fraction'],simplify_indexes=simplify_indexes)

    F['Identity']= 1- F.Distance


    return F



def to_graph(F,attributes=None,**kws):

    df= F.copy()

    df['Genome1']= df.index.get_level_values(0)
    df['Genome2']= df.index.get_level_values(1)

    G= nx.from_pandas_edgelist(df,'Genome1','Genome2',attributes,**kws)

    return G




def get_connected_components(Graph):

    try:
        next(Graph.selfloop_edges())
    except StopIteration:
        warings.warn("The graph does't contain self loops, cluster of size 1 will not be retreaved")


    return list(nx.connected_components(Graph))


def map_to_best_genome(G,quality_score):


    CC= list(nx.connected_components(G))

    Mapping=pd.Series(index=quality_score.index)

    for c in CC:
        representative= quality_score.loc[c].idxmax()
        Mapping.loc[c]=representative

    single_genomes=[n for n,d in G.degree() if d==0 ]
    if len(single_genomes)>0:
        Mapping.loc[single_genomes]=single_genomes


    return Mapping

def internalize_index(index, genomefolder,extension='.fasta'):
    """adapts indexes of genomes to represent the full path to genome, this is for internal use
    """


    if extension in index:
        raise Exception(f"Index has already an extension prpbbly, you don't want to internalize_indexes. {i}")

    return os.path.join(genomefolder,index+extension)

def extarnilize_index(indexes,genomefolder,extension='.fasta'):

    return index.replace(extension,'').replace(genomefolder,'')





def clustermap(Identity):
    import seaborn as sns
    import scipy.spatial as sp, scipy.cluster.hierarchy as hc

    D= Identity.unstack()
    linkage = hc.linkage(1-D.fillna(1), method='single')

    cg= sns.clustermap(D.fillna(0), row_linkage=linkage, col_linkage=linkage)
    return cg
