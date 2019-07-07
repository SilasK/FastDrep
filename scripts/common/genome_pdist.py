import pandas as pd


import networkx as nx
import warnings
import os


def simplify_indexes(df):
    "assumes indexes are path of files, removes extesnion and dirname"

    path= df.index[0][0]
    extension= os.path.splitext(path)[-1]
    dir_name= os.path.dirname(path)+'/'

    def simplify(index):
        return index.replace(extension,'').replace(dir_name,'')

    return df.rename(index=simplify)

def load_ani_table_(dist_file,header=None,simplify_names=False):

    F = pd.read_csv(dist_file,sep='\t',header=None,index_col=[0,1])

    F.columns= header
    F.index.names=['Genome1','Genome2']

    if simplify_names:

        F= simplify_indexes(F)




    return F








def load_fastani(dist_file,simplify_names=True):
    """Loads fastANI output calculates overlap.
    Outputs a table with ['Genome1','Genome2','Identity','Nmapped','Ntotal','Overlap' ] in header"""

    F = load_ani_table_(dist_file,['Identity','Nmapped','Ntotal'],simplify_names=simplify_names)
    F.loc[:,'Overlap']= F.Nmapped.values/F.Ntotal.values
    F['Identity']/=100.

    return F

def load_mash(dist_file,simplify_names=True):
    """Loads fastANI output calculates overlap.
    Outputs a table with ['Genome1','Genome2','Distance','Pvalue','Fraction','Identity'] in header"""
    F = load_ani_table_(dist_file,['Distance','Pvalue','Fraction'],simplify_names=simplify_names)

    F['Identity']= 1- F.Distance
    F['Nmapped']=F.Fraction.map(lambda s: int(s.split('/')[0]))
    F['Ntotal']=F.Fraction.map(lambda s: int(s.split('/')[1]))
    F['Fraction']=F.Nmapped/F.Ntotal


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



def clustermap(DistanceMatrix,linkage_method='average',**kws):
    import seaborn as sns
    import scipy.spatial as sp, scipy.cluster.hierarchy as hc

    linkage = hc.linkage(sp.distance.squareform(DistanceMatrix), method=linkage_method)

    cg=sns.clustermap(1-DistanceMatrix,
               row_linkage= linkage, col_linkage=linkage,
               **kws
                )

    return cg
