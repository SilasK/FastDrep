import pandas as pd


import networkx as nx
import warnings
import os
import scipy.spatial as sp
import scipy.cluster.hierarchy as hc
from sklearn.metrics  import silhouette_score
import numpy as np

def simplify_index(index):
    "assumes single index are path of files, removes extesnion and dirname"


    path= index[0]
    extension= os.path.splitext(path)[-1]
    dir_name= os.path.dirname(path)+'/'

    return index.str.replace(extension,'').str.replace(dir_name,'')


def load_ani_table_(dist_file,header=None,simplify_names=False):

    F = pd.read_csv(dist_file,sep='\t',header=None,index_col=[0,1])

    F.columns= header
    F.index.names=['Genome1','Genome2']

    if simplify_names:

        F.index =pd.MultiIndex(levels= [simplify_index(F.index.levels[0]),
                              simplify_index(F.index.levels[1])],codes= F.index.codes  )




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


def evaluate_clusters(labels,Dist):

        try:
            Silhouette_score= silhouette_score(Dist, metric='precomputed', labels=labels)
            N_clusters = np.unique(labels).shape[0]
        except ValueError:
            Silhouette_score, N_clusters = np.nan, np.nan

        return Silhouette_score, N_clusters

def evaluate_clusters_range(N_range,Dist,linkage_method='average',criterion='maxclust'):

    linkage = hc.linkage(sp.distance.squareform(Dist), method=linkage_method)

    def get_clusters(N):
        labels= hc.fcluster(linkage,N,criterion=criterion)
        return evaluate_clusters(labels,Dist)

    Scores= pd.DataFrame([get_clusters(t) for t in N_range],index= N_range, columns= ['Silhouette_score','N_clusters'])

    return Scores


def evaluate_clusters_tresholds(tresholds,Dist,linkage_method='average',criterion='distance'):

    linkage = hc.linkage(sp.distance.squareform(Dist), method=linkage_method)

    def get_clusters(t):
        labels= hc.fcluster(linkage,1-t,criterion=criterion)
        return evaluate_clusters(labels,Dist)

    Scores= pd.DataFrame([get_clusters(t) for t in tresholds],index= tresholds, columns= ['Silhouette_score','N_clusters'])

    return Scores


def plot_scores(Scores,xlabel='Treshold'):

    import matplotlib.pylab as plt

    f,axe= plt.subplots(2,1,sharex=True,figsize=(6,5))
    Scores.Silhouette_score.plot(marker='.',ax=axe[0])

    axe[0].set_ylabel('Silhouette score')
    Scores.N_clusters.plot(marker='.',ax=axe[1])
    axe[1].set_ylabel('N clusters')

    axe[1].set_xlabel(xlabel)


    return f,axe




def group_species_linkage(M,threshold = 0.95,fillna=0.8,linkage_method='average',square=False):

    assert threshold>0.3, "threshold is an identity value"

    cutoff = (1-threshold)


    ID= M.Identity.unstack()
    all_index=ID.index.union(ID.columns)
    ID= ID.reindex(index=all_index,columns=all_index)

    #take smaler of both comparisons (fastANI)
    # ID= ID+(ID.T-ID).applymap(lambda s: min(s,0))
    # ID.values[np.eye(ID.shape[0],dtype=bool)]=1


    Dist= (1-ID.fillna(fillna))
    if square:
        cutoff= cutoff**2
        Dist=Dist.pow(2)

    linkage = hc.linkage(sp.distance.squareform(Dist), method=linkage_method)
    labels= pd.Series(hc.fcluster(linkage,cutoff,criterion='distance'), index= Dist.index)



    return labels

def load_quality(checkm_file):
    Q= pd.read_csv(checkm_file, index_col=0,sep='\t')
    Q= Q.rename(columns={'strain heterogeneity':'strain_heterogeneity'})
    Q.index= Q.index.str.replace('.fasta','')

    return Q

def best_genome_from_table(Grouping,quality_score):

    Mapping = pd.Series(index=Grouping.index)

    for group in Grouping.unique():
        genomes= Grouping.index[Grouping==group]
        representative= quality_score.loc[genomes].idxmax()
        Mapping.loc[genomes]=representative

    return Mapping




def clustermap(DistanceMatrix,linkage_method='average',**kws):
    import seaborn as sns
    import scipy.spatial as sp, scipy.cluster.hierarchy as hc

    linkage = hc.linkage(sp.distance.squareform(DistanceMatrix), method=linkage_method)

    cg=sns.clustermap(1-DistanceMatrix,
               row_linkage= linkage, col_linkage=linkage,
               **kws
                )

    return cg
