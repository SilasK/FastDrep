import pandas as pd

import scipy.spatial as sp
import scipy.cluster.hierarchy as hc
from sklearn.metrics  import silhouette_score
import numpy as np
from common import genome_pdist as gd




def automatic_cluster_species(Dist,seed_tresholds= [0.9,0.95],linkage_method='average'):

    linkage = hc.linkage(sp.distance.squareform(Dist), method=linkage_method)

    def get_Nclusters(treshold):
        labels= hc.fcluster(linkage,(1-treshold),criterion='distance')
        return max(labels)



    N_range= [get_Nclusters(t) for t in seed_tresholds]



    assert (N_range[1]-N_range[0])< 60, "Need to evaluate more than 60 tresholds"

    assert ~np.isnan(N_range).any(), "N range is not defined"

    Scores= gd.evaluate_clusters_range(np.arange(min(N_range),max(N_range)+1),Dist)


    if N_range[0]==N_range[1]:
        labels= hc.fcluster(linkage,(1-seed_tresholds[0]),criterion='distance')
    else:

        N_species= Scores.Silhouette_score.idxmax()



        labels= hc.fcluster(linkage,N_species,criterion='maxclust')

    print(f"Identified { max(labels)} species")



    return Scores,labels


if __name__=='__main__':

    M= gd.load_mash(snakemake.input.dists)

    ID= M.Identity.unstack()
    ID= ID.loc[ID.index,ID.index]
    Dist= 1-ID.fillna(0.8)


    Scores,labels= automatic_cluster_species(Dist)

    labels= pd.Series(labels,name='Species')
    labels.index.name='genome'
    n_leading_zeros= len(str(max(labels)))
    format_int='Species{:0'+str(n_leading_zeros)+'d}'
    labels=labels.apply(format_int.format)


    labels.to_csv(snakemake.output.cluster_file,sep='\t',header=False)
    Scores.to_csv(snakemake.output.scores,sep='\t')
