import pandas as pd

import scipy.spatial as sp
import scipy.cluster.hierarchy as hc

import numpy as np
from common import genome_pdist as gd




if __name__=='__main__':

    linkage_method= snakemake.params.linkage_method
    treshold = 0.995
    quality_score_formula = snakemake.config['quality_score']
    M= pd.read_csv(snakemake.input.dists,index_col=[0,1],sep='\t')
    mag2species= pd.read_csv(snakemake.input.mag2species,index_col=0,sep='\t')


    fraction_below_treshold= (M.ANI<treshold).groupby(M.Species).sum() / M.groupby('Species').size()
    have_strain_gap= (fraction_below_treshold >0.05)

    Strains= mag2species.copy()

    for species,data in M.groupby('Species'):

        if have_strain_gap[species]:

            ID= data.ANI.unstack()

            all_indexes= ID.index.union(ID.columns)
            ID= ID.reindex(index=all_indexes, columns=all_indexes).fillna(0)
            ID= ID+ID.T
            ID.values[np.eye(ID.shape[0],dtype=bool)]=1
            Dist= 1-ID.fillna(0.8)
            linkage = hc.linkage(sp.distance.squareform(Dist), method=linkage_method)



            labels= hc.fcluster(linkage,(1-treshold)/2,criterion='distance')

            Nmax= max(labels)

            assert Nmax< 60, "Need to evaluate more than 60 tresholds"

            assert ~np.isnan(Nmax), "N range is not defined"


            Scores= gd.evaluate_clusters_range(np.arange(2,Nmax+1),Dist,linkage_method=linkage_method)

            N_strains= Scores.Silhouette_score.idxmax()

            if ~np.isnan(N_strains):
                labels= hc.fcluster(linkage,N_strains,criterion='maxclust')

            Strains.loc[ID.index,'StrainNr']= labels

    Strains['Strain']='Strain'+Strains.SpeciesNr.astype(int).astype('str')+'_'+Strains.StrainNr.dropna().astype(int).astype('str')
    Strains.loc[Strains.Strain.isnull(),'Strain']= Strains.loc[Strains.Strain.isnull(),'Species']


    print(f"Identified { Strains.Strain.unique().size} strains")



    Q= gd.load_quality(snakemake.input.quality)
    quality_score= Q.eval(quality_score_formula)

    Strains['Representative_Strain']=gd.best_genome_from_table(Strains.Strain,quality_score)

    Strains.to_csv(snakemake.output.mag2strain,sep='\t')
