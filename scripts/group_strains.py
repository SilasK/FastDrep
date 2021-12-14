import pandas as pd

import scipy.spatial as sp
import scipy.cluster.hierarchy as hc

import numpy as np
from common import genome_pdist as gd




if __name__ == "__main__":

    linkage_method = snakemake.params.linkage_method
    treshold = 0.995
    quality_score_formula = snakemake.config["quality_score"]



    M = gd.load_parquet(snakemake.input.dists)

    mag2species = pd.read_csv(snakemake.input.mag2species, index_col=0, sep="\t")
    Strains = mag2species.copy()

    M["Species"] = mag2species.loc[M.index.get_level_values(0), "Species"].values
    M["Species2"] = mag2species.loc[M.index.get_level_values(1), "Species"].values
    M = M.query("Species==Species2")



    strains_threshold= snakemake.config.get("strains_threshold",'auto')

    # simple threshold based spliting
    if strains_threshold != "auto":

        for species, data in M.groupby("Species"):

            ID = data.Identity.unstack()

            all_indexes = ID.index.union(ID.columns)
            ID = ID.reindex(index=all_indexes, columns=all_indexes).fillna(0)
            ID = ID + ID.T
            ID.values[np.eye(ID.shape[0], dtype=bool)] = 1
            Dist = 1 - ID.fillna(0.8)
            Dist.clip(0, 1, inplace=True)
            linkage = hc.linkage(sp.distance.squareform(Dist), method=linkage_method)

            labels = hc.fcluster(linkage, (1 - float(strains_threshold)) / 2, criterion="distance")

            Strains.loc[ID.index, "StrainNr"] = labels


    # Detect automatically threshold
    else:

        # detect stain gap
        fraction_below_treshold = (M.Identity < treshold).groupby(
            M.Species
        ).sum() / M.groupby("Species").size()
        have_strain_gap = fraction_below_treshold > 0.05



        for species, data in M.groupby("Species"):

            if have_strain_gap[species]:

                ID = data.Identity.unstack()

                all_indexes = ID.index.union(ID.columns)
                ID = ID.reindex(index=all_indexes, columns=all_indexes).fillna(0)
                ID = ID + ID.T
                ID.values[np.eye(ID.shape[0], dtype=bool)] = 1
                Dist = 1 - ID.fillna(0.8)
                Dist.clip(0, 1, inplace=True)
                linkage = hc.linkage(sp.distance.squareform(Dist), method=linkage_method)

                labels = hc.fcluster(linkage, (1 - treshold) / 2, criterion="distance")

                Nmax = max(labels)

                assert Nmax < 500, "Need to evaluate more than 500 tresholds"

                assert ~np.isnan(Nmax), "N range is not defined"

                Scores = gd.evaluate_clusters_range(
                    np.arange(2, Nmax + 1), Dist, linkage_method=linkage_method
                )

                N_strains = Scores.Silhouette_score.idxmax()

                if ~np.isnan(N_strains):
                    labels = hc.fcluster(linkage, N_strains, criterion="maxclust")

                Strains.loc[ID.index, "StrainNr"] = labels


    # Finalize naming

    Strains["Strain"] = (
        "Strain"
        + Strains.SpeciesNr.astype(int).astype("str")
        + "_"
        + Strains.StrainNr.dropna().astype(int).astype("str")
    )
    Strains.loc[Strains.Strain.isnull(), "Strain"] = Strains.loc[
        Strains.Strain.isnull(), "Species"
    ]

    print(f"Identified { Strains.Strain.unique().size} strains")

    Q = pd.read_csv(snakemake.input.genome_info,sep='\t',index_col=0)
    quality_score = Q.eval(quality_score_formula)

    Strains["Representative_Strain"] = gd.best_genome_from_table(
        Strains.Strain, quality_score
    )

    Strains.to_csv(snakemake.output.mag2strain, sep="\t")


    for col in ["Strain","Representative_Strain"]:
        Q[col] = Strains[col]

    Q.to_csv(snakemake.input.genome_info,sep='\t')
