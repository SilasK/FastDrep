import pandas as pd


import networkx as nx
import warnings
import os
import scipy.spatial as sp
import scipy.cluster.hierarchy as hc
from sklearn.metrics import silhouette_score
import numpy as np


def simplify_index(index):
    "assumes single index are path of files, removes extesnion and dirname"

    path = index[0]
    filename, extension = os.path.splitext(path)

    if extension == ".gz":
        extension = os.path.splitext(filename)[-1] + extension

    dir_name = os.path.dirname(path) + "/"

    return pd.Index(index).str.replace(extension, "").str.replace(dir_name, "")


def load_ani_table_(dist_file, header=None, simplify_names=False):

    F = pd.read_csv(dist_file, sep="\t", header=None, index_col=[0, 1])

    if header is not None:
        F.columns = header
    F.index.names = ["Genome1", "Genome2"]

    if simplify_names:

        F.index = pd.MultiIndex(
            levels=[
                simplify_index(F.index.levels[0]),
                simplify_index(F.index.levels[1]),
            ],
            codes=F.index.codes,
        )

    return F


def load_bbsketch(dist_file, format=3, simplify_names=True):
    """ reads output of sendsketch.sh
        format=3 [query,ref,ANI..]
        format=2 Table for one query
            parses parameters in first line returns df,params
    """

    if format == 3:
        bbs = pd.read_csv(dist_file, index_col=[0, 1], sep="\t")
        bbs.index.names = ["Genome1", "Genome2"]
        if (bbs.QTaxID == -1).all():
            bbs.drop(["QTaxID", "RTaxID"], axis=1, inplace=True)

        bbs["Identity"] = bbs.iloc[:, 0] / 100.0

        if "SSU" in bbs:
            bbs["SSU"] = bbs.SSU.replace(".", np.nan)

        if simplify_names:

            bbs.index = pd.MultiIndex(
                levels=[
                    simplify_index(bbs.index.levels[0]),
                    simplify_index(bbs.index.levels[1]),
                ],
                codes=bbs.index.codes,
            )

        return bbs
    elif format == 2:

        f = open(send_sketch_file)
        f.readline()  # trash empty line
        comment_line = f.readline().strip()
        params = dict(key_value.split(":") for key_value in comment_line.split("\t"))

        df = pd.read_csv(f, sep="\t")

        convert_percentages(df)

        return df, params
    else:
        raise NotImplementedError(
            "I don't know how to parse other formats than 2,3 of bbsketch"
        )


def load_fastani(dist_file, simplify_names=True):
    """Loads fastANI output calculates overlap.
    Outputs a table with ['Genome1','Genome2','Identity','Nmapped','Ntotal','Overlap' ] in header"""

    F = load_ani_table_(
        dist_file, ["Identity", "Nmapped", "Ntotal"], simplify_names=simplify_names
    )
    F.loc[:, "Overlap"] = F.Nmapped.values / F.Ntotal.values
    F["Identity"] /= 100.0

    return F


def load_mash(dist_file, simplify_names=True):
    """Loads fastANI output calculates overlap.
    Outputs a table with ['Genome1','Genome2','Distance','Pvalue','Fraction','Identity'] in header"""
    F = load_ani_table_(
        dist_file, ["Distance", "Pvalue", "Fraction"], simplify_names=simplify_names
    )

    F["Identity"] = 1 - F.Distance
    F["Nmapped"] = F.Fraction.map(lambda s: int(s.split("/")[0]))
    F["Ntotal"] = F.Fraction.map(lambda s: int(s.split("/")[1]))
    F["Fraction"] = F.Nmapped / F.Ntotal

    return F


def load_parquet(parquet_file):

    M= pd.read_parquet(parquet_file,columns=["Distance"])
    M['Identity']= 1-M.Distance
    return M

def load_bindash(dist_file, simplify_names=True):
    """Loads bindash output.
    Outputs a table with
    ['Genome1','Genome2','Distance','Pvalue','Fraction','Nmapped','Ntotal','Identity']
    in header.

    Bindash tables are not necessarily simetrical.
    """
    F = load_ani_table_(
        dist_file, ["Distance", "Pvalue", "Fraction"], simplify_names=simplify_names
    )

    F["Nmapped"] = F.Fraction.map(lambda s: int(s.split("/")[0])).astype(int)
    F["Ntotal"] = F.Fraction.map(lambda s: int(s.split("/")[1])).astype(int)
    F["Fraction"] = F.Nmapped / F.Ntotal
    F["Identity"] = 1 - F.Distance

    return F


def load_mummer(dist_file):

    M = pd.read_csv(dist_file, sep="\t", index_col=[0, 1])
    M["Identity"] = M.ANI
    return M


def load_minimap(dist_file):

    M = pd.read_csv(dist_file, sep="\t", index_col=[0, 1])
    assert "Identity" in M.columns
    return M


def to_graph(F, attributes=None, **kws):

    df = F.copy()

    df["Genome1"] = df.index.get_level_values(0)
    df["Genome2"] = df.index.get_level_values(1)

    G = nx.from_pandas_edgelist(df, "Genome1", "Genome2", attributes, **kws)

    return G


def evaluate_clusters(labels, Dist):

    try:
        Silhouette_score = silhouette_score(Dist, metric="precomputed", labels=labels)
        N_clusters = np.unique(labels).shape[0]
    except ValueError:
        Silhouette_score, N_clusters = np.nan, np.nan

    return Silhouette_score, N_clusters


def evaluate_clusters_range(
    N_range, Dist, linkage_method="average", criterion="maxclust"
):

    linkage = hc.linkage(sp.distance.squareform(Dist), method=linkage_method)

    def get_clusters(N):
        labels = hc.fcluster(linkage, N, criterion=criterion)
        return evaluate_clusters(labels, Dist)

    Scores = pd.DataFrame(
        [get_clusters(t) for t in N_range],
        index=N_range,
        columns=["Silhouette_score", "N_clusters"],
    )

    return Scores


def evaluate_clusters_thresholds(
    thresholds, Dist, linkage_method="average", criterion="distance"
):

    linkage = hc.linkage(sp.distance.squareform(Dist), method=linkage_method)

    def get_clusters(t):
        labels = hc.fcluster(linkage, 1 - t, criterion=criterion)
        return evaluate_clusters(labels, Dist)

    Scores = pd.DataFrame(
        [get_clusters(t) for t in thresholds],
        index=thresholds,
        columns=["Silhouette_score", "N_clusters"],
    )

    return Scores


def plot_scores(Scores, xlabel="Treshold"):

    import matplotlib.pylab as plt

    f, axe = plt.subplots(2, 1, sharex=True, figsize=(6, 5))
    Scores.Silhouette_score.plot(marker=".", ax=axe[0])

    axe[0].set_ylabel("Silhouette score")
    Scores.N_clusters.plot(marker=".", ax=axe[1])
    axe[1].set_ylabel("N clusters")

    axe[1].set_xlabel(xlabel)

    return f, axe


def group_species_linkage(
    M, threshold=0.95, fillna=0.8, linkage_method="average", square=False
):

    assert threshold > 0.3, "threshold is an identity value"

    cutoff = 1 - threshold

    ID = M.Identity.unstack()
    all_index = ID.index.union(ID.columns)
    ID = ID.reindex(index=all_index, columns=all_index)

    # take smaler of both comparisons (fastANI)
    # ID= ID+(ID.T-ID).applymap(lambda s: min(s,0))
    # ID.values[np.eye(ID.shape[0],dtype=bool)]=1

    Dist = 1 - ID.fillna(fillna)
    if square:
        cutoff = cutoff ** 2
        Dist = Dist.pow(2)

    linkage = hc.linkage(sp.distance.squareform(Dist), method=linkage_method)
    labels = pd.Series(
        hc.fcluster(linkage, cutoff, criterion="distance"), index=Dist.index
    )

    return labels


def load_quality(checkm_file):
    Q = pd.read_csv(checkm_file, index_col=0, sep="\t")
    Q = Q.rename(
        columns={
            "Strain heterogeneity": "strain_heterogeneity",
            "strain heterogeneity": "strain_heterogeneity",
            "Contamination": "contamination",
            "Completeness": "completeness",
        }
    )
    Q.index = Q.index.str.replace(".fasta", "")

    return Q


def best_genome_from_table(Grouping, quality_score):

    Mapping = pd.Series(index=Grouping.index)

    for group in Grouping.unique():
        genomes = Grouping.index[Grouping == group]
        representative = quality_score.loc[genomes].idxmax()
        Mapping.loc[genomes] = representative

    return Mapping


def clustermap(DistanceMatrix, linkage_method="average", **kws):
    import seaborn as sns
    import scipy.spatial as sp, scipy.cluster.hierarchy as hc

    linkage = hc.linkage(sp.distance.squareform(DistanceMatrix), method=linkage_method)

    cg = sns.clustermap(
        1 - DistanceMatrix, row_linkage=linkage, col_linkage=linkage, **kws
    )

    return cg


def pairewise2matrix(M, column="Identity", fillna=np.nan):
    """
        This functions turns a pairewise genome distance table [genome1, genome2, column...]
        In to a matrix [genome1 genome2] of the values of column.
        When ANI values are symetrical (with minimal error),
        usually only one halve of NxN possibilities values are calculated.


        Diagonal values are set to 1

    """

    ID = M[column].unstack()

    all_indexes = ID.index.union(ID.columns)

    ID = ID.reindex(index=all_indexes, columns=all_indexes)
    ID = ID.fillna(0)
    ID = ID + ID.T
    ID.values[np.eye(ID.shape[0], dtype=bool)] = 1
    return ID.replace(0, fillna)
