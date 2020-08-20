import sys, os
import networkx as nx


sys.stdout = open(snakemake.log[0], "w")
sys.stderr = open(snakemake.log[0], "a")

from common import genome_pdist as gd
import pandas as pd

Q= gd.load_quality(snakemake.input.quality)
Q['quality_score']= Q.eval(snakemake.config['quality_score'])

Msh= gd.load_mash(snakemake.input.dists)
Msh_high= Msh.query(f"Identity>={snakemake.params.treshold}")
G= gd.to_graph(Msh_high)

#map genomes to clusters by single linkage
Clustering= {}
for i,cc in enumerate(gd.nx.connected_components(G)):
    for g in cc:
          Clustering[g]=i
Clustering =gd.best_genome_from_table(pd.Series(Clustering),Q.quality_score)

# don't foget  Genomes qith no connection
lonly_genomes= Q.index.difference(Clustering.index)
Clustering= Clustering.append(pd.Series(lonly_genomes,lonly_genomes))


#save clustering
Clustering.name='Precluster_representative'
Clustering.index.name='Genome'
Clustering.to_csv(snakemake.output.preclustering,sep='\t')

#Subset_pairwise distances
Cluster_representatives= Clustering.unique()

Msh_sub= Msh.loc[ (
Msh.index.levels[0].intersection(Cluster_representatives),
Msh.index.levels[1].intersection(Cluster_representatives)
 ),:].query(f"Identity>={snakemake.params.min_identity}")


# save list of all comparison to perform alignment on them
G= gd.to_graph(Msh_sub)
G.remove_edges_from(nx.selfloop_edges(G))

N_comparisons= G.number_of_edges()

assert N_comparisons>=1,"No connection correspond the criteria for preclustering"

nx.write_edgelist(G,snakemake.output.edgelist,delimiter='\t',data=False,comments=None)

with open(snakemake.log.stats,'w') as logfile:
    logfile.write(f"From {Clustering.shape[0]} genomes {Cluster_representatives.shape[0]} "
      f"({Cluster_representatives.shape[0]/Clustering.shape[0]*100:.2f}%) are representatives.\n"
      f"This decreases the number of interaction to {N_comparisons:d} ({N_comparisons/Msh.shape[0]:.2g})\n"
     )
