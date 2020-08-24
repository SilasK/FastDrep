import sys, os
import networkx as nx


sys.stdout = open(snakemake.log[0], "w")
sys.stderr = open(snakemake.log[0], "a")

from common import genome_pdist as gd
import pandas as pd

Q= gd.load_quality(snakemake.input.quality)
Q['quality_score']= Q.eval(snakemake.config['quality_score'])

Msh= gd.load_mash(snakemake.input.dists)
Msh= Msh.query(f"Identity>={snakemake.params.treshold}")
G= gd.to_graph(Msh)

#map genomes to clusters by single linkage
Clustering= {}
for i,cc in enumerate(gd.nx.connected_components(G)):
    for g in cc:
          Clustering[g]=i
Clustering =gd.best_genome_from_table(pd.Series(Clustering),Q.quality_score)

# don't foget  Genomes with no connection
lonly_genomes= Q.index.difference(Clustering.index)
Clustering= Clustering.append(pd.Series(lonly_genomes,lonly_genomes))


#save clustering
Clustering.name='Precluster_representative'
Clustering.index.name='Genome'
Clustering.to_csv(snakemake.output.preclustering,sep='\t')

#Subset_pairwise distances
Cluster_representatives= Clustering.unique()

with open(snakemake.log.stats,'w') as logfile:
    logfile.write(f"From {Clustering.shape[0]} genomes {Cluster_representatives.shape[0]} "
      f"({Cluster_representatives.shape[0]/Clustering.shape[0]*100:.2f}%) are representatives.\n"
     )

assert len(Cluster_representatives)>1, "less than two precluster representatives. All your genomes might belong to the same strain/subspecies."

# write list of genomes of Cluster representatives

fasta_files=[os.path.join(snakemake.params.genome_folder,f"{genome}{snakemake.config['fasta_extension']}") for genome in Cluster_representatives]

with open(snakemake.output.cluster_list,'w') as f:
    f.write('\n'.join(fasta_files)+'\n' )
