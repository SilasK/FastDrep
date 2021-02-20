import sys, os
from snakemake import execute_snakemake


sys.stdout = open(snakemake.log[0], "w")
sys.stderr = open(snakemake.log[0], "a")

execute_snakemake( snakemake.params.snakefile,
          config={comparison_list=snakemake.input.comparison_list,
              genome_folder=snakemake.input.genome_folder,
              subset=snakemake.wildcards.subset,
              genome_stats=snakemake.input.genome_stats,
              tmpfolder=snakemake.config['tmpfolder']
              },
          cores=snakemake.threads,
          lock=False,
          force_incomplete=True
          )
