#!/usr/bin/env python
import os, shutil
import pandas as pd
import numpy as np
from glob import glob
import yaml


try:
    from tqdm import tqdm
except ImportError:
    tqdm=list


def list_genomes(folder,fasta_extension):
    "Return list of genome names basenames of files in folder and gven extension"

    Nchar_extension= len(fasta_extension)
    return [f[:-Nchar_extension] for f in os.listdir(folder) if f.endswith(fasta_extension)]


def validate_db(folder,genome_info,config):

    fasta_list=set(list_genomes(folder,config['fasta_extension']))


    assert genome_info.index.is_unique, "genome table index are not unique"


    missing_in_table= fasta_list.difference(genome_info.index)
    if len(missing_in_table) >0:
        raise Exception(f"{len(missing_in_table)} files are not in the table"
                        )

    missing_fasta= genome_info.index.difference(fasta_list)
    if len(missing_fasta) >0:
        raise Exception(f"I miss the fasta file for {len( missing_fasta)} genomes."

                    )
    quality_formula = config['quality_score']
    try:
        quality_score= genome_info.eval(quality_formula).isnull()
    except Exception as e:
        raise Exception(f'Cannot calcualte quality score according to formaula "{quality_formula}" \n') from e

    assert quality_score.isnull().any()==False, "Quality score is NA for some genomes"

    print("Database verified.")


def gen_names_for_range_(N,prefix='',start=1):
    """generates a range of IDS with leading zeros so sorting will be ok"""
    n_leading_zeros= len(str(N))
    format_int=prefix+'{:0'+str(n_leading_zeros)+'d}'
    return [format_int.format(i) for i in range(start,N+start)]



def append_to_db(genome_info_table,folder,genome_info_table2,folder2,config):

    genome_info = pd.read_csv(genome_info_table,index_col=0,sep='\t')
    genome_info_new = pd.read_csv(genome_info_table2,index_col=0,sep='\t')



    validate_db(folder,genome_info,config)
    try:
        validate_db(folder2,genome_info_new,config)
    except Exception as e:
        raise Exception("Error in validation of external database") from e

    # Length or length
    if not (('Length' in genome_info_new.columns) and ('N50' in genome_info_new.columns )):
        print("Didn't found genome information in new genome_info"
              " Calculate genome statistics. This can take some time")

        from .genome_stats import get_many_genome_stats

        genome_stats_new= get_many_genome_stats(list_genomes(folder2,config['fasta_extension']), threads=config.get('threads',8))
        genome_info_new= genome_info_new.join(genome_stats_new)

        genome_info_new.to_csv('genome_info_new_with_genomestats.tsv',sep='\t')


    if ('rename_method' is config) or (config['rename_method']=='None'):
        raise NotADirectoryError("I don't know what to do with this for now")
    elif config['rename_method']=='prefix':
        prefix= config['mag_prefix']


        # verify if last genome in db correspond to N genomes in total.
        Ngenomes_in_db= genome_info.shape[0]
        last_index= genome_info.index.values[-1]
        if not int(last_index.replace(prefix,'')) == Ngenomes_in_db:
            raise Exception("Last genome in DB doesn't correspond to N of genomes. \n"
                            f"Found {last_index}, expect prefix {prefix} and number {Ngenomes_in_db}")

        # verify if new genomes can be appended without increasing the leading zeros
        number_of_digits= len(last_index.replace(prefix,''))
        Max_number = 10**number_of_digits -1

        if Ngenomes_in_db + genome_info_new.shape[0] > Max_number:

            raise Exception("I cannot append the new genome to the existing one"
                            " because the naming doesn't allow to increase the number of genomes:\n"
                            f"{prefix}{Max_number} < {Ngenomes_in_db } + {genome_info_new.shape[0]}"
                            )

        #rename new genomes
        genome_info_new['Original_name'] = genome_info_new.index
        genome_info_new['Genome'] = np.arange(genome_info_new.shape[0])+Ngenomes_in_db+1

        format_int=prefix+'{:0'+str(number_of_digits)+'d}'
        genome_info_new['Genome'] = genome_info_new['Genome'].map(format_int.format)


        Mapping= genome_info_new.Genome.copy()

        genome_info_new.index= genome_info_new.Genome
        genome_info_new.drop('Genome',axis=1, inplace=True)


    else:
        raise NotImplementedError("config['rename_method'] should be one of 'None', 'prefix'")




    #copy genomes
    print("Copy fasta files into database. This can take some time")
    fasta_extension= config['fasta_extension']
    for old_name,new_name in tqdm(Mapping.iteritems(),total=Mapping.shape[0]):

        shutil.copy(f"{folder2}/{old_name}{fasta_extension}",
                    f"{folder1}/{new_name}{fasta_extension}"
                   )

    genome_info.to_csv(genome_info_table.replace('.tsv','_old.tsv'),sep='\t')
    genome_info= genome_info.append(genome_info_new,sort=False)
    genome_info.to_csv(genome_info_table.replace('.tsv','_new.tsv'),sep='\t')

    validate_db(folder,genome_info,config)
    shutil.move(genome_info_table.replace('.tsv','_new.tsv'),genome_info_table)

    print("Everithing worked fine!")



#config= yaml.safe_load(open('config.yaml'))
