#!/usr/bin/env python3
from snakemake.shell import shell
from common.io import simplify_path
import os
from glob import glob

from multiprocessing.pool import Pool


def bbsketch(fasta,command):
    name0=simplify_path(fasta)
    #print(command.format(fasta=fasta,name0=name0))
    shell(command)




def sketch(input,sketch=None,translate=True,k=None,threads=1,
           extra=""):

    if len(input)==1 & os.path.isdir(input[0]):
        fasta_files=glob(f"{input[0]}/*")
    else:
        fasta_files= input

    if sketch is None:
        sketch=f"bbsketch_{'AA' if translate else 'nt'}.sketch.gz"

        print(f"Save sketch in {sketch}")

    if os.path.exists(sketch):
        raise IOError("Output sketch already exists!")



    command= "bbsketch.sh " \
             "in={fasta} name0={name0} overwrite=t "\
             f"translate={'t' if translate else 'f'} " \

    output_file= "stdout.sketch"
    if sketch.endswith('.gz'):
        output_file+='.gz'
    command +=f"out={output_file} "

    if k is not None:
        command+=f"k={k} "
    command+=f"{extra} >> {sketch} 2> /dev/null "




    pool = Pool(threads)

    pool.starmap(bbsketch,[(f,command) for f in fasta_files])





if __name__ == "__main__":

    import argparse

    p = argparse.ArgumentParser()
    p.add_argument("-o","--sketch", default=None)
    p.add_argument("-k")
    p.add_argument("--translate",type=bool,default=True)
    p.add_argument(dest='input',nargs='+')
    p.add_argument("--threads", required = False, default = 4, type = int)
    args = vars(p.parse_args())
    sketch(**args)
