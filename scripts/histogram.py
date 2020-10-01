#!/usr/bin/env python

import pandas as pd
import numpy as np
import os,sys

if len(sys.argv)>1:

    file=sys.argv[1]
else:
    file='tables/bindash_dists.tsv'

out_file= file.replace('.tsv','.parquet')




def simplify_index(index):
    "assumes single index are path of files, removes extesnion and dirname"

    path = index[0]
    filename, extension = os.path.splitext(path)

    if extension == ".gz":
        extension = os.path.splitext(filename)[-1] + extension

    dir_name = os.path.dirname(path) + "/"

    return index.str.replace(extension, "").str.replace(dir_name, "")


if os.path.exists(out_file):

    M= pd.read_parquet(out_file)
else:

        

    data=dict(Genome1=[],Genome2=[],Distance=[],P=[],Matching=[],Total=[])

    with open(file) as f:
        for line in f:

            g1,g2,dist,p,fract = line.strip().split()
            
            matching,total= fract.split('/')

            data['Genome1'].append(g1)
            data['Genome2'].append(g2)
            data['Distance'].append(float(dist))
            data['P'].append(float(p))
            data['Matching'].append(int(matching))
            data['Total'].append(int(total))

    M=pd.DataFrame(data)
    del data

    M['Genome1'] = simplify_index(M.Genome1)
    M['Genome2'] = simplify_index(M.Genome2)

    M.to_parquet(out_file,index=False)
    print('saved to parquet')

Msub=M.query('Distance<0.1')


counts, bins=np.histogram(Msub.Distance.values,1000)
hist=pd.Series(counts,index=bins[1:],name='Counts')
hist.index.name='Distance'


hist.to_csv('histogram_dist.tsv',sep='\t',header=True)


counts, bins=np.histogram(np.log10(M.Distance.values),1000,range=(-6,1))
hist=pd.Series(counts,index=bins[1:],name='Counts')

hist.index.name="log10dist"


hist.to_csv('histogram_logdist.tsv',sep='\t',header=True)
