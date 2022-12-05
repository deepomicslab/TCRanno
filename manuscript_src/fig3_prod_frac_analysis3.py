import pandas as pd
import os,sys
from numpy import *

def parse_infile(infile):
    df = pd.read_csv(infile,sep='\t',low_memory=False,usecols=['aminoAcid','frequencyCount (%)','sequenceStatus'])
    df1 = df[df.sequenceStatus=='In']
    prod_frac = sum(df1['frequencyCount (%)'].tolist())
    df2 = df1[df1['frequencyCount (%)']<1e-4]
    low_freq = sum(df2['frequencyCount (%)'].tolist())
    df3 = df1[df1['frequencyCount (%)']>=1e-2]
    high_freq = sum(df3['frequencyCount (%)'].tolist())
    return prod_frac,low_freq,high_freq

if __name__ == '__main__': 
    indir = 'ImmuneAccess_full'
    sample_dir = sys.argv[1]
    out_prefix = sys.argv[2]
    all_files = [x.replace('_1e-4','') for x in os.listdir(sample_dir)]
    med = []
    prod = []
    low = []
    high = []
    for each in all_files:
        infile = indir+'/'+each
        prod_frac, low_freq, high_freq = parse_infile(infile)
        med_freq = prod_frac - low_freq - high_freq
        prod.append(prod_frac)
        low.append(low_freq)
        med.append(med_freq)
        high.append(high_freq)
    savez(out_prefix+'.npz',prod=array(prod),low=array(low),med=array(med),high=array(high))
