import pandas as pd
import os,sys
from numpy import *

def parse_infile(infile): ##CDR3_aa  Count   Frequency
    df = pd.read_csv(infile,sep='\t',low_memory=False)
    prod_frac = sum(df['Frequency'].tolist())
    #nonprod_frac = 1-prod_frac
    df2 = df[df['Frequency']<1e-4]
    low_freq = sum(df2['Frequency'].tolist())
    df3 = df[df['Frequency']>=1e-2]
    high_freq = sum(df3['Frequency'].tolist())
    return prod_frac,low_freq,high_freq

if __name__ == '__main__':
    indir = sys.argv[1] # SLE_processed_data
    out_prefix = sys.argv[2]
    all_files = os.listdir(indir)
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
