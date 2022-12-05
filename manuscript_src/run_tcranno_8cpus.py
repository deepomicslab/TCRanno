from tcranno import *
import sys,os
import argparse

# input repertoire file(s) under indir; specify a output directory as outdir
# this script assumes the input repertoires are .tsv format

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Inputs')
    parser.add_argument('--indir', type=str, required=True, default='', help='input directory where input repertoire files are located')
    parser.add_argument('--outdir', type=str, required=True, default='', help='output directory')
    args = parser.parse_args()
    indir = args.indir
    outdir = args.outdir
    all_files = [x for x in os.listdir(indir) if x.endswith('tsv')] # take all .tsv files as input, use with caution
    model = model_predict.load_encoder()
    DB, DB_VDJ, AO_map = core_analysis.load_DB(ref_DB='IEDB')
    num_threads = 8
    for i in range(len(all_files)):
        infile=indir+'/'+all_files[i]
        outprefix = outdir+'/'+all_files[i].split('.tsv')[0]
        core_analysis.tcr2tcr(infile=infile, outprefix=outprefix, encoder=model, DB=DB, DB_VDJ=DB_VDJ, AO_map=AO_map, header=True, cdr3_aa_col=0, frequency=True, frequency_col=1, sep='\t', k=10, t=num_threads)
        tcr2tcr_output = outprefix+'_tcr2tcr_output.tsv'
        repertoire_analysis.tcr2ept(infile=tcr2tcr_output,outprefix=outprefix,AO_map=AO_map)
        repertoire_analysis.tcr2ag(infile=tcr2tcr_output,outprefix=outprefix,AO_map=AO_map)
        repertoire_analysis.tcr2org(infile=tcr2tcr_output,outprefix=outprefix,AO_map=AO_map)
        