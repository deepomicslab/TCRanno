from tcranno import *
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Inputs')
    parser.add_argument('--infile', type=str, required=True, default='', help='input file')
    parser.add_argument('--outprefix', type=str, required=True, default='', help='prefix of the tcr2tcr output file')
    parser.add_argument('--sep', type=str, default='\s+',help='delimitor') # only specify when using non-whitespace delimitor
    parser.add_argument('--cdr3_aa_col', type=str, default=None,help='cdr3_aa column name')
    parser.add_argument('--frequency_col', type=str, default=None,help='frequency column name')
    parser.add_argument('--count_col', type=str, default=None,help='count column name')
    parser.add_argument('--k', type=int, default=10,help='output top K closest TCR sequences for each query')
    parser.add_argument('--t', type=int, default=1,help='number of threads')
    parser.add_argument('--limit', type=str, default=None,help='frequency filter')
    parser.add_argument('--ref_DB', type=str, default='IEDB',help='path to reference database (.pkl)')
    parser.add_argument('--ref_AO_map', type=str, default=None,help='path to reference epitope-antigen-organism-map-dict (.pkl)')
    parser.add_argument('--model', type=str, default=None,help='path to model for embedding generation')
    args = parser.parse_args()
    infile = args.infile
    outprefix = args.outprefix
    sep = args.sep
    cdr3_aa_col = args.cdr3_aa_col
    frequency_col = args.frequency_col
    count_col = args.count_col
    k = args.k
    t = args.t
    limit = args.limit
    if cdr3_aa_col is None:
        header = False
    else:
        header = True
    if frequency_col is None:
        frequency = False
    else:
        frequency = True
    if count_col is None:
        count = False
    else:
        count = True
    if not limit is None:
        limit = float(limit)
    ref_DB = args.ref_DB
    ref_AO_map = args.ref_AO_map
    model_path = args.model
    if model_path is None:
        encoder1 = model_predict.load_encoder()
    else:
        encoder1 = model_predict.load_encoder(model_path)
    DB, DB_VDJ, AO_map = core_analysis.load_DB(ref_DB=ref_DB,ref_AO_map=ref_AO_map)
    core_analysis.tcr2tcr(infile, outprefix, encoder1, DB, DB_VDJ, AO_map, header=header, cdr3_aa_col=cdr3_aa_col, frequency=frequency, frequency_col=frequency_col, count=count, count_col=count_col, sep=sep, k=k, thread=t, limit=limit)
