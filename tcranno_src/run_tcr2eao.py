from tcranno import *
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Inputs')
    parser.add_argument('--infile', type=str, required=True, default='', help='input file')
    parser.add_argument('--is_tcr2tcr', type=str, required=True, default='True', help='whether the input file is a tcr2tcr output')
    parser.add_argument('--outprefix', type=str, required=True, default='', help='prefix of the tcr2ept output file')
    parser.add_argument('--anno_type', type=str, default='all',help='annotation type, one of tcr2ept, tcr2ag, tcr2org, all')
    parser.add_argument('--limit', type=float, default=1e-4,help='filter clonotypes with fraction lower than the specified limit')
    parser.add_argument('--sep', type=str, default='\s+',help='delimitor') # only specify when using non-whitespace delimitor
    parser.add_argument('--cdr3_aa_col', type=str, default=None,help='cdr3_aa column name')
    parser.add_argument('--frequency_col', type=str, default=None,help='frequency column name')
    parser.add_argument('--count_col', type=str, default=None,help='count column name')
    parser.add_argument('--k', type=int, default=10,help='output top K epitopes/antigens/organisms with highest TCR fraction')
    parser.add_argument('--t', type=int, default=1,help='number of threads')
    parser.add_argument('--ref_DB', type=str, default='IEDB',help='path to reference database (.pkl)')
    parser.add_argument('--ref_AO_map', type=str, default=None,help='path to reference epitope-antigen-organism-map-dict (.pkl)')
    parser.add_argument('--model', type=str, default=None,help='path to model for embedding generation')
    args = parser.parse_args()
    infile = args.infile
    is_tcr2tcr = args.is_tcr2tcr
    outprefix = args.outprefix
    limit = args.limit
    ref_DB = args.ref_DB
    ref_AO_map = args.ref_AO_map
    DB, DB_VDJ, AO_map = core_analysis.load_DB(ref_DB=ref_DB,ref_AO_map=ref_AO_map)
    anno_type = args.anno_type
    k = args.k
    ept = False
    ag = False
    org = False
    if anno_type == 'all':
        ept = True
        ag = True
        org = True
    elif anno_type == 'tcr2ept':
        ept = True
    elif anno_type == 'tcr2ag':
        ag = True
    elif anno_type == 'tcr2org':
        org = True
    if is_tcr2tcr.startswith('T') or is_tcr2tcr.startswith('t'): #True, true
        if ept:
            repertoire_analysis.tcr2ept(infile, outprefix, AO_map=AO_map, k=k, limit=limit)
        if ag:
            repertoire_analysis.tcr2ag(infile, outprefix, AO_map=AO_map, k=k, limit=limit)
        if org:
            repertoire_analysis.tcr2org(infile, outprefix, AO_map=AO_map, k=k, limit=limit)
    elif is_tcr2tcr.startswith('F') or is_tcr2tcr.startswith('f'): #False, false
        sep = args.sep
        cdr3_aa_col = args.cdr3_aa_col
        frequency_col = args.frequency_col
        count_col = args.count_col
        t = args.t
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
        tcr2tcr_out = outprefix+'_tcr2tcr_output.tsv'
        model_path = args.model
        if model_path is None:
            encoder1 = model_predict.load_encoder()
        else:
            encoder1 = model_predict.load_encoder(model_path)
        core_analysis.tcr2tcr(infile, outprefix, encoder1, DB, DB_VDJ, AO_map, header=header, cdr3_aa_col=cdr3_aa_col, frequency=frequency, frequency_col=frequency_col, count=count, count_col=count_col, sep=sep, k=1, thread=t, limit=limit)
        if ept:
            repertoire_analysis.tcr2ept(tcr2tcr_out, outprefix, AO_map=AO_map, k=k, limit=limit)
        if ag:
            repertoire_analysis.tcr2ag(tcr2tcr_out, outprefix, AO_map=AO_map, k=k, limit=limit)
        if org:
            repertoire_analysis.tcr2org(tcr2tcr_out, outprefix, DB=DB, AO_map=AO_map, k=k, limit=limit)
    else:
        print('Unrecognized is_tcr2tcr parameter:',is_tcr2tcr)
