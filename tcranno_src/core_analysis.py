import pandas as pd
import multiprocessing as mp
import os,sys,pickle
import time
from tcranno.model_predict import *
from tcranno.FindClosestSeq import *
from math import ceil

def load_VJ_map():
    d = os.path.dirname(sys.modules['tcranno'].__file__)
    vj_map = 'data/VJ_map.pkl'
    path = os.path.join(d,vj_map)
    with open(path,'rb') as g: # VJ_map
        print("VJ_map loading finished")
        VJ_map = pickle.load(g)
    return VJ_map

def check_valid_cdr3(cdr3, aa=set(['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V'])):
    if len(cdr3)>=4 and len(cdr3)<=30:
        if set(cdr3)<=aa:
            return True
    return False

def parse_file(infile, header=False, cdr3_aa_col=None, frequency=False, frequency_col=None, count=False, count_col=None, sep='\s+'):
    cdr3s = []
    fractions = []
    if header==False:
        data = open(infile)
        for lines in data:
            lines = lines.rstrip().upper()
            if check_valid_cdr3(lines):
                cdr3s.append(lines)
            else:
                print(lines,'is not a valid CDR3 sequence, please check!')
    else:
        if cdr3_aa_col==None:
            print("No cdr3_aa column given!");exit(0)
        df = pd.read_csv(infile, sep=sep)
        cols = df.columns.tolist()
        if str(cdr3_aa_col).isnumeric():
            cdr3_aa_col = cols[int(cdr3_aa_col)]
        cdr3_aa = df[cdr3_aa_col].tolist()
        if frequency:
            if frequency_col==None:
                print("frequency=True but no frequency column given!");exit(0)
            if str(frequency_col).isnumeric():
                frequency_col = cols[int(frequency_col)]
            fracs = df[frequency_col].tolist()
        elif count:
            if count_col==None:
                print("count=True but no count column given!");exit(0)
            if str(count_col).isnumeric():
                count_col = cols[int(count_col)]
            fracs = df[count_col].tolist()
            total = sum(fracs)
            fracs = [x/total for x in fracs]
        else:
            fracs = [1]*len(cdr3_aa)
            total = sum(fracs)
            fracs = [x/total for x in fracs]
        for i,cdr3 in enumerate(cdr3_aa):
            cdr3 = cdr3.upper()
            if check_valid_cdr3(cdr3):
                cdr3s.append(cdr3)
                fractions.append(fracs[i])
            else:
                print(cdr3,'is not a valid CDR3 sequence, please check!')
    return cdr3s, fractions

def load_DB(ref_DB='IEDB', ref_AO_map=None):
    if ref_DB =='IEDB' or ref_DB == 'DB_FULL':
        d = os.path.dirname(sys.modules['tcranno'].__file__)
        resource1 = 'data/'+ref_DB+'.pkl'
        resource2 = 'data/'+ref_DB+'_VDJ.pkl'
        resource3 = 'data/DB_FULL_AG_ORG_key_map.pkl'
        path1 = os.path.join(d,resource1)
        path2 = os.path.join(d,resource2)
        path3 = os.path.join(d,resource3)
        with open(path1,'rb') as a: # VJ_map
            DB = pickle.load(a)
        with open(path2,'rb') as b:
            DB_VDJ = pickle.load(b)
        with open(path3,'rb') as c:
            AO_map = pickle.load(c)
    else:
        with open(ref_DB,'rb') as a:
            DB = pickle.load(a)
        DB_VDJ={}
        VJ_map = load_VJ_map()
        for cdr3 in DB:
            DB_VDJ[cdr3] = FindVDJ(cdr3,VJ_map)
        if ref_AO_map==None:
            d = os.path.dirname(sys.modules['tcranno'].__file__)
            resource3 = 'data/DB_FULL_AG_ORG_key_map.pkl'
            path3 = os.path.join(d,resource3)
        else:
            path3 = ref_AO_map
        with open(path3,'rb') as b:
            AO_map = pickle.load(b)
    return DB, DB_VDJ, AO_map
    
def tcr2tcr(infile, outprefix, encoder, DB, DB_VDJ, AO_map, header=False, cdr3_aa_col=None, frequency=False, frequency_col=None, count=False, count_col=None, sep='\s+', k=10, thread=1, batch=200):
    #start = time.time()
    outfile = outprefix+'_tcr2tcr_output.tsv'
    query_cdr3s, fracs = parse_file(infile,header,cdr3_aa_col,frequency,frequency_col,count,count_col,sep)
    if len(fracs)==len(query_cdr3s):
        fractions = fracs
    else:
        f = 1.0/len(query_cdr3s)
        fractions = [f]*len(query_cdr3s) #None
    query_latent = get_norm_latent(query_cdr3s,encoder)
    db_cdr3s = sorted(list(DB.keys()))
    db_latent = get_norm_latent(db_cdr3s,encoder)
    VJ_map = load_VJ_map()
    output = 'Index\tRank\tCDR3_sequence\tFrequency\tMatched_Epitope(s) [Matched_Antigen->Matched_Organism]\tWeighted_Score\tc1|c2|c3\tw1|w2|w3\n'
    #batch = 200
    t = int(ceil(len(query_cdr3s)/batch))
    if thread == -1:
        thread = len(os.sched_getaffinity(0)) # no. of cpus available; not the total no. of physical cpus
    if thread > 1 and t>1:
        print("No. of CPUs available:",thread)
        pool = mp.Pool(thread)
        items = [(query_cdr3s[i*batch:(i+1)*batch], i*batch+1, fractions[i*batch:(i+1)*batch], db_cdr3s, db_latent, query_latent[i*batch:(i+1)*batch], DB, DB_VDJ, AO_map, VJ_map, k) for i in range(t)]
        outputs = pool.starmap_async(FindClosestSeq_batch, items)
        for o in outputs.get():
            output+=o
    else:
        output += FindClosestSeq_batch(query_cdr3s, 1, fractions, db_cdr3s, db_latent, query_latent, DB, DB_VDJ, AO_map, VJ_map, k)
    complete_matches = output.count('Complete_Match')
    heads = '### tcr2tcr output format\n## total input sequences: '+str(len(query_cdr3s))+' completely matched sequences: '+str(complete_matches)+'\n\n'
    pipeline = open(outfile,'w')
    pipeline.write(heads+output)
    pipeline.close()
    #end = time.time()
    #elapsed = end - start
    #print('Elapsed time: %.2f seconds.' % elapsed)

def tcr2tcr_df(infile, encoder, DB, DB_VDJ, AO_map, header=False, cdr3_aa_col=None, frequency=False, frequency_col=None, count=False, count_col=None, sep='\s+', k=10, thread=1, batch=200):
    query_cdr3s, fracs = parse_file(infile,header,cdr3_aa_col,frequency,frequency_col,count,count_col,sep)
    if len(fracs)==len(query_cdr3s):
        fractions = fracs
    else:
        f = 1.0/len(query_cdr3s)
        fractions = [f]*len(query_cdr3s) #None
    query_latent = get_norm_latent(query_cdr3s,encoder)
    db_cdr3s = sorted(list(DB.keys()))
    db_latent = get_norm_latent(db_cdr3s,encoder) 
    VJ_map = load_VJ_map()
    t = int(ceil(len(query_cdr3s)/batch))
    if thread == -1:
        thread = len(os.sched_getaffinity(0)) # no. of cpus available; not the total no. of physical cpus
    if thread > 1 and t>1:
        print("No. of CPUs available:",thread)
        pool = mp.Pool(thread)
        items = [(query_cdr3s[i*batch:(i+1)*batch], i*batch+1, fractions[i*batch:(i+1)*batch], db_cdr3s, db_latent, query_latent[i*batch:(i+1)*batch], DB, DB_VDJ, AO_map, VJ_map, k) for i in range(t)]
        outputs = pool.starmap_async(FindClosestSeq_batch_lst, items)
        output = []
        for ops in outputs.get():
            for o in ops:
                output.append(o)
    else:
        output = FindClosestSeq_batch_lst(query_cdr3s, 1, fractions, db_cdr3s, db_latent, query_latent, DB, DB_VDJ, AO_map, VJ_map, k)
    colname = ['Index','Rank','CDR3_sequence','Frequency','Matched_Epitope(s) [Matched_Antigen->Matched_Organism]','Weighted_Score','c1|c2|c3','w1|w2|w3']
    df = pd.DataFrame(output,columns=colname)
    return df

