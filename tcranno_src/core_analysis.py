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

def clonotype_stats(freqs):
    prod = 0
    low=0
    med=0
    high=0
    for f in freqs:
        if f < 1e-4:
            low+=f
        elif f < 1e-2:
            med+=f
        else:
            high+=f
        prod += f
    return prod,low,med,high

def parse_file(infile, header=False, cdr3_aa_col=None, frequency=False, frequency_col=None, count=False, count_col=None, sep='\s+', limit=None, perform_stats=False):
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
        if perform_stats:
            prod, low, med, high = clonotype_stats(fracs)
        if limit is None:
            for i,cdr3 in enumerate(cdr3_aa):
                cdr3 = cdr3.upper()
                if check_valid_cdr3(cdr3):
                    cdr3s.append(cdr3)
                    fractions.append(fracs[i])
                else:
                    print(cdr3,'is not a valid CDR3 sequence, please check!')
        else:
            filtered = 0
            for i,cdr3 in enumerate(cdr3_aa):
                if fracs[i]<limit:
                    filtered +=1
                    continue
                cdr3 = cdr3.upper()
                if check_valid_cdr3(cdr3):
                    cdr3s.append(cdr3)
                    fractions.append(fracs[i])
                else:
                    print(cdr3,'is not a valid CDR3 sequence, please check!')
            print(filtered,'sequences have been filitered by the frequency filter:',limit)
    if perform_stats:
        return cdr3s, fractions, prod, low, med, high
    else:
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

def tcr2tcr(infile, outprefix, encoder, DB, DB_VDJ, AO_map, header=False, cdr3_aa_col=None, frequency=False, frequency_col=None, count=False, count_col=None, sep='\s+', k=10, thread=1, batch=200, perform_stats=True, limit=None):
    #start = time.time()
    outfile = outprefix+'_tcr2tcr_output.tsv'
    if perform_stats and header:
        query_cdr3s, fracs, prod, low, med, high = parse_file(infile,header,cdr3_aa_col,frequency,frequency_col,count,count_col,sep,limit,perform_stats)
    else:
        query_cdr3s, fracs = parse_file(infile,header,cdr3_aa_col,frequency,frequency_col,count,count_col,sep,limit)
    if len(fracs)==len(query_cdr3s):
        fractions = fracs
    else:
        f = 1.0/len(query_cdr3s)
        fractions = [f]*len(query_cdr3s) #None
    query_latent = get_norm_latent(query_cdr3s,encoder)
    db_cdr3s = sorted(list(DB.keys()))
    db_latent = get_norm_latent(db_cdr3s,encoder)
    VJ_map = load_VJ_map()
    output = 'Index\tRank\tCDR3_sequence\tFrequency\tMatched_Epitope(s) [Matched_Antigen->Matched_Organism]\tWeighted_Score\ts_T|s_D|s_E\tw_T|w_D|w_E\n'
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
    if limit is None:
        limit = 'None'
    else:
        limit = str(limit)
    if perform_stats:
        complete_matches = output.count('Complete_Match')
        heads = '### tcr2tcr output format\n## total input sequences: '+str(len(query_cdr3s))+', completely matched sequences: '+str(complete_matches)+'\n# frequency filter: '+str(limit)+'; productive fraction: '+str(round(prod,4))+'; low-freq fraction: '+str(round(low,4))+'; medium-freq fraction: '+str(round(med,4))+'; high-freq fraction: '+str(round(high,4))+'; effective fraction: '+str(round(med+high,4))+'\n'
    else:
        heads = '### tcr2tcr output format\n## frequency filter: '+str(limit)+'\n\n'
    pipeline = open(outfile,'w')
    pipeline.write(heads+output)
    pipeline.close()
    #end = time.time()
    #elapsed = end - start
    #print('Elapsed time: %.2f seconds.' % elapsed)

def tcr2tcr_df(infile, encoder, DB, DB_VDJ, AO_map, header=False, cdr3_aa_col=None, frequency=False, frequency_col=None, count=False, count_col=None, sep='\s+', k=10, thread=1, batch=200, limit=None):
    print('frequency filter:',limit)
    query_cdr3s, fracs = parse_file(infile,header,cdr3_aa_col,frequency,frequency_col,count,count_col,sep,limit)
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
    colname = ['Index','Rank','CDR3_sequence','Frequency','Matched_Epitope(s) [Matched_Antigen->Matched_Organism]','Weighted_Score','s_T|s_D|s_E','w_T|w_D|w_E']
    df = pd.DataFrame(output,columns=colname)
    return df

