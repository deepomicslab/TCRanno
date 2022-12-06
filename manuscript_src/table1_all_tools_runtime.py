# random generation of test repertoire of size N
import sys, time, pickle, math, subprocess
import pandas as pd
from matplotlib import pyplot as plt
from numpy import *
from tcranno import *

from tcrdist.rep_funcs import _pw
import pwseqdist as pw
from pwseqdist.matrices import seq2vec
from pwseqdist.nb_metrics import nb_tcrdist

def compute_tcrdist_matrix(metricb,test_cdr3s,db_cdr3s,ncpu,matrix):
    i=0
    while i+100 <= len(test_cdr3s):
        cdr3 = test_cdr3s[i:i+100]
        row = _pw(metric=metricb,seqs1=cdr3,seqs2=db_cdr3s,ncpus=ncpu,use_numba=True,distance_matrix=matrix)
        if i==0:
            dmat = row
        else:
            dmat = vstack((dmat,row))
        i+=100
    if len(test_cdr3s) > i:
        cdr3=test_cdr3s[i:]
        row = _pw(metric=metricb,seqs1=cdr3,seqs2=db_cdr3s,ncpus=ncpu,use_numba=True,distance_matrix=matrix)
        dmat = vstack((dmat,row))
    return dmat

def compute_nwhamming_matrix(metricb,test_cdr3s,db_cdr3s,ncpu):
    i=0
    while i+100 <= len(test_cdr3s):
        cdr3 = test_cdr3s[i:i+100]
        row = _pw(metric=metricb,seqs1=cdr3,seqs2=db_cdr3s,ncpus=ncpu)
        if i==0:
            dmat = row
        else:
            dmat = vstack((dmat,row))
        i+=100
    if len(test_cdr3s) > i:
        cdr3=test_cdr3s[i:]
        row = _pw(metric=metricb,seqs1=cdr3,seqs2=db_cdr3s,ncpus=ncpu)
        dmat = vstack((dmat,row))
    return dmat
    
def run_tcrdist(method,infile,outfile,db_cdr3s):
    num_cores = 8
    with open(infile,'r') as a:
        lines = a.readlines()
    test_cdr3s = [x.rstrip() for x in lines]
    query_seqs=[test_cdr3s[i*1000:(i+1)*1000] for i in range(int(math.ceil(len(test_cdr3s)/1000)))]
    if method == 'nwhamming':
        metricb = pw.metrics.nw_hamming_metric
        for i,seqs in enumerate(query_seqs):
            mat = compute_nwhamming_matrix(metricb,seqs,db_cdr3s,num_cores)
            if i==0:
                dmat = mat
            else:
                dmat = vstack((dmat,mat))
    elif method == 'tcrdist':
        metricb = pw.metrics.nb_vector_tcrdist
        matrix=pw.matrices.tcr_nb_distance_matrix
        for i,seqs in enumerate(query_seqs):
            mat = compute_tcrdist_matrix(metricb,seqs,db_cdr3s,num_cores,matrix)
            if i==0:
                dmat = mat
            else:
                dmat = vstack((dmat,mat))
    elif method == 'levenshtein':
        metricb = pw.metrics.nb_vector_editdistance
        matrix = pw.matrices.identity_nb_distance_matrix
        for i,seqs in enumerate(query_seqs):
            mat = compute_tcrdist_matrix(metricb,seqs,db_cdr3s,num_cores,matrix)
            if i==0:
                dmat = mat
            else:
                dmat = vstack((dmat,mat))
    pipeline=open(outfile,'w')
    for i in range(len(test_cdr3s)):
        row = dmat[i]
        sorted_index_array = argsort(row) # smallest rank first
        topk = sorted_index_array[:10]
        candidates = array(db_cdr3s)[topk]
        scores = row[topk]
        for j in range(10):
            pipeline.write(test_cdr3s[i]+'\t'+candidates[j]+'\t'+str(scores[j])+'\n')
    pipeline.close()
    
def run_tcrmatch(infile,outfile):
    cmd = './tcrmatch -i '+infile+' -t 8 -s 0.84 > '+outfile
    print(cmd)
    subprocess.call(cmd, shell=True)
    tcrmatch_output = []
    with open(outfile, "r") as infile:
        prev_query = ''
        records = {}
        count = 0
        for i,line in enumerate(infile):
            line_content = line.rstrip().split(" ")
            seq = line_content[0]
            match = line_content[1]
            score = float(line_content[2])
            if i!=0 and seq!=prev_query:
                count += 1
                top = sorted(records.items(),reverse=True,key=lambda x:x[1])
                for j in range(min(10,len(records))):
                    tcrmatch_output.append((prev_query,top[j][0],str(top[j][1])))
                records = {}
                records[match]=score
                prev_query=seq
            else:
                records[match]=score
                prev_query = seq
        count += 1
    #print(prev_query,count)
        top = sorted(records.items(),reverse=True,key=lambda x:x[1])
        for j in range(min(10,len(records))):
            tcrmatch_output.append((seq,top[j][0],str(top[j][1])))
    pipeline=open(outfile,'w')     
    for records in tcrmatch_output:
        pipeline.write((' ').join(records)+'\n')
    pipeline.close()

if __name__ == '__main__':
    size = int(sys.argv[1])
    infile = '/mnt/d/TCR_data/healthy587_first10000tcr_pools.tsv'
    df = pd.read_csv(infile,sep='\t')
    df2 = df.sample(n = size)
    testfile = 'random_repertoire_size_'+str(size)+'.csv'
    tcr2tcr_out = 'runtime.size.'+str(size)+'_tcr2tcr_output.tsv'
    tcrmatch_out = 'tcrmatch.runtime.size.'+str(size)+'.tsv'
    tcrdist_out = 'tcrdist.runtime.size.'+str(size)+'.tsv'
    nw_out = 'nwhamming.runtime.size.'+str(size)+'.tsv'
    lev_out = 'levenshtein.runtime.size.'+str(size)+'.tsv'
    df2.to_csv(testfile,columns=['cdr3_aa'],index=False,header=False)
    
    start = time.time()
    encoder1 = model_predict.load_encoder()
    IEDB, IEDB_VDJ, AO_map = core_analysis.load_DB(ref_DB='IEDB')
    core_analysis.tcr2tcr(testfile, 'runtime.size.'+str(size)+', encoder1, IEDB, IEDB_VDJ, AO_map, header=False, cdr3_aa_col=None, frequency=False, frequency_col=None, count=False, count_col=None, sep=None, k=10, thread=-1)
    end = time.time()
    elapsed = end - start
    print('TCRanno Runtime: %.2f seconds.' % elapsed)
    
    start = time.time()
    run_tcrmatch(testfile,tcrmatch_out)
    end = time.time()
    elapsed = end - start
    print('TCRMatch Runtime: %.2f seconds.' % elapsed)
    
    db_cdr3s = sorted(list(IEDB.keys()))
    
    start = time.time()
    run_tcrdist('nwhamming',testfile,nw_out,db_cdr3s)
    end = time.time()
    elapsed = end - start
    print('nwhamming Runtime: %.2f seconds.' % elapsed)
    
    start = time.time()
    run_tcrdist('tcrdist',testfile,tcrdist_out,db_cdr3s)
    end = time.time()
    elapsed = end - start
    print('tcrdist Runtime: %.2f seconds.' % elapsed)
    
    start = time.time()
    run_tcrdist('levenshtein',testfile,lev_out,db_cdr3s)
    end = time.time()
    elapsed = end - start
    print('levenshtein Runtime: %.2f seconds.' % elapsed)

