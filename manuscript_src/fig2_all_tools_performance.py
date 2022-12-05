# run under TCRMatch directory; replace TCRMatch's orginal IEDB_data.tsv with benchmark/IEDB_data_benchmark_TCRMatch.tsv
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

def compute_nwhamming_matrix(metricb,test_cdr3s,db_cdr3s,ncpu,matrix):
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
        matrix = pw.matrices.identity_nb_distance_matrix
        for i,seqs in enumerate(query_seqs):
            mat = compute_nwhamming_matrix(metricb,seqs,db_cdr3s,num_cores,matrix)
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

def parse_tcr2tcr_results(infile, TEST, REF_DB):
    query = ''
    data=open(infile)
    SCORE={}
    for lines in data:
        if lines.startswith('Index') or lines.startswith('#') or lines.rstrip()=='':
            continue
        content = lines.split('\t')
        rank = content[1]
        cdr3 = content[2]
        if rank == 'Input': 
            query = cdr3
            answer = TEST[query]
            SCORE[query]=zeros(10)
            k = 0
            continue
        else:
            epts = REF_DB[cdr3]
            if len(epts.intersection(answer))>0:
                SCORE[query][k]=1
            k += 1
    data.close()
    counts = zeros(10)      
    for cdr3 in SCORE:
        for i in range(10):
            if SCORE[cdr3][i] == 1:
                for j in range(i,10):
                    counts[j]+=1
                break
    rates = counts/len(TEST)
    return rates

def parse_other_results(infile, TEST, REF_DB):
    query = ''
    data=open(infile)
    SCORE={}
    prev_query=''
    for lines in data:
        content = lines.split()
        query = content[0]
        cdr3 = content[1]
        if query!=prev_query: 
            answer = TEST[query]
            SCORE[query]=zeros(10)
            k = 0
            prev_query=query
            epts = REF_DB[cdr3]
            if len(epts.intersection(answer))>0:
                SCORE[query][k]=1
        else:
            epts = IEDB[cdr3]
            k += 1
            if len(epts.intersection(answer))>0:
                SCORE[query][k]=1
        
    data.close()
    counts = zeros(10)      
    for cdr3 in SCORE:
        for i in range(10):
            if SCORE[cdr3][i] == 1:
                for j in range(i,10):
                    counts[j]+=1
                break
    rates = counts/len(TEST)
    return rates

def plot_topk_hit_rate(fname,list_of_lists,test_name,size,k=10):
    tcrmatch = list_of_lists[0]
    tcrdist = list_of_lists[1]
    levenshtein = list_of_lists[2]
    nwhamming = list_of_lists[3]
    tcr2tcr = list_of_lists[4]
    labels = [1,2,3,4,5,6,7,8,9,10]
    x = arange(1,k+1)
    plt.figure(figsize=(9,9)) #plt.figure(figsize=(10,8))
    plt.plot(x,nwhamming, color='C0', linestyle='dashed', marker='^',label='nwhamming',markersize=10)
    plt.plot(x,levenshtein,color='C1', linestyle='dashed', marker='*',label='levenshtein',markersize=10)
    plt.plot(x,tcrdist,color='C2', linestyle='dashed', marker='+',label='tcrdist',markersize=10)
    plt.plot(x,tcrmatch,'kv--',label='tcrmatch',markersize=10)
    plt.plot(x,tcr2tcr,'ro--',label='tcr2tcr',markersize=10)
    plt.xlabel('K',fontsize=12)
    plt.ylabel('Hit rate',fontsize=12)
    plt.xticks(x, labels)
    plt.title(test_name+' (n='+str(size)+')',fontsize=20)
    plt.legend(prop={'size': 15})
    plt.savefig(fname, dpi=300)
    plt.show()
    plt.close()

if __name__ == '__main__':
    infile = sys.argv[1]
    test_dict = sys.argv[2]
    outprefix = sys.argv[3]
    with open(test_dict,'rb') as a:
        TEST = pickle.load(a)
    encoder1 = model_predict.load_encoder()
    IEDB, IEDB_VDJ, AO_key, AO_map = core_analysis.load_DB(ref_DB='IEDB')
    core_analysis.tcr2tcr(infile, outprefix+'.tcr2tcr.out', encoder1, IEDB, IEDB_VDJ, AO_key, AO_map, header=False, cdr3_aa_col=None, frequency=False, frequency_col=None, count=False, count_col=None, sep=None, k=10, thread=-1)
    tcr2tcr_res = parse_tcr2tcr_results(outprefix+'.tcr2tcr.out',TEST,IEDB)
    run_tcrmatch(infile,outprefix+'.tcrmatch.out')
    tcrmatch_res = parse_other_results(outprefix+'.tcrmatch.out',TEST,IEDB)
    db_cdr3s = sorted(list(IEDB.keys()))
    run_tcrdist('nwhamming',infile,outprefix+'.nw.out',db_cdr3s)
    nw_res = parse_other_results(outprefix+'.nw.out',TEST,IEDB)
    run_tcrdist('tcrdist',infile,outprefix+'.tcrdist.out',db_cdr3s)
    tcrdist_res = parse_other_results(outprefix+'.tcrdist.out',TEST,IEDB)
    run_tcrdist('levenshtein',infile,outprefix+'.lev.out',db_cdr3s)
    lev_res = parse_other_results(outprefix+'.lev.out',TEST,IEDB)
    print(tcr2tcr_res,tcrmatch_res,tcrdist_res,lev_res,nw_res)
    plot_topk_hit_rate(outprefix+'.all_performance.png',[tcrmatch_res,tcrdist_res,lev_res,nw_res,tcr2tcr_res],outprefix,len(TEST))
