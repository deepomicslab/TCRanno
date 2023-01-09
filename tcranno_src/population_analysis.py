import pandas as pd
import os, sys
from numpy import *

def get_topk_output(indict, f, limit, k):
    j = 0
    counted = []
    output = []
    while j < k and len(indict)>0:
        m = max(indict.values())
        if m < limit or m*100/f < 0.1: # only output >= 0.1%
            break
        for keys in indict:
            if indict[keys]==m:
                j += 1
                counted.append(keys)
                out = keys+'('+str(round(m,5))+', '+str(round(m*100/f,2))+'%)'
                output.append(out)
        for each in counted:
            indict.pop(each)
        counted = []
    return output #('; ').join(output)

def population_stats(all_files, ftype, limit=1e-4, k=10):
    n = len(all_files)
    if ftype not in ['tcr2ept','tcr2ag','tcr2org']:
        print('Unrecognized feature type: must be one of {tcr2ept,tcr2ag,tcr2org}');exit(0)
    if ftype == 'tcr2ept':
        CM_name = 'CM:Top_Epitopes (mean_sum_of_frequency)'
        PM_name = 'PM:Top_Epitopes (mean_sum_of_frequency)'
    elif ftype == 'tcr2ag':
        CM_name = 'CM:Top_Antigens (mean_sum_of_frequency)'
        PM_name = 'PM:Top_Antigens (mean_sum_of_frequency)'
    else:
        CM_name = 'CM:Top_Organisms (mean_sum_of_frequency)'
        PM_name = 'PM:Top_Organisms (mean_sum_of_frequency)'
    comp_match = {}
    part_match = {}
    cm_freqs = []
    pm_freqs = []
    cm_count = 0
    for infile in all_files:
        #print(infile)
        df = pd.read_csv(infile,sep='\t')
        colnames = df.columns.tolist()
        #print(colnames)
        col1 = df[colnames[1]].tolist()
        col2 = df[colnames[2]].tolist()
        col3 = df[colnames[3]].tolist()
        cm_freq = float(df.iloc[0]['Complete_Match (CM)'])
        pm_freq = float(df.iloc[0]['Predicted_Match (PM)'])
        cm_freqs.append(cm_freq)
        pm_freqs.append(pm_freq)
        cm_top = [x for x in col1[4:] if x!='-']
        pm_top = [x for x in col2[4:] if x!='-']
        if cm_freq > 0:#len(cm_top) >= 1:
            cm_count+=1
            for t in cm_top:
                t = t.strip(' ')
                if t.count('(')==0:
                    print(t)
                idx = t.rindex('(')
                item = t[:idx].strip(' ')
                frac = t[idx+1:]
                f = float(frac.split(',')[0])
                if item not in comp_match:
                    comp_match[item] = f#*cm_freq
                else:
                    comp_match[item] += f#*cm_freq
        if len(pm_top)<=1:
            continue
        for t in pm_top:
            t = t.strip(' ')
            idx = t.rindex('(')
            item = t[:idx].strip(' ')
            frac = t[idx+1:]
            f = float(frac.split(',')[0])
            if item not in part_match:
                part_match[item] = f#*pm_freq
            else:
                part_match[item] += f#*pm_freq
    mean_cm_freq = mean(cm_freqs)
    mean_pm_freq = mean(pm_freqs)
    for items in comp_match:
        comp_match[items]/=n #(n*mean_cm_freq)
    for items in part_match:
        part_match[items]/=n #(n*mean_pm_freq)
    cm = get_topk_output(comp_match, mean_cm_freq, limit, k)
    pm = get_topk_output(part_match, mean_pm_freq, limit, k)
    maxk = max(len(cm),len(pm))
    if len(cm) < maxk:
        for j in range(maxk - len(cm)):
            cm.append('.')
    cm.append('-'*40)
    cm.append(mean_cm_freq)
    if len(pm) < maxk:
        for j in range(maxk - len(pm)):
            pm.append('.')
    pm.append('-'*40)
    pm.append(mean_pm_freq)
    #lst = [['Complete Match',mean_cm_freq,cm_out],['Partial Match',mean_pm_freq,pm_out]]
    ranks = [i+1 for i in range(k)]
    ranks.append('---')
    ranks.append('Total')
    #df2 = pd.DataFrame(lst, columns =['Match Type', 'Mean_Sum_of_Frequency', df.columns[-1]])
    df2 = pd.DataFrame(list(zip(ranks, cm, pm)),columns =['Rank', CM_name, PM_name])
    print("total no. of samples:",n)
    print("no. of samples with complete matches:",cm_count)
    return df2
    
def get_files(indirs, file_type): # {'tcr2ept','tcr2ag','tcr2org'}
    all_files=[]
    for indir in indirs:
        all_files += [indir+'/'+x for x in os.listdir(indir) if x.count(file_type)>0]
    return all_files
    
if __name__=='__main__':
    indirs = sys.argv[1].split(',')
    ftype = sys.argv[2]
    topk = int(sys.argv[3])
    all_files = get_files(indirs, ftype)
    df = population_stats(all_files, ftype, k=topk)
    if len(sys.argv) > 4:
        outfile = sys.argv[4]
        df.to_csv(outfile,sep='\t',index=False)
    else:
        #pd.display(df.to_string(index=False))
        print(df.to_markdown(index=False))
