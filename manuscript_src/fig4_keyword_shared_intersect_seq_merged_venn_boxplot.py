import pandas as pd
from numpy import *
import random as rd
import matplotlib
from matplotlib import pyplot as plt
import os,sys
from scipy.stats import mannwhitneyu
from statannot import add_stat_annotation
from matplotlib_venn import venn2,venn2_circles
import seaborn as sns

def parse_infile(infile, keyword, match_type='CM'):
    df = pd.read_csv(infile,sep='\t',low_memory=False,skiprows=3)
    ept_ag_org_column = df.columns[4]
    freq_column = df.columns[3]
    cdr3_column = df.columns[2]
    if match_type=='CM':
        m = df[df['Index'].str.contains('-0') & df[ept_ag_org_column].str.contains(keyword)]
    else:
        m = df[df['Index'].str.endswith('-1') & df[ept_ag_org_column].str.contains(keyword)]
    cdr3s = m[cdr3_column].tolist()
    fracs = m[freq_column].tolist()
    FRACS = {}
    for i,cdr3 in enumerate(cdr3s):
        if cdr3 not in FRACS:
            FRACS[cdr3]=[]
        FRACS[cdr3].append(fracs[i])
    return FRACS
    
def parse_infile2(infile, selected_cdr3s, match_type='CM'):
    df = pd.read_csv(infile,sep='\t',low_memory=False,skiprows=3)
    ept_ag_org_column = df.columns[4]
    freq_column = df.columns[3]
    cdr3_column = df.columns[2]
    if match_type=='CM':
        m = df[df['Index'].str.contains('-0') & df[cdr3_column].isin(selected_cdr3s)]
    else:
        m = df[df['Index'].str.endswith('-1') & df[cdr3_column].isin(selected_cdr3s)]
    fracs = m[freq_column].tolist()
    return sum(fracs)
    
def get_files(indirs, file_type): # {'tcr2ept','tcr2ag','tcr2org', 'tcr2tcr'}
    all_files=[]
    for indir in indirs:
        all_files += [indir+'/'+x for x in os.listdir(indir) if x.count(file_type)>0]
    return all_files

def get_seqs_for_venn(all_files, keyword, match_type='CM'):
    SEQS = {}
    for infile in all_files:
        FRACS = parse_infile(infile, keyword, match_type)
        for keys in FRACS:
            if keys not in SEQS:
                SEQS[keys]=[]
            SEQS[keys]+=FRACS[keys]
    return SEQS
    
def get_seqs_for_boxplot(all_files, selected_cdr3s, match_type='CM'):
    SEQS = []
    for infile in all_files:
        SEQS.append(parse_infile2(infile, selected_cdr3s, match_type))
    return SEQS
    
def plot_venn(pos,neg,intersect,keyword,labels,match_type='CM'):
    font = {'family' : 'sans',
        'weight' : 'bold',
        'size'   : 22}
    matplotlib.rc('font', **font)
    plt.figure(figsize=(5,4))
    #colors = ['royalblue','orangered','forestgreen','gold']
    if match_type=='CM':
        venn2(subsets = (pos, neg, intersect), set_labels = labels, set_colors=('royalblue','orangered'),alpha=1.0)
    else:
        venn2(subsets = (pos, neg, intersect), set_labels = labels, set_colors=('forestgreen','gold'),alpha=1.0)
    venn2_circles(subsets=(pos, neg, intersect));plt.title(keyword)
    plt.savefig('venn_'+keyword+'_'+match_type+'.png',dpi=300)
    plt.close()
    
def get_df_p(frac1,frac2,frac3,frac4,keyword,match_type='CM'):
    labels = ['CMV+ only','Shared (CMV+)','Shared (CMV-)','CMV- only']
    df = pd.DataFrame(list(zip(frac1+frac2+frac3+frac4,[labels[0]]*len(frac1)+[labels[1]]*len(frac2)+[labels[2]]*len(frac3)+[labels[3]]*len(frac4),[keyword]*(len(frac1)+len(frac2)+len(frac3)+len(frac4)))),columns = ['Fraction','Group','Keyword'])
    U, p = mannwhitneyu(frac2, frac3)
    return df,p
    
def plot_boxplot(df,ps,keywords,match_type='CM'):
    matplotlib.style.use('default')
    plt.figure(figsize=(8,4))
    colors = ['royalblue','skyblue','lightsalmon','orangered'];sns.set_palette(array([matplotlib.colors.to_rgba(x) for x in colors]))
    #colors = ['forestgreen','yellowgreen','lemonchiffon','gold'];sns.set_palette(array([matplotlib.colors.to_rgba(x) for x in colors]))
    labels = ['CMV+ only','Shared (CMV+)','Shared (CMV-)','CMV- only']
    ax = sns.boxplot(x=df['Keyword'],y=df['Fraction'],hue=df['Group'],showfliers = False)
    box_pairs = [((k,labels[1]),(k,labels[2])) for k in df.Keyword.unique()]
    add_stat_annotation(ax,data=df,x=df['Keyword'],y=df['Fraction'],hue=df['Group'],box_pairs=box_pairs,perform_stat_test=False,pvalues=ps,text_format='simple',loc='outside', verbose=1)
    ax.set(xlabel=None);ax.legend()
    plt.savefig('merged_boxplot_seq_fraction_color1.'+match_type+'.png',dpi=300)
    plt.close()
    
def get_boxplot_inputs(all_files_cmv_pos,all_files_cmv_neg, keyword):
    label = ('CMV+','CMV-')
    SEQS_POS = get_seqs_for_venn(all_files_cmv_pos, keyword, match_type='CM')
    SEQS_NEG = get_seqs_for_venn(all_files_cmv_neg, keyword, match_type='CM')
    pos_seqs = set(SEQS_POS.keys())
    neg_seqs = set(SEQS_NEG.keys())
    shared = pos_seqs.intersection(neg_seqs)
    only_pos = pos_seqs - shared
    only_neg = neg_seqs - shared
    plot_venn(len(only_pos),len(only_neg),len(shared),keyword,label,match_type='CM')
    POS = get_seqs_for_boxplot(all_files_cmv_pos, only_pos, match_type='CM')
    NEG = get_seqs_for_boxplot(all_files_cmv_neg, only_neg, match_type='CM')
    SPOS = get_seqs_for_boxplot(all_files_cmv_pos, shared, match_type='CM')
    SNEG = get_seqs_for_boxplot(all_files_cmv_neg, shared, match_type='CM')
    cmdf,cmp = get_df_p(POS,SPOS,SNEG,NEG,keyword,match_type='CM') 

    SEQS_POS = get_seqs_for_venn(all_files_cmv_pos, keyword, match_type='PM')
    SEQS_NEG = get_seqs_for_venn(all_files_cmv_neg, keyword, match_type='PM')
    pos_seqs = set(SEQS_POS.keys())
    neg_seqs = set(SEQS_NEG.keys())
    shared = pos_seqs.intersection(neg_seqs)
    only_pos = pos_seqs - shared
    only_neg = neg_seqs - shared
    plot_venn(len(only_pos),len(only_neg),len(shared),keyword,label,match_type='PM')
    POS = get_seqs_for_boxplot(all_files_cmv_pos, only_pos, match_type='PM')
    NEG = get_seqs_for_boxplot(all_files_cmv_neg, only_neg, match_type='PM')
    SPOS = get_seqs_for_boxplot(all_files_cmv_pos, shared, match_type='PM')
    SNEG = get_seqs_for_boxplot(all_files_cmv_neg, shared, match_type='PM')
    pmdf,pmp = get_df_p(POS,SPOS,SNEG,NEG,keyword,match_type='PM') 
    return cmdf, cmp, pmdf, pmp
    
if __name__=='__main__':
    base = '/mnt/c/Users/jqluo6/Desktop/tcranno_results/tcranno_immuneaccess_analysis/'
    indirs = [base+'CMV_786_1e-4_tcranno_output/CMV_pos',base+'CMV_786_1e-4_tcranno_output/CMV_neg']
    keywords = ['CMV','IE1','pp65']
    all_files_cmv_pos = get_files([indirs[0]], 'tcr2tcr')
    all_files_cmv_neg = get_files([indirs[1]], 'tcr2tcr')
    for i,keyword in enumerate(keywords):
        cmdf, cmp, pmdf, pmp = get_boxplot_inputs(all_files_cmv_pos,all_files_cmv_neg, keyword)
        if i == 0:
            cm = cmdf
            cmps = [cmp]
            pm = pmdf
            pmps = [pmp]
        else:
            cm = pd.concat([cm,cmdf])
            pm = pd.concat([pm,pmdf])
            cmps.append(cmp)
            pmps.append(pmp)
    #plot_boxplot(cm,cmps,keywords,match_type='CM')
    #plot_boxplot(pm,pmps,keywords,match_type='PM')