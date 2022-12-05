import pandas as pd
from numpy import *
import random as rd
import matplotlib
from matplotlib import pyplot as plt
import os,sys
from scipy.stats import mannwhitneyu
from statannot import add_stat_annotation
import seaborn as sns

def parse_infile(infile, keyword, match_type='CM'):
    complete_match = []
    df = pd.read_csv(infile,sep='\t',low_memory=False,skiprows=4)
    if match_type=='CM':
        m = df[df.columns[1]].tolist()
    else:
        m = df[df.columns[2]].tolist()
    #pm = df[df.columns[2]].tolist()
    #x[:x.rfind(' (')]:float(x[x.rfind(',')+2:x.rfind('%')])/100
    M = {x[:x.rfind(' (')]:float(x[x.rfind('(')+1:x.rfind(',')]) for x in m if x!='-'}
    frac = 0
    for each in list(M.keys()):
        if keyword in each:
            frac+=M[each]
    if frac == 0:
        return None
    return frac
    
def get_files(indirs, file_type): # {'tcr2ept','tcr2ag','tcr2org'}
    all_files=[]
    for indir in indirs:
        all_files += [indir+'/'+x for x in os.listdir(indir) if x.count(file_type)>0]
    return all_files
    
def get_df(dir1,dir2,ftype,keyword, match_type):
    fractions1 = []
    fractions2 = []
    all_files1 = get_files(dir1, ftype)
    all_files2 = get_files(dir2, ftype)
    for infile in all_files1:
        frac = parse_infile(infile, keyword, match_type)
        if frac is None:
            continue
        fractions1.append(frac)
    for infile in all_files2:
        frac = parse_infile(infile, keyword, match_type)
        if frac is None:
            continue
        fractions2.append(frac)
    df = pd.DataFrame(list(zip(fractions1+fractions2,[label1]*len(fractions1)+[label2]*len(fractions2),[keyword]*(len(fractions1)+len(fractions2)))), columns = ['Fraction','Group','Keyword'])
    U1, p = mannwhitneyu(fractions1, fractions2); print(keyword,p)
    return df,p 
    
if __name__=='__main__':
    indirs1 = sys.argv[1].split(',')
    label1 = sys.argv[2]
    indirs2 = sys.argv[3].split(',')
    label2 = sys.argv[4]
    keywords = sys.argv[5].split(',') # ftype1:keyword1,ftype2:keyword2 e.g. tcr2org:SARS-CoV2,tcr2ag:surface_glycoprotein
    match_type = sys.argv[6]
    if match_type.startswith('C') or match_type.startswith('c'):
        mt = 'Complete Matches'
    else:
        mt = 'Predicted Matches'
    ps = []
    for i,each in enumerate(keywords):
        ftype = each.split(':')[0]
        keyword = each.split(':')[1].replace('_',' ')
        if i == 0:
            df,p = get_df(indirs1,indirs2,ftype,keyword, match_type)
            ps.append(p)
        else:
            x, p = get_df(indirs1,indirs2,ftype,keyword, match_type)
            df = pd.concat([df,x])
            ps.append(p)
    plt.figure(figsize=(5,4))#;plt.title(match_type)
    colors = ['forestgreen','gold'];sns.set_palette(array([matplotlib.colors.to_rgba(x) for x in colors]))
    #colors = ['royalblue','orangered'];sns.set_palette(array([matplotlib.colors.to_rgba(x) for x in colors]))
    ax = sns.boxplot(x=df['Keyword'],y=df['Fraction'], hue=df['Group'],width = 0.6, showfliers = False)
    ax.set(xlabel=None);ax.set(ylabel=None);ax.set_title(mt, y=1.0, pad=-24)
    box_pairs = [((k,label1),(k,label2)) for k in df.Keyword.unique()]
    #box_pairs = [(("CMV", "CMV+"), ("CMV", "CMV-")),(("IE1", "CMV+"), ("IE1", "CMV-")),(("pp65", "CMV+"), ("pp65", "CMV-"))]
    add_stat_annotation(ax,data=df,x=df['Keyword'],y=df['Fraction'],hue=df['Group'],box_pairs=box_pairs,perform_stat_test=False,pvalues=ps,text_format='simple',loc='outside', verbose=1)
    ax.legend()
    plt.savefig('Merged_boxplot_'+label1+'_vs_'+label2+'_fraction.'+mt+'.png',dpi=300)
