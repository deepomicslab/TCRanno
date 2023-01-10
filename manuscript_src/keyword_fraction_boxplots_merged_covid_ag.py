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
    cm = df[df.columns[1]].tolist()
    pm = df[df.columns[2]].tolist()
    #pm = df[df.columns[2]].tolist()
    #x[:x.rfind(' (')]:float(x[x.rfind(',')+2:x.rfind('%')])/100
    CM = {x[:x.rfind(' (')]:float(x[x.rfind('(')+1:x.rfind(',')]) for x in cm if x!='-'}
    PM = {x[:x.rfind(' (')]:float(x[x.rfind('(')+1:x.rfind(',')]) for x in pm if x!='-'}
    frac = 0
    for each in list(CM.keys()):
        if keyword in each:
            frac+=CM[each]
    if match_type=='CM':
        return frac
    elif match_type=='PM':
        frac = 0
    for each in list(PM.keys()):
        if keyword in each:
            frac+=PM[each]
    return frac
    
def get_files(indirs, file_type): # {'tcr2ept','tcr2ag','tcr2org'}
    all_files=[]
    for indir in indirs:
        all_files += [indir+'/'+x for x in os.listdir(indir) if x.count(file_type)>0]
    return all_files
    
def get_frac(indir,ftype,match_type,keyword):
    fraction = []
    all_files = get_files(indir, ftype)
    for infile in all_files:
        frac = parse_infile(infile, keyword, match_type)
        if frac is None:
            continue    
        fraction.append(frac)
    return fraction

def get_df(dir1,dir2,dir3,dir4,dir5,ftype,keyword,match_type,labels):
    fractions1 = get_frac([dir1],ftype,match_type,keyword)
    fractions2 = get_frac([dir2],ftype,match_type,keyword)
    fractions3 = get_frac([dir3],ftype,match_type,keyword)
    fractions4 = get_frac([dir4],ftype,match_type,keyword)
    fractions5 = get_frac([dir5],ftype,match_type,keyword)
    df = pd.DataFrame(list(zip(fractions1+fractions2+fractions3+fractions4+fractions5,[labels[0]]*len(fractions1)+[labels[1]]*len(fractions2)+[labels[2]]*len(fractions3)+[labels[3]]*len(fractions4)+[labels[4]]*len(fractions5),[keyword]*(len(fractions1)+len(fractions2)+len(fractions3)+len(fractions4)+len(fractions5)))), columns = ['Fraction','Group','Keyword'])
    U1, p1 = mannwhitneyu(fractions1, fractions2) # pre-vs-post-vaccination
    U2, p2 = mannwhitneyu(fractions1, fractions5) # pre-vs-infected
    U3, p3 = mannwhitneyu(fractions3, fractions4) # conval vs deceased
    return df, [p1,p2,p3]
    
def plotting(match_type,labels,keywords,indirs):
    if match_type.startswith('C') or match_type.startswith('c'):
        mt = 'Complete Matches'
    elif match_type.startswith('P') or match_type.startswith('p'):
        mt = 'Predicted Matches'
    else:
        mt = 'Total (CM+PM)'
    #print(indirs[0],indirs[1],indirs[2],indirs[3],indirs[4])
    ps = []
    for i,each in enumerate(keywords):
        ftype = each.split(':')[0]
        keyword = each.split(':')[1].replace('_',' ')
        if i == 0:
            df,p = get_df(indirs[0],indirs[1],indirs[2],indirs[3],indirs[4],ftype,keyword,match_type,labels)
            for e in p:
                ps.append(e)
        else:
            x, p = get_df(indirs[0],indirs[1],indirs[2],indirs[3],indirs[4],ftype,keyword,match_type,labels)
            df = pd.concat([df,x])
            for e in p:
                ps.append(e)
    plt.figure(figsize=(18,9)) #(5,4)
    colors = ['royalblue','orangered','forestgreen','yellowgreen','gold'];sns.set_palette(array([matplotlib.colors.to_rgba(x) for x in colors]))#sns.set_palette(array(sns.color_palette("tab10"))[idx:idx+2])
    ax = sns.boxplot(x=df['Keyword'],y=df['Fraction'], hue=df['Group'],width = 0.6, showfliers = False)
    ax.set(xlabel=None);ax.set_title(mt, y=1.0, pad=-24, fontsize = 20)
    plt.ylabel('Fraction', fontsize=12);ax.tick_params(axis='both', which='major', labelsize=10)
    #plt.title(mt,fontsize = 20,)
    box_pair1 = [((k,labels[0]),(k,labels[1])) for k in df.Keyword.unique()]
    box_pair2 = [((k,labels[0]),(k,labels[4])) for k in df.Keyword.unique()]
    box_pair3 = [((k,labels[2]),(k,labels[3])) for k in df.Keyword.unique()]
    box_pairs = []
    for i in range(len(box_pair1)):
        box_pairs.append(box_pair1[i])
        box_pairs.append(box_pair2[i])
        box_pairs.append(box_pair3[i])
    #box_pairs = [((x.split(':')[1].replace('_',' '), "pre-vaccination"), (x.split(':')[1].replace('_',' '), "post-vaccination")),(x.split(':')[1].replace('_',' '), "pre-vaccination"), (x.split(':')[1].replace('_',' '), "infected")),((x.split(':')[1].replace('_',' '), "convalescent"), (x.split(':')[0], "deceased")) for x in keywords]
    add_stat_annotation(ax,data=df,x=df['Keyword'],y=df['Fraction'],hue=df['Group'],box_pairs=box_pairs,perform_stat_test=False,pvalues=ps,text_format='simple',loc='outside', verbose=1)
    plt.legend(prop={'size': 15}) #
    plt.savefig('fig5.Merged_boxplot_COVID_ag_fraction.sat95.'+mt+'.clean.png',dpi=300, bbox_inches='tight')
    
if __name__=='__main__':
    base = '/mnt/c/Users/jqluo6/Desktop/tcranno_results/tcranno_immuneaccess_analysis/'
    indirs = [base+'covid_vaccine_1e-4_tcranno_output/pre-vaccination',base+'covid_vaccine_1e-4_tcranno_output/post-vaccination',base+'covid_conval_tcranno_output',base+'covid_deceased_tcranno_output',base+'covid1485_1e-4_tcranno_output']
    labels = ['pre-vaccination','post-vaccination','convalescent','deceased','exposed/infected']
    keywords = ['tcr2ag:ORF1ab','tcr2ag:surface_glycoprotein','tcr2ag:ORF7b','tcr2ag:ORF3a','tcr2ag:nucleocapsid_phosphoprotein','tcr2ag:membrane_glycoprotein','tcr2ag:ORF10']
    plotting('CM',labels,keywords,indirs)
    plotting('PM',labels,keywords,indirs)
    plotting('Total',labels,keywords,indirs)

