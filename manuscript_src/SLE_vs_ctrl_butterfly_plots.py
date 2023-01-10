import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
from numpy import *

def parse_infile(infile):
    percentage = {}
    fraction = {}
    df = pd.read_csv(infile,sep='\t',low_memory=False)
    cm = df[df.columns[1]].tolist()
    CM = {x[:x.rfind(' (')]:float(x[x.rfind(',')+2:x.rfind('%')]) for x in cm if x.count('(')>0}
    CMf = {x[:x.rfind(' (')]:float(x[x.rfind('(')+1:x.rfind(',')]) for x in cm if x.count('(')>0}
    explained = 0
    i = 0
    all_items = list(CM.keys())
    while i < len(all_items):
        item = all_items[i]
        percent = CM[item]
        frac = CMf[item]
        if percent < 1:
            break
        percentage[item[:item.rfind('(')]]=percent
        fraction[item[:item.rfind('(')]]=frac
        i += 1
    return percentage,fraction
	
SLE_p,SLE_f = parse_infile('../suppl_files/SLE.pop.tcr2ag.30.tsv')
C_p,C_f = parse_infile('../suppl_files/CMV_all.pop.tcr2ag.30.tsv')
all_keys = set(SLE_p.keys()).intersection(set(C_p.keys()))
for keys in list(SLE_p.keys()):
    if keys not in all_keys:
        SLE_p.pop(keys)
        SLE_f.pop(keys)
for keys in list(C_p.keys()):
    if keys not in all_keys:
        C_p.pop(keys)
        C_f.pop(keys)
percent = pd.DataFrame(list(zip(list(SLE_p.keys()),list(SLE_p.values()),['SLE']*len(SLE_p)))+
list(zip(list(C_p.keys()),[x*(-1) for x in list(C_p.values())],['Control']*len(C_p))),columns = ['Ag','Percent','Keyword'])
fraction = pd.DataFrame(list(zip(list(SLE_f.keys()),list(SLE_f.values()),['SLE']*len(SLE_f)))+
list(zip(list(C_f.keys()),[x*(-1) for x in list(C_f.values())],['Control']*len(C_f))),columns = ['Ag','Fraction','Keyword'])

plt.figure(figsize=(8,8), dpi= 600)
sns.set_palette(sns.color_palette("Spectral_r", n_colors=len(percent)//2))
group_col = 'Keyword'
order_of_bars = percent.Ag.unique()
alphas = [1,0.5]

plt.subplot(1, 2, 1)
for a, group in zip(alphas, percent[group_col].unique()):
    sns.barplot(x='Percent', y='Ag', data=percent.loc[percent[group_col]==group, :], order=order_of_bars, label=group, 
                saturation=a) 
plt.xlabel("Fraction Percentage (%)",fontsize=8)
plt.ylabel("",fontsize=1)
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)
plt.title("Control vs SLE (Antigens)", fontsize=10)

plt.subplot(1, 2, 2)
for a, group in zip(alphas, fraction[group_col].unique()):
    sns.barplot(x='Fraction', y='Ag', data=fraction.loc[fraction[group_col]==group, :], order=order_of_bars, label=group, 
                saturation=a) 
plt.xlabel("Fraction",fontsize=8)
#plt.ylabel("Organisms",fontsize=6)
plt.xticks(fontsize=8)
plt.yticks([])
plt.ylabel("",fontsize=1)
plt.title("Control vs SLE (Antigens)", fontsize=10)
plt.subplots_adjust(left=0.05,bottom=0.05,right=0.95,top=0.95,wspace=0.05)
plt.savefig('SLE_vs_Control_ag.png', bbox_inches='tight')

SLE_p,SLE_f = parse_infile('../suppl_files/SLE.pop.tcr2org.tsv')
C_p,C_f = parse_infile('../suppl_files/CMV_all.pop.tcr2org.tsv')

all_keys = set(SLE_p.keys()).intersection(set(C_p.keys()))
for keys in list(SLE_p.keys()):
    if keys not in all_keys:
        SLE_p.pop(keys)
        SLE_f.pop(keys)
for keys in list(C_p.keys()):
    if keys not in all_keys:
        C_p.pop(keys)
        C_f.pop(keys)

percent = pd.DataFrame(list(zip(list(SLE_p.keys()),list(SLE_p.values()),['SLE']*len(SLE_p)))+
list(zip(list(C_p.keys()),[x*(-1) for x in list(C_p.values())],['Control']*len(C_p))),columns = ['Org','Percent','Keyword'])

fraction = pd.DataFrame(list(zip(list(SLE_f.keys()),list(SLE_f.values()),['SLE']*len(SLE_f)))+
list(zip(list(C_f.keys()),[x*(-1) for x in list(C_f.values())],['Control']*len(C_f))),columns = ['Org','Fraction','Keyword'])

plt.figure(figsize=(8,4), dpi= 600)

group_col = 'Keyword'
order_of_bars = percent.Org.unique()
#color = ['royalblue','orangered','forestgreen','gold','crimson','violet','peru','greenyellow']
alphas = [1,0.55]
#sns.set_palette(array([matplotlib.colors.to_rgba(x) for x in color]))
sns.set_palette(sns.color_palette("Spectral_r", n_colors=len(percent)//2))

plt.subplot(1, 2, 1)
for a, group in zip(alphas, percent[group_col].unique()):
    sns.barplot(x='Percent', y='Org', data=percent.loc[percent[group_col]==group, :], order=order_of_bars, label=group, 
                saturation=a) 
plt.xlabel("Fraction Percentage (%)",fontsize=8)
plt.ylabel("",fontsize=1)
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)
plt.title("Control vs SLE (Organisms)", fontsize=10)

plt.subplot(1, 2, 2)
for a, group in zip(alphas, fraction[group_col].unique()):
    sns.barplot(x='Fraction', y='Org', data=fraction.loc[fraction[group_col]==group, :], order=order_of_bars, label=group, 
                saturation=a) 
plt.xlabel("Fraction",fontsize=8)
plt.ylabel("",fontsize=1)
plt.xticks(fontsize=8)
plt.yticks([])
plt.title("Control vs SLE (Organisms)", fontsize=10)
plt.subplots_adjust(left=0.1,bottom=0.1,right=0.9,top=0.9,wspace=0.05)
plt.savefig('SLE_vs_Control_org.png', bbox_inches='tight')