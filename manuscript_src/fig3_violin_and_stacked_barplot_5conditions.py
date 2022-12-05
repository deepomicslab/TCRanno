from numpy import *
from matplotlib import pyplot as plt
import seaborn as sns
import matplotlib, os, palettable
import pandas as pd

def violinplot(box_1, box_2, box_3, box_4, box_5, labels, title):
    my_pal = {"Healthy": 'royalblue', "Vaccine": "crimson", "SARS-CoV2":"forestgreen", "Cancer": "orange", "SLE":"gold"}
    df = pd.DataFrame(list(zip(box_1, box_2, box_3, box_4, box_5)), columns = labels)
    plt.figure(figsize=(4,3))
    plt.title(title)
    ax = sns.violinplot(data=df, order=labels, cut=0, width=0.9,palette=my_pal,saturation=0.95)
    plt.tick_params(axis='both', which='major', labelsize=8)
    plt.savefig(title+'.volin.png',dpi=300)

def parse_tcr2tcr(infile):
    df = pd.read_csv(infile,sep='\t',low_memory=False,skiprows=3)
    cm = df[df['Index'].str.contains('-0')]
    pm = df[df['Index'].str.contains('-1')]
    if len(cm) > 0:
        cm_frac = sum(cm['Frequency'].tolist())
        pm_frac = sum(pm['Frequency'].tolist())
        non = 1 - cm_frac - pm_frac
    else:
        cm_frac = 0
        pm_frac = sum(pm['Frequency'].tolist())
        non = 1 - pm_frac
    return cm_frac,pm_frac,non
    
def get_files(indirs, file_type): # {'tcr2ept','tcr2ag','tcr2org', 'tcr2tcr'}
    all_files=[]
    for indir in indirs:
        all_files += [indir+'/'+x for x in os.listdir(indir) if x.count(file_type)>0]
    return all_files

def frac_df(infiles,group):
    cm_fracs = []
    pm_fracs = []
    nons = []
    groups = []
    nc = 0
    for infile in infiles:
        cm_frac,pm_frac,non = parse_tcr2tcr(infile)
        if cm_frac+pm_frac == 0:
            continue
        if cm_frac == 0:
            nc+=1
        cm_fracs.append(cm_frac/(cm_frac+pm_frac))
        pm_fracs.append(pm_frac/(cm_frac+pm_frac))
        nons.append(non)
    groups = [group]*len(nons)
    df = pd.DataFrame(list(zip(cm_fracs,pm_fracs,groups)),columns = ['CM','PM','Group'])
    #print(nc)
    return df

# draw violin plots
CMV = load('CMV_786_new.npz')
covid1485 = load('covid1485_new.npz')
cancer = load('all_cancers_new.npz')
vaccine = load('covid_vaccine_new.npz')
SLE = load('SLE_new.npz')

colors = ['royalblue','crimson','forestgreen','orange','gold']#,'orangered'] #dodgerblue
sns.set_palette(array([matplotlib.colors.to_rgba(x,alpha=1.0) for x in colors]))

violinplot(CMV['prod'],vaccine['prod'],covid1485['prod'],cancer['prod'],SLE['prod'],['Healthy','Vaccine','SARS-CoV2','Cancer','SLE'],'Productive fraction')
violinplot(CMV['low'],vaccine['low'],covid1485['low'],cancer['low'],SLE['low'],['Healthy','Vaccine','SARS-CoV2','Cancer','SLE'],'Low-frequency clonotype fraction')
violinplot(CMV['med'],vaccine['med'],covid1485['med'],cancer['med'],SLE['med'],['Healthy','Vaccine','SARS-CoV2','Cancer','SLE'],'Medium-frequency clonotype fraction')
violinplot(CMV['high'],vaccine['high'],covid1485['high'],cancer['high'],SLE['high'],['Healthy','Vaccine','SARS-CoV2','Cancer','SLE'],'High-frequency clonotype fraction')
'''
# draw relative fraction stacked barplot
base = '/mnt/c/Users/jqluo6/Desktop/tcranno_results/tcranno_immuneaccess_analysis/'
cmv_dirs = [base+'CMV_786_1e-4_tcranno_output/CMV_pos',base+'CMV_786_1e-4_tcranno_output/CMV_neg',base+'CMV_786_1e-4_tcranno_output/CMV_unknown']
cmvs = get_files(cmv_dirs,'tcr2tcr')
covid1485s = get_files([base+'covid1485_1e-4_tcranno_output'],'tcr2tcr')
cancers = get_files([base+'cancer_1e-4_tcranno_output'],'tcr2tcr')
covidvac_dirs = [base+'covid_vaccine_1e-4_tcranno_output/pre-vaccination',base+'covid_vaccine_1e-4_tcranno_output/post-vaccination']
covidvacs = get_files(covidvac_dirs,'tcr2tcr')
sles = get_files([base+'SLE_1e-4_tcranno_output'],'tcr2tcr')

cmv_df = frac_df(cmvs,'Healthy')
covidvac_df = frac_df(covidvacs,'Vaccine')
cancer_df = frac_df(cancers,'Cancer')
covid1485_df = frac_df(covid1485s,'SARS-CoV2')
sle_df = frac_df(sles,'SLE')
total_df = pd.concat([cmv_df,covidvac_df,cancer_df,covid1485_df,sle_df])

fig, ax = plt.subplots(figsize=(20, 4))
bold_10 = palettable.cartocolors.qualitative.Bold_10.mpl_colors
stack1 = total_df[total_df['Group']=='Healthy']
stack2 = total_df[total_df['Group']=='Vaccine']
stack3 = total_df[total_df['Group']=='SARS-CoV2']
stack4 = total_df[total_df['Group']=='Cancer']
stack5 = total_df[total_df['Group']=='SLE']
stacks = [stack1,stack2,stack3,stack4,stack5]
prev = 0
colors = ['royalblue','crimson','forestgreen','orange','gold']
for i,stack in enumerate(stacks):
    ax.stackplot(arange(prev, prev+len(stack)),stack['PM'].tolist(),stack['CM'].tolist(),colors=[colors[i],bold_10[i*2]],alpha=0.85)
    prev += len(stack)
    #print(len(stack))
ax.set_xticks([])
ax.set_ylabel('Relative Fraction')
ax.set_title('Complete Matches vs Predicted Matches')
fig.savefig('stacked_barplot_cm_vs_pm_5conds.png',dpi=300)

print(mean(total_df['CM'].tolist()),mean(total_df['PM'].tolist()))
'''