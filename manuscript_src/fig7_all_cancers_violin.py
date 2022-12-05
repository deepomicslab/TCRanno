from numpy import *
from matplotlib import pyplot as plt
import seaborn as sns
import matplotlib, os, palettable
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.manifold import TSNE
from sklearn.feature_selection import f_classif

def violinplot(boxes, labels, title): #list(zip(box_1, box_2, box_3, box_4, box_5, box_6, box_7, box_8, box_9))
    my_pal = {"Control": 'purple', "Urothelial": "crimson", "Breast":"orange", "Lung": "gold", "Others":"yellowgreen", "Melanoma":"forestgreen","CMV+":"orangered","CMV-":"darkseagreen","SLE":"royalblue","Healthy":"coral"}
    df = pd.DataFrame(boxes, columns = labels)
    plt.figure(figsize=(7,2)) #plt.figure(figsize=(7,2.5))
    plt.title(title)
    ax = sns.violinplot(data=df, order=labels, cut=0, width=1.0,palette=my_pal,saturation=0.95)
    plt.tick_params(axis='both', which='major', labelsize=8)
    plt.savefig('fig7.'+title+'.9pop.bigger.volin.png',dpi=300)
    
def get_eff_frac(arr):
    med = arr['med']
    high = arr['high']
    ef = [med[i]+high[i] for i in range(len(med))]
    return ef

def get_files(indirs, file_type): # {'tcr2ept','tcr2ag','tcr2org', 'tcr2tcr'}
    all_files=[]
    for indir in indirs:
        all_files += sorted([indir+'/'+x for x in os.listdir(indir) if x.count(file_type)>0])
    return all_files

def parse_infile(infile, keyword, match='total'):
    complete_match = []
    df = pd.read_csv(infile,sep='\t',low_memory=False,skiprows=4)
    cm = df[df.columns[1]].tolist()
    #print(cm)
    pm = df[df.columns[2]].tolist()
    CM = {x[:x.rfind(' (')]:float(x[x.rfind('(')+1:x.rfind(',')]) for x in cm if x!='-'}
    PM = {x[:x.rfind(' (')]:float(x[x.rfind('(')+1:x.rfind(',')]) for x in pm if x!='-'}
    frac = 0
    for each in list(CM.keys()):
        if keyword in each:
            frac+=CM[each]
    if match=='CM':
        return frac
    for each in list(PM.keys()):
        if keyword in each:
            frac+=PM[each]      
    return frac
    
def parse_infile_multi_keywords(infile,keywords, match='total'):
    frac = empty(len(keywords))
    for i,key in enumerate(keywords):
        frac[i] = parse_infile(infile, key, match)
    return frac

def get_frac_multi_keywords(indir,med,high,orgs,ags):
    fraction = []
    all_files = get_files(indir, 'tcr2org')
    #print(len(all_files),all_files[0])
    for i,infile in enumerate(all_files):
        ag_file = infile.replace('tcr2org','tcr2ag')
        frac = parse_infile_multi_keywords(infile, orgs)
        frac2 = parse_infile_multi_keywords(ag_file, ags)
        med_frac = med[i]
        high_frac = high[i]
        #eff_frac = med[i]+high[i]
        fraction.append(list(frac)+list(frac2)+[med_frac,high_frac])
    return array(fraction)
    
def make_feat_df(indirs,label,med,high,orgs,ags,cols):
    df = pd.DataFrame(get_frac_multi_keywords(indirs,med,high,orgs,ags),columns=cols)
    df['Group']=label
    return df
  
# draw violin plots
base = 'tcranno_results/tcranno_immuneaccess_analysis/'
melanoma_dir = base+'melanoma_1e-4_tcranno_output'
lc_dir = base+'lung_cancer_1e-4_tcranno_output'
bc_dir = base+'breast_cancer_1e-4_tcranno_output'
uc_dir = base+'urothelial_cancer_1e-4_tcranno_output'
other_dir = base+'cancer_1e-4_tcranno_output'
cmv_pd = base+'CMV_786_1e-4_tcranno_output/CMV_pos'
cmv_nd = base+'CMV_786_1e-4_tcranno_output/CMV_neg'
ctrl_dir = base+'covid_vaccine_1e-4_tcranno_output/pre-vaccination'
sle_dir = base+'SLE_1e-4_tcranno_output'

melanoma = load(base+'melanoma_sorted.npz')
lc = load(base+'lung_cancer_sorted.npz')
cancer = load(base+'other_cancer_sorted.npz')
uc = load(base+'urothelial_cancer_sorted.npz')
bc = load(base+'breast_cancer_sorted.npz')
ctrl = load(base+'pre_vac_sorted.npz')
cmv_pos = load(base+'CMV_pos_sorted.npz')
cmv_neg = load(base+'CMV_neg_sorted.npz')
sle = load(base+'SLE_sorted.npz')

violinplot(list(zip(cmv_neg['med'],ctrl['med'],cmv_pos['med'],uc['med'],bc['med'],lc['med'],cancer['med'],melanoma['med'],sle['med'])),["CMV-","Control","CMV+","Urothelial","Breast","Lung","Others","Melanoma","SLE"],'Medium-frequency fraction')
violinplot(list(zip(cmv_neg['high'],ctrl['high'],cmv_pos['high'],uc['high'],bc['high'],lc['high'],cancer['high'],melanoma['high'],sle['high'])),["CMV-","Control","CMV+","Urothelial","Breast","Lung","Others","Melanoma","SLE"],'High-frequency fraction')

violinplot(list(zip(cmv_neg['med'],ctrl['med'],cmv_pos['med'],uc['med'],bc['med'],lc['med'],cancer['med'],melanoma['med'],sle['med'])),["CMV-","Control","CMV+","Urothelial","Breast","Lung","Others","Melanoma","SLE"],'Medium-frequency fraction')
violinplot(list(zip(cmv_neg['high'],ctrl['high'],cmv_pos['high'],uc['high'],bc['high'],lc['high'],cancer['high'],melanoma['high'],sle['high'])),["CMV-","Control","CMV+","Urothelial","Breast","Lung","Others","Melanoma","SLE"],'High-frequency fraction')
