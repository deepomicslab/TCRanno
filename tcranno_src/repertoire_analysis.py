import os
from tcranno.core_analysis import *

def load_AO(ref_AO_map=None):
    if ref_AO_map==None:
        d = os.path.dirname(sys.modules['tcranno'].__file__)
        resource3 = 'data/DB_FULL_AG_ORG_key_map.pkl'
        path3 = os.path.join(d,resource3)
    else:
        path3 = ref_AO_map
    with open(path3,'rb') as c:
        AO_map = pickle.load(c)
    return AO_map

def ept_frequency_calculations(dataframe, limit, weight=1.0):
    EPTS = {}
    all_records = [x.split('; ') for x in dataframe['Matched_Epitope(s) [Matched_Antigen->Matched_Organism]'].tolist()]
    all_epts = []
    for rec in all_records:
        epts = []
        for each in rec:
            ep = each.split(' [')[0]
            epts.append(ep)
        all_epts.append(epts)
    frequencies = dataframe['Frequency'].tolist()
    for i,epts in enumerate(all_epts):
        freq = frequencies[i]
        for e in epts:
            if e not in EPTS:
                EPTS[e]=weight*freq
            else:
                EPTS[e]+=weight*freq
    all_epts = list(EPTS.keys())
    for e in all_epts:
        if EPTS[e] < limit:
            EPTS.pop(e)
    return EPTS, sum(frequencies)
    
def ag_frequency_calculations(dataframe, AO_map, limit, weight=1.0):
    EPTS = {}
    AGS = {} #{x[0]:0 for x in AO_map.values()}
    all_records = [x.split('; ') for x in dataframe['Matched_Epitope(s) [Matched_Antigen->Matched_Organism]'].tolist()]
    all_epts = []
    for rec in all_records:
        epts = []
        for each in rec:
            ep = each.split(' [')[0]
            epts.append(ep)
        all_epts.append(epts)   
    frequencies = dataframe['Frequency'].tolist()
    for i,epts in enumerate(all_epts):
        freq = frequencies[i]
        counted_ag = []
        for e in epts:
            if e not in EPTS:
                EPTS[e]=weight*freq
            else:
                EPTS[e]+=weight*freq
            ag = AO_map[e][0]
            org = AO_map[e][1]
            ag = ag + '->' + org 
            if ag not in counted_ag:
                counted_ag.append(ag)
                if ag not in AGS:
                    AGS[ag]=0
                AGS[ag]+=weight*freq
    all_ags = list(AGS.keys())
    for ag in all_ags:
        if AGS[ag] < limit:
            AGS.pop(ag)
    return AGS, sum(frequencies)
    
def org_frequency_calculations(dataframe, AO_map, limit, weight=1.0):
    EPTS = {}
    ORGS = {} #{x[1]:0 for x in AO_map.values()}
    all_records = [x.split('; ') for x in dataframe['Matched_Epitope(s) [Matched_Antigen->Matched_Organism]'].tolist()]
    all_epts = []
    for rec in all_records:
        epts = []
        for each in rec:
            ep = each.split(' [')[0]
            epts.append(ep)
        all_epts.append(epts)
    frequencies = dataframe['Frequency'].tolist()
    for i,epts in enumerate(all_epts):
        freq = frequencies[i]
        counted_org = []
        for e in epts:
            if e not in EPTS:
                EPTS[e]=weight*freq
            else:
                EPTS[e]+=weight*freq
            org = AO_map[e][1]
            if org not in counted_org:
                counted_org.append(org)
                if org not in ORGS:
                    ORGS[org]=0
                ORGS[org]+=weight*freq
    all_orgs = list(ORGS.keys())
    for org in all_orgs:
        if ORGS[org] < limit:
            ORGS.pop(org)
    return ORGS, sum(frequencies)
    
def output_topk_epts(EPTS,sfreq,AO_map,k,weight=1.0):
    counted=[]
    j = 0
    output = []
    while j < k and len(EPTS)>0:
        m = max(EPTS.values())
        for keys in EPTS:
            if EPTS[keys]==m:
                counted.append(keys)
        j += 1
        scounted = sorted(counted)
        all_keys = (',').join(scounted)
        key = scounted[0]
        out = all_keys+' ['+AO_map[key][0]+'->'+AO_map[key][1]+'] ('+str(round(m/weight,5))+', '+str(round(100*m/(sfreq*weight),2))+'%)'
        output.append(out)
        for each in counted:
            EPTS.pop(each)
        counted = []
    return output
    
def output_topk_ags(AGS,sfreq,k,weight=1.0):
    counted=[]
    j = 0
    output = []
    while j < k and len(AGS)>0:
        m = max(AGS.values())
        for keys in AGS:
            if AGS[keys]==m:
                j += 1
                counted.append(keys)
                out = keys+' ('+str(round(m/weight,5))+', '+str(round(100*m/(sfreq*weight),2))+'%)'
                output.append(out)
        for each in counted:
            AGS.pop(each)
        counted = []
    return output
    
def output_topk_orgs(ORGS,sfreq,k,weight=1.0):
    counted=[]
    j = 0
    output = []
    while j < k and len(ORGS)>0:
        m = max(ORGS.values())
        for keys in ORGS:
            if ORGS[keys]==m:
                j += 1
                counted.append(keys)
                out = keys+' ('+str(round(m/weight,5))+', '+str(round(100*m/(sfreq*weight),2))+'%)'
                output.append(out)
        for each in counted:
            ORGS.pop(each)
        counted = []
    return output #('; ').join(output)

def tcr2ept(infile, outprefix, AO_map=None, k=30, hit_rate=0.2, limit=1e-4):
    df = pd.read_csv(infile,sep='\t',low_memory=False,skiprows=3)
    if AO_map is None:
        AO_map = load_AO()
    complete_matches, sfreq_complete = ept_frequency_calculations(df[df['Index'].str.endswith('-0')], limit, 1.0)
    partial_matches, sfreq_partial = ept_frequency_calculations(df[df['Index'].str.endswith('-1')], limit, hit_rate)
    all_epts = set(complete_matches.keys()).union(set(partial_matches.keys()))
    weighted_sum = {e:0 for e in all_epts}
    total_sum = {e:0 for e in all_epts}
    for e in all_epts:
        if e in complete_matches:
            weighted_sum[e]+=complete_matches[e]
            total_sum[e]+=complete_matches[e]
        if e in partial_matches:
            weighted_sum[e]+=partial_matches[e]
            total_sum[e]+=partial_matches[e]/hit_rate    
    outfile = outprefix+'_tcr2ept.tsv'
    pipeline = open(outfile,'w')
    pipeline.write('Features\tComplete_Match (CM)\tPredicted_Match (PM)\tTotal (CM+PM)\tWeighted_Sum (WS)\n')
    cm = output_topk_epts(complete_matches,sfreq_complete,AO_map,k)
    pm = output_topk_epts(partial_matches,sfreq_partial,AO_map,k,hit_rate)
    ts = output_topk_epts(total_sum,sfreq_complete+sfreq_partial,AO_map,k)
    weighted_sum_freq = sfreq_complete+sfreq_partial*hit_rate
    ws = output_topk_epts(weighted_sum, weighted_sum_freq,AO_map,k)
    pipeline.write('Sum_of_Frequency\t'+str(sfreq_complete)+'\t'+str(sfreq_partial)+'\t'+str(sfreq_complete+sfreq_partial)+'\t'+str(weighted_sum_freq)+'\n')
    pipeline.write('Weight\t1.0\t0.2\t-\t-\n')
    pipeline.write(('\t').join(['----------']*5)+'\n') 
    pipeline.write('Rank\tCM:Top_Epitopes (sum_of_freq, %)\tPM:Top_Epitopes (sum_of_freq, %)\tTOTAL:Top_Epitopes (sum_of_freq, %)\tWS:Top_Epitopes (sum_of_freq, %)\n')
    for j in range(min(k,max(len(cm),len(pm),len(ts),len(ws)))):
        line = ['-','-','-','-']
        if len(cm) > j:
            line[0] = cm[j]
        if len(pm) > j:
            line[1] = pm[j]
        if len(ts) > j:
            line[2] = ts[j]
        if len(ws) > j:
            line[3] = ws[j]
        rank = j+1
        pipeline.write(str(rank)+'\t'+('\t').join(line)+'\n')
    pipeline.close()

def tcr2ag(infile,outprefix,AO_map=None, k=20, hit_rate=0.2, limit=1e-4):
    df = pd.read_csv(infile,sep='\t',low_memory=False,skiprows=3)
    if AO_map is None:
        AO_map = load_AO()
    complete_matches, sfreq_complete = ag_frequency_calculations(df[df['Index'].str.endswith('-0')], AO_map, limit, 1.0)
    partial_matches, sfreq_partial = ag_frequency_calculations(df[df['Index'].str.endswith('-1')], AO_map, limit, hit_rate)
    all_ags = set(complete_matches.keys()).union(set(partial_matches.keys()))
    weighted_sum = {e:0 for e in all_ags}
    total_sum = {e:0 for e in all_ags}
    for e in all_ags:
        if e in complete_matches:
            weighted_sum[e]+=complete_matches[e]
            total_sum[e]+=complete_matches[e]
        if e in partial_matches:
            weighted_sum[e]+=partial_matches[e]
            total_sum[e]+=partial_matches[e]/hit_rate
    outfile = outprefix+'_tcr2ag.tsv'
    pipeline = open(outfile,'w')
    pipeline.write('Features\tComplete_Match (CM)\tPredicted_Match (PM)\tTotal (CM+PM)\tWeighted_Sum (WS)\n')
    cm = output_topk_ags(complete_matches,sfreq_complete,k)
    pm = output_topk_ags(partial_matches,sfreq_partial,k,hit_rate)
    ts = output_topk_ags(total_sum,sfreq_complete+sfreq_partial,k)
    weighted_sum_freq = sfreq_complete+sfreq_partial*hit_rate
    ws = output_topk_ags(weighted_sum,weighted_sum_freq,k)
    pipeline.write('Sum_of_Frequency\t'+str(sfreq_complete)+'\t'+str(sfreq_partial)+'\t'+str(sfreq_complete+sfreq_partial)+'\t'+str(weighted_sum_freq)+'\n')
    pipeline.write('Weight\t1.0\t0.2\t-\t-\n')
    pipeline.write(('\t').join(['----------']*5)+'\n') 
    pipeline.write('Rank\tCM:Top_Antigens (sum_of_freq, %)\tPM:Top_Antigens (sum_of_freq, %)\tTOTAL:Top_Antigens (sum_of_freq, pct%)\tWS:Top_Antigens (sum_of_freq, %)\n')
    for j in range(min(k,max(len(cm),len(pm),len(ts),len(ws)))):
        line = ['-','-','-','-']
        if len(cm) > j:
            line[0] = cm[j]
        if len(pm) > j:
            line[1] = pm[j]
        if len(ts) > j:
            line[2] = ts[j]
        if len(ws) > j:
            line[3] = ws[j]
        rank = j+1
        pipeline.write(str(rank)+'\t'+('\t').join(line)+'\n')
    pipeline.close()
    
def tcr2org(infile, outprefix, AO_map=None, k=10, hit_rate=0.2, limit=1e-4):
    df = pd.read_csv(infile,sep='\t',low_memory=False,skiprows=3)
    if AO_map is None:
        AO_map = load_AO()
    complete_matches, sfreq_complete = org_frequency_calculations(df[df['Index'].str.contains('-0')], AO_map, limit, 1.0)
    partial_matches, sfreq_partial = org_frequency_calculations(df[df['Index'].str.endswith('-1')], AO_map, limit, hit_rate)
    all_orgs = set(complete_matches.keys()).union(set(partial_matches.keys()))
    weighted_sum = {e:0 for e in all_orgs}
    total_sum = {e:0 for e in all_orgs}
    for e in all_orgs:
        if e in complete_matches:
            weighted_sum[e]+=complete_matches[e]
            total_sum[e]+=complete_matches[e]
        if e in partial_matches:
            weighted_sum[e]+=partial_matches[e]
            total_sum[e]+=partial_matches[e]/hit_rate
    outfile = outprefix+'_tcr2org.tsv'
    pipeline = open(outfile,'w')
    pipeline.write('Features\tComplete_Match (CM)\tPredicted_Match (PM)\tTotal (CM+PM)\tWeighted_Sum (WS)\n')
    cm = output_topk_orgs(complete_matches,sfreq_complete,k)
    pm = output_topk_orgs(partial_matches,sfreq_partial,k,hit_rate)
    ts = output_topk_orgs(total_sum,sfreq_complete+sfreq_partial,k)
    weighted_sum_freq = sfreq_complete+sfreq_partial*hit_rate
    ws = output_topk_orgs(weighted_sum,weighted_sum_freq,k)
    pipeline.write('Sum_of_Frequency\t'+str(sfreq_complete)+'\t'+str(sfreq_partial)+'\t'+str(sfreq_complete+sfreq_partial)+'\t'+str(weighted_sum_freq)+'\n')
    pipeline.write('Weight\t1.0\t0.2\t-\t-\n')
    pipeline.write(('\t').join(['----------']*5)+'\n')
    pipeline.write('Rank\tCM:Top_Organisms (sum_of_freq, %)\tPM:Top_Organisms (sum_of_freq, %)\tTOTAL:Top_Organisms (sum_of_freq, %)\tWS:Top_Organisms (sum_of_freq, %)\n')
    for j in range(min(k,max(len(cm),len(pm),len(ts),len(ws)))):
        line = ['-','-','-','-']
        if len(cm) > j:
            line[0] = cm[j]
        if len(pm) > j:
            line[1] = pm[j]
        if len(ts) > j:
            line[2] = ts[j]
        if len(ws) > j:
            line[3] = ws[j]
        rank = j+1
        pipeline.write(str(rank)+'\t'+('\t').join(line)+'\n')
    pipeline.close()
