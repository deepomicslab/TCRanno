from numpy import *
import Levenshtein

def FindVDJ(cdr3,VJ_map):
    if cdr3.startswith('CC'):
        cdr3 = cdr3[1:]
    best_len = 0
    best_v = ''
    best_j = ''
    best_d = ''    
    for v in VJ_map:
        if not cdr3.startswith(v):
            continue
        for j in VJ_map[v]:
            if len(v)+len(j) <= best_len:
                break
            if not cdr3[len(v):].endswith(j):
                continue
            if len(v)+len(j) > best_len:
                best_len = len(v)+len(j)
                best_v = v
                best_j = j
                best_d = cdr3[len(v):len(cdr3)-len(j)]
                break            
    if best_len==0:
        if not cdr3.startswith('C'):
            cdr3 = 'C'+cdr3
        if not cdr3.endswith('F'):
            cdr3 = cdr3+'F'
        return FindVDJ(cdr3,VJ_map)
    return best_v, best_d, best_j

def D_segment_levenshtein_score(query, candidates, candidates_vdj, VJ_map):
    scores = []
    qv,qd,qj = FindVDJ(query,VJ_map)
    for target in candidates:
        tv,td,tj = candidates_vdj[target]
        scores.append(Levenshtein.distance(qd,td))
    maxs = max(scores)
    norm_scores = [1-s/maxs for s in scores]
    return array(norm_scores)

def levenshtein_score(query,candidates):
    scores = []
    if len(query) <= 10:
        trim = 0
    else:
        trim = math.ceil((len(query)-10)/4)
    if trim == 0:
        query_seq = query
    else:
        query_seq = query[trim:-trim]
    for target in candidates:
        if trim==0:
            target_seq = target
        else:
            target_seq = target[trim:-trim]
        scores.append(1-Levenshtein.distance(query_seq,target_seq)/len(query))
    return array(scores)
    #maxs = max(scores)#len(query_seq)
    #norm_scores = [1-s/maxs for s in scores]
    #return array(norm_scores)

def pooling(query_seq, candidates, n=100):
    norm_score = levenshtein_score(query_seq, candidates)
    sorted_index_array = argsort(norm_score)[::-1]
    if norm_score[sorted_index_array[0]]==norm_score[sorted_index_array[n-1]]:
        j=n
        while norm_score[sorted_index_array[j]]==norm_score[sorted_index_array[0]]:
            j+=1
        pool_idxs = sorted_index_array[:j]
    else:
        pool_idxs = sorted_index_array[:n]
    #selected_cdr3s = array(candidates)[pool_idxs]
    #norm_score = levenshtein_score(query_seq, selected_cdr3s,1)
    return pool_idxs, norm_score[pool_idxs]
    
def FindClosestSeq_batch(query_cdr3s, index_start, fractions, db_cdr3s, db_latent, query_latent, DB, DB_VDJ, AO_map, VJ_map, k):
    output = ''
    index = index_start
    for i in range(len(query_cdr3s)):
        query = query_cdr3s[i]
        if fractions is None:
            output += str(index)+'\tInput\t'+query+'\t1\n'
        else:
            output += str(index)+'\tInput\t'+query+'\t'+str(fractions[i])+'\n'    
        if query in db_cdr3s:
            epts = DB[query]
            ept_ag_orgs = []
            for e in epts: 
                ept_ag_orgs.append(e+' ['+AO_map[e][0]+'->'+AO_map[e][1]+']')
            if fractions is None:
                output += str(index)+'-0\tComplete_Match\t'+query+'\t1\t'+('; ').join(ept_ag_orgs)+'\t1.0\t1.0|1.0|1.0\t0.333333|0.333333|0.333333\n'
            else:
                output += str(index)+'-0\tComplete_Match\t'+query+'\t'+str(fractions[i])+'\t'+('; ').join(ept_ag_orgs)+'\t1.0\t1.0|1.0|1.0\t0.333333|0.333333|0.333333\n'
            index += 1
            continue
        pool_idxs, norm_score1 = pooling(query, db_cdr3s)
        pool = list(array(db_cdr3s)[pool_idxs])
        pool_latent = db_latent[pool_idxs]
        norm_score2 = D_segment_levenshtein_score(query, pool, DB_VDJ, VJ_map)
        dotproduct = dot(pool_latent,query_latent[i])
        var1 = var(norm_score1)
        var2 = var(norm_score2) 
        var3 = var(dotproduct) 
        s=var1+var2+var3
        if var1==0:
            t1 = 0
        else:
            t1 = s/var1
        t2 = s/var2
        t3 = s/var3
        st = (t1+t2+t3)
        w1 = t1/st
        w2 = t2/st
        w3 = t3/st        
        final_scores = array([w1*norm_score1[j]+w2*norm_score2[j]+w3*dotproduct[j] for j in range(len(dotproduct))])
        sorted_index_array = argsort(final_scores)[::-1]
        topk_idx = sorted_index_array[:k]
        topk_final_scores = final_scores[sorted_index_array[:k]]
        topk_norm_score1 = norm_score1[sorted_index_array[:k]]
        topk_norm_score2 = norm_score2[sorted_index_array[:k]]
        topk_dotproduct = dotproduct[sorted_index_array[:k]]
        for j,idxs in enumerate(topk_idx):
            rank = j+1
            weighted_score = topk_final_scores[j]
            c1 = topk_norm_score1[j]
            c2 = topk_norm_score2[j]
            c3 = topk_dotproduct[j]
            candidate_cdr3 = pool[idxs]
            epts = DB[candidate_cdr3]
            ept_ag_orgs = []
            for e in epts:
                ept_ag_orgs.append(e+' ['+AO_map[e][0]+'->'+AO_map[e][1]+']')
            if fractions is None:
                output += str(index)+'-'+str(rank)+'\t'+str(rank)+'\t'+candidate_cdr3+'\t1\t'+('; ').join(ept_ag_orgs)+'\t'+str(weighted_score)+'\t'+('|').join([str(c1),str(c2),str(c3)])+'\t'+('|').join([str(w1),str(w2),str(w3)])+'\n'
            else:
                output += str(index)+'-'+str(rank)+'\t'+str(rank)+'\t'+candidate_cdr3+'\t'+str(fractions[i])+'\t'+('; ').join(ept_ag_orgs)+'\t'+str(weighted_score)+'\t'+('|').join([str(c1),str(c2),str(c3)])+'\t'+('|').join([str(w1),str(w2),str(w3)])+'\n'  
        index += 1
    return output
    
def FindClosestSeq_batch_lst(query_cdr3s, index_start, fractions, db_cdr3s, db_latent, query_latent, DB, DB_VDJ, AO_map, VJ_map, k):
    output = []
    index = index_start
    for i in range(len(query_cdr3s)):
        query = query_cdr3s[i]
        if fractions is None:
            output.append([str(index),'Input',query,1])
        else:
            output.append([str(index),'Input',query,fractions[i]])    
        if query in db_cdr3s:
            epts = DB[query]
            ept_ag_orgs = []
            for e in epts: 
                ept_ag_orgs.append(e+' ['+AO_map[e][0]+'->'+AO_map[e][1]+']')
            if fractions is None:
                output.append([str(index)+'-0','Complete_Match',query,1,('; ').join(ept_ag_orgs),1.0,'1.0|1.0|1.0','0.333333|0.333333|0.333333'])
            else:
                output.append([str(index)+'-0','Complete_Match',query,fractions[i],('; ').join(ept_ag_orgs),1.0,'1.0|1.0|1.0','0.333333|0.333333|0.333333'])
            index += 1
            continue
        pool_idxs, norm_score1 = pooling(query, db_cdr3s)
        pool = list(array(db_cdr3s)[pool_idxs])
        pool_latent = db_latent[pool_idxs]
        norm_score2 = D_segment_levenshtein_score(query, pool, DB_VDJ, VJ_map)
        dotproduct = dot(pool_latent,query_latent[i])
        var1 = var(norm_score1)
        var2 = var(norm_score2) 
        var3 = var(dotproduct) 
        s=var1+var2+var3
        if var1==0:
            t1 = 0
        else:
            t1 = s/var1
        t2 = s/var2
        t3 = s/var3
        st = (t1+t2+t3)
        w1 = t1/st
        w2 = t2/st
        w3 = t3/st        
        final_scores = array([w1*norm_score1[j]+w2*norm_score2[j]+w3*dotproduct[j] for j in range(len(dotproduct))])
        sorted_index_array = argsort(final_scores)[::-1]
        topk_idx = sorted_index_array[:k]
        topk_final_scores = final_scores[sorted_index_array[:k]]
        topk_norm_score1 = norm_score1[sorted_index_array[:k]]
        topk_norm_score2 = norm_score2[sorted_index_array[:k]]
        topk_dotproduct = dotproduct[sorted_index_array[:k]]
        for j,idxs in enumerate(topk_idx):
            rank = j+1
            weighted_score = topk_final_scores[j]
            c1 = topk_norm_score1[j]
            c2 = topk_norm_score2[j]
            c3 = topk_dotproduct[j]
            candidate_cdr3 = pool[idxs]
            epts = DB[candidate_cdr3]
            ept_ag_orgs = []
            for e in epts:
                ept_ag_orgs.append(e+' ['+AO_map[e][0]+'->'+AO_map[e][1]+']')
            if fractions is None:
                output.append([str(index)+'-'+str(rank),str(rank),candidate_cdr3,1,('; ').join(ept_ag_orgs),weighted_score,('|').join([str(c1),str(c2),str(c3)]),('|').join([str(w1),str(w2),str(w3)])])
            else:
                output.append([str(index)+'-'+str(rank),str(rank),candidate_cdr3,fractions[i],('; ').join(ept_ag_orgs),weighted_score,('|').join([str(c1),str(c2),str(c3)]),('|').join([str(w1),str(w2),str(w3)])])   
        index += 1
    return output

def FindClosestSeq(query, index, fractions, db_cdr3s, db_latent, query_latent, DB, DB_VDJ, AO_map, VJ_map, k):
    output = ''
    if fractions is None:
        output += str(index)+'\tInput\t'+query+'\t1\n'
    else:
        output += str(index)+'\tInput\t'+query+'\t'+str(fractions)+'\n'  
    if query in db_cdr3s:
        epts = DB[query]
        ept_ag_orgs = []
        for e in epts: 
            ept_ag_orgs.append(e+' ['+AO_map[e][0]+'->'+AO_map[e][1]+']')
        if fractions is None:
            output += str(index)+'-0\tComplete_Match\t'+query+'\t1\t'+('; ').join(ept_ag_orgs)+'\t1.0\t1.0|1.0|1.0\t0.333333|0.333333|0.333333\n'
        else:
            output += str(index)+'-0\tComplete_Match\t'+query+'\t'+str(fractions)+'\t'+('; ').join(ept_ag_orgs)+'\t1.0\t1.0|1.0|1.0\t0.333333|0.333333|0.333333\n'    
        return output
    pool_idxs, norm_score1 = pooling(query, db_cdr3s)
    pool = list(array(db_cdr3s)[pool_idxs])
    pool_latent = db_latent[pool_idxs]
    norm_score2 = D_segment_levenshtein_score(query, pool, DB_VDJ, VJ_map)
    dotproduct = dot(pool_latent,query_latent)
    var1 = var(norm_score1)
    var2 = var(norm_score2) 
    var3 = var(dotproduct) 
    s=var1+var2+var3
    if var1==0:
        t1 = 0
    else:
        t1 = s/var1
    t2 = s/var2
    t3 = s/var3
    st = (t1+t2+t3)
    w1 = t1/st
    w2 = t2/st
    w3 = t3/st        
    final_scores = array([w1*norm_score1[j]+w2*norm_score2[j]+w3*dotproduct[j] for j in range(len(dotproduct))])
    sorted_index_array = argsort(final_scores)[::-1]
    topk_idx = sorted_index_array[:k]
    topk_final_scores = final_scores[sorted_index_array[:k]]
    topk_norm_score1 = norm_score1[sorted_index_array[:k]]
    topk_norm_score2 = norm_score2[sorted_index_array[:k]]
    topk_dotproduct = dotproduct[sorted_index_array[:k]]
    for j,idxs in enumerate(topk_idx):
        rank = j+1
        weighted_score = topk_final_scores[j]
        c1 = topk_norm_score1[j]
        c2 = topk_norm_score2[j]
        c3 = topk_dotproduct[j]
        candidate_cdr3 = pool[idxs]
        epts = DB[candidate_cdr3]
        ept_ag_orgs = []
        for e in epts:
            ept_ag_orgs.append(e+' ['+AO_map[e][0]+'->'+AO_map[e][1]+']')
        if fractions is None:
            output += str(index)+'-'+str(rank)+'\t'+str(rank)+'\t'+candidate_cdr3+'\t1\t'+('; ').join(ept_ag_orgs)+'\t'+str(weighted_score)+'\t'+('|').join([str(c1),str(c2),str(c3)])+'\t'+('|').join([str(w1),str(w2),str(w3)])+'\n'
        else:
            output += str(index)+'-'+str(rank)+'\t'+str(rank)+'\t'+candidate_cdr3+'\t'+str(fractions)+'\t'+('; ').join(ept_ag_orgs)+'\t'+str(weighted_score)+'\t'+('|').join([str(c1),str(c2),str(c3)])+'\t'+('|').join([str(w1),str(w2),str(w3)])+'\n'            
    return output
    
def FindClosestSeq_lst(query, index, fractions, db_cdr3s, db_latent, query_latent, DB, DB_VDJ, AO_map, VJ_map, k):
    output = []
    if fractions is None:
        output.append([str(index),'Input',query,1])
    else:
        output.append([str(index),'Input',query,fractions])
    if query in db_cdr3s:
        epts = DB[query]
        ept_ag_orgs = []
        for e in epts: 
            ept_ag_orgs.append(e+' ['+AO_map[e][0]+'->'+AO_map[e][1]+']')
        if fractions is None:
            output.append([str(index)+'-0','Complete_Match',query,1,('; ').join(ept_ag_orgs),1.0,'1.0|1.0|1.0','0.333333|0.333333|0.333333'])
        else:
            output.append([str(index)+'-0','Complete_Match',query,fractions,('; ').join(ept_ag_orgs),1.0,'1.0|1.0|1.0','0.333333|0.333333|0.333333'])   
        return output
    pool_idxs, norm_score1 = pooling(query, db_cdr3s)
    pool = list(array(db_cdr3s)[pool_idxs])
    pool_latent = db_latent[pool_idxs]
    norm_score2 = D_segment_levenshtein_score(query, pool, DB_VDJ, VJ_map)
    dotproduct = dot(pool_latent,query_latent)
    var1 = var(norm_score1)
    var2 = var(norm_score2) 
    var3 = var(dotproduct) 
    s=var1+var2+var3
    if var1==0:
        t1 = 0
    else:
        t1 = s/var1
    t2 = s/var2
    t3 = s/var3
    st = (t1+t2+t3)
    w1 = t1/st
    w2 = t2/st
    w3 = t3/st        
    final_scores = array([w1*norm_score1[j]+w2*norm_score2[j]+w3*dotproduct[j] for j in range(len(dotproduct))])
    sorted_index_array = argsort(final_scores)[::-1]
    topk_idx = sorted_index_array[:k]
    topk_final_scores = final_scores[sorted_index_array[:k]]
    topk_norm_score1 = norm_score1[sorted_index_array[:k]]
    topk_norm_score2 = norm_score2[sorted_index_array[:k]]
    topk_dotproduct = dotproduct[sorted_index_array[:k]]
    for j,idxs in enumerate(topk_idx):
        rank = j+1
        weighted_score = topk_final_scores[j]
        c1 = topk_norm_score1[j]
        c2 = topk_norm_score2[j]
        c3 = topk_dotproduct[j]
        candidate_cdr3 = pool[idxs]
        epts = DB[candidate_cdr3]
        ept_ag_orgs = []
        for e in epts:
            ept_ag_orgs.append(e+' ['+AO_map[e][0]+'->'+AO_map[e][1]+']')
        if fractions is None:
            output.append([str(index)+'-'+str(rank),str(rank),candidate_cdr3,1,('; ').join(ept_ag_orgs),weighted_score,('|').join([str(c1),str(c2),str(c3)]),('|').join([str(w1),str(w2),str(w3)])])
        else:
            output.append([str(index)+'-'+str(rank),str(rank),candidate_cdr3,fractions,('; ').join(ept_ag_orgs),weighted_score,('|').join([str(c1),str(c2),str(c3)]),('|').join([str(w1),str(w2),str(w3)])])   
    return output
