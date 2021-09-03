def calculate_common_nucleotides(x, y):
    assert(len(x) == len(y))
    return sum([1 if x[i]==y[i] else 0 for i in range(len(x))])


def calculate_matches(seq_i, seq_j):
    if len(seq_i)>=len(seq_j):
        sig_delsi=seq_i
        sig_kratsi=seq_j
    else:
        sig_delsi=seq_j
        sig_kratsi=seq_i
    
    LD=len(sig_delsi)
    LK=len(sig_kratsi)
    
    score=[]
    
    for a in range(2, LK):
        x=sig_kratsi[-a:];
        y=sig_delsi[:a];
        score.append(calculate_common_nucleotides(x, y))

    score[0:9]=[0]*9;
    
    if LD!=LK:
        for b in range(abs(LD-LK)):
            x=sig_kratsi;
            y=sig_delsi[b:b+LK];
            score.append(calculate_common_nucleotides(x, y))
    
    for c in range(2, LK):
        x=sig_kratsi[0:len(sig_kratsi)-c+1]
        y=sig_delsi[-LK+c-1:]
        score.append(calculate_common_nucleotides(x, y))
        
    score[-9:] = [0]*9
    
    return max(score), score.index(max(score))
