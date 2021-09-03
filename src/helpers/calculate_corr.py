from scipy.stats import pearsonr


def calculate_corr(seq_i, seq_j):
    if len(seq_i)>=len(seq_j):
        longer_signal=seq_i
        shorter_signal=seq_j
    else:
        longer_signal=seq_j
        shorter_signal=seq_i
    
    LD=len(longer_signal)
    LK=len(shorter_signal)
    
    corr=[]
    
    for a in range(2, LK):
        x=shorter_signal[-a:];
        y=longer_signal[:a];
        corr.append(pearsonr(x, y)[0])

    corr[0:9]=[0]*9;
    
    if LD!=LK:
        for b in range(abs(LD-LK)):
            x=shorter_signal;
            y=longer_signal[b:b+LK];
            corr.append(pearsonr(x, y)[0])
    
    for c in range(2, LK):
        x=shorter_signal[0:len(shorter_signal)-c+1]
        y=longer_signal[-LK+c-1:]
        corr.append(pearsonr(x, y)[0])
        
    corr[-9:] = [0]*9
    
    return max(corr), corr.index(max(corr))

