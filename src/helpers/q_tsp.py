import numpy as np
import cmath
import multiprocessing
from Bio import SeqIO
from tqdm import tqdm
from joblib import Parallel, delayed
import itertools
import sys
from calculate_corr import calculate_corr


def read_from_fasta(filename):
    reads = []
    for record in SeqIO.parse(filename, "fasta"):
        reads.append(str(record.seq))
    return reads


def reads2signal(reads):
    conv_dict = {'A': complex(1, 1),
                 'C': complex(-1, -1),
                 'G': complex(-1, 1),
                 'T': complex(1, -1)}
    return [np.cumsum(np.array(list(map(lambda i: np.angle(conv_dict[i]), read)))) for read in reads]


def reads_pair_corr_coeff(i, j):
    [correlation_matrix[i, j], lags_matrix[i, j]] = calculate_corr(reads_signal[i], reads_signal[j])


if __name__=='__main__':
    fasta_filename = sys.argv[1]
    reads = read_from_fasta(fasta_filename)
    reads_signal = reads2signal(reads)

    N = len(reads_signal)
    correlation_matrix = np.zeros((N, N))
    lags_matrix = np.zeros((N, N))

    num_cores = multiprocessing.cpu_count()
    inputs = tqdm(filter(lambda x: x[0] < x[1], list(itertools.product(range(N), range(N)))))

    Parallel(n_jobs=num_cores, require='sharedmem')(delayed(reads_pair_corr_coeff)(i, j) for i, j in inputs)

    correlation_matrix = np.array(correlation_matrix)
    correlation_matrix = correlation_matrix + correlation_matrix.T
    lags_matrix = np.array(lags_matrix)
    lags_matrix = lags_matrix + lags_matrix.T

    with open(fasta_filename + ".corr", 'wb') as f:
        np.save(f, correlation_matrix)
        
    np.set_printoptions(precision=2)
    print(correlation_matrix)

