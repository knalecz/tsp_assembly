import numpy as np
from Bio import SeqIO

from metric import Metric


class MatrixScallingUtil():

    def __scale_min_max(self, v: np.array) -> np.array:
        vmin = v.min()
        np.fill_diagonal(v, -100)
        vmax = v.max()
        return (v - vmin) / (vmax - vmin)

    def __count_max_read_len(self, reads_full_path: str):
        max_read_len = 0
        for record in SeqIO.parse(reads_full_path, "fasta"):
            if len(record.seq) > max_read_len:
                max_read_len = len(record.seq)
        return max_read_len

    def read_and_scale_matrix(self, metric: Metric, matrix_full_path: str):
        reads_full_path = matrix_full_path.replace(f'.{metric.name}', '')
        read_length = self.__count_max_read_len(reads_full_path)
        distance_matrix = np.load(matrix_full_path)

        if metric == Metric.corr:
            distance_matrix = (self.__scale_min_max(
                1-distance_matrix)*1000).astype(int)
        else:
            distance_matrix = (read_length-distance_matrix).astype(int)

        distance_matrix = np.insert(distance_matrix, 0, values=0, axis=1)
        distance_matrix = np.insert(distance_matrix, 0, values=0, axis=0)
        np.fill_diagonal(distance_matrix, 0)

        return distance_matrix


