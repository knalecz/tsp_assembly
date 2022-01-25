import abc
from typing import List

import numpy as np
from Bio import SeqIO
from recordclass import recordclass
from scipy import stats


class OverlapDetector(metaclass=abc.ABCMeta):
    def __init__(self, distance_matrix: np.array) -> None:
        self.distance_matrix = distance_matrix

    @abc.abstractmethod
    def detect_false_overlaps(self, path) -> List[int]:
        pass

    def name(self):
        return type(self).__name__


class CovBasedOverlapDetector(OverlapDetector):
    def __init__(self, distance_matrix: np.array, coverage: int) -> None:
        super().__init__(distance_matrix)
        self.coverage = coverage

    def __get_cov_based_threshold(self, node):
        # +2 for two zeros in each row of the distance matrix (artificial zero and zero on the diagonal)
        k = self.coverage + 2
        return np.partition(self.distance_matrix[node], k)[k]

    def detect_false_overlaps(self, path) -> List[int]:
        false_overlaps = []
        for node, next_node in zip(path[:-1], path[1:]):
            threshold = self.__get_cov_based_threshold(node)
            if self.distance_matrix[node][next_node] > threshold:
                false_overlaps.append((node, next_node))
        return false_overlaps


class GroundTruthOverlapDetector(OverlapDetector):
    def __init__(self, distance_matrix: np.array, reads_filename: str, initial_seq_len: int) -> None:
        super().__init__(None)
        self.reads_filename = reads_filename
        self.initial_seq_len = initial_seq_len

    Read = recordclass('Read', 'id seq start_pos end_pos')

    def __get_reads(self):
        reads = {}
        for record in SeqIO.parse(self.reads_filename, "fasta"):
            idx = int(record.id.split('_')[-2])
            pos = record.id.split('_')[-1].split(':')[1].split("-")
            read_start, read_end = int(pos[0]), int(pos[1])
            reads[idx] = self.Read(id=idx, seq=str(record.seq),
                                   start_pos=read_start, end_pos=read_end)
        return reads

    def __is_overlapping_interval(self, a_start, a_end, b_start, b_end):
        return not ((b_start >= a_end) or (a_start >= b_end))

    def __is_overlap(self, read_a, read_b):

        if (read_a.end_pos > self.initial_seq_len) and (read_b.end_pos > self.initial_seq_len):
            return True

        if (read_a.end_pos > self.initial_seq_len):
            return self.__is_overlapping_interval(read_a.start_pos, self.initial_seq_len, read_b.start_pos, read_b.end_pos) or \
                self.__is_overlapping_interval(
                    0, read_a.end_pos % self.initial_seq_len, read_b.start_pos, read_b.end_pos)

        if (read_b.end_pos > self.initial_seq_len):
            return self.__is_overlapping_interval(read_b.start_pos, self.initial_seq_len, read_a.start_pos, read_a.end_pos) or \
                self.__is_overlapping_interval(
                    0, read_b.end_pos % self.initial_seq_len, read_a.start_pos, read_a.end_pos)

        return self.__is_overlapping_interval(read_a.start_pos, read_a.end_pos, read_b.start_pos, read_b.end_pos)

    def detect_false_overlaps(self, path) -> List[int]:
        reads = self.__get_reads()

        false_overlaps = []
        for node, next_node in zip(path[:-1], path[1:]):
            if not self.__is_overlap(reads[node], reads[next_node]):
                false_overlaps.append((node, next_node))
        return false_overlaps
