from solver import DWaveLeapHybridSolver, DWaveVRPSolver, GoogleORToolsSolver, Solver
from overlap_detector import CovBasedOverlapDetector, GroundTruthOverlapDetector, OverlapDetector
from metric import Metric
from matrix_scaling_util import MatrixScallingUtil
import json
import os
import sys
from typing import List

import numpy as np
from recordclass import recordclass


path_to_data = 'TODO'
path_to_results = 'TODO'


class Experiment:
    def __init__(self, matrix_filename: str, metric: Metric, coverage: int, error: float, solvers: List[Solver],
                 overlaps_detectors: List[OverlapDetector], sequence_name: str, sequence_length: int) -> None:
        self.matrix_filename = matrix_filename
        self.metric = metric
        self.coverage = coverage
        self.error = error
        self.solvers = solvers
        self.overlaps_detectors = overlaps_detectors
        self.sequence_name = sequence_name
        self.sequence_length = sequence_length
        self.results = dict.fromkeys([s.name() for s in self.solvers])

    ContigSet = recordclass('ContingSet', 'path cost contigs')

    def run(self):
        for solver in self.solvers:
            path, cost = solver.solveTSP()
            print(f'path = {path}')
            print(f'cost = {cost}')
            contig_set = self.ContigSet(path=path,
                                        cost=cost,
                                        contigs=dict.fromkeys([od.name() for od in self.overlaps_detectors]))
            for overlaps_detector in self.overlaps_detectors:
                contig_set.contigs[overlaps_detector.name(
                )] = overlaps_detector.detect_false_overlaps(path[1:-1])
            print(contig_set)
            self.results[solver.name()] = contig_set

    def save(self):
        json_filename = os.path.basename(
            f'exp.{matrix_filename[:matrix_filename.index("fasta")]}{self.metric.name}.json')

        exp_json = {}
        exp_json["matrix_filename"] = self.matrix_filename
        exp_json["sequence_name"] = self.sequence_name
        exp_json["sequence_name"] = self.sequence_name
        exp_json["coverage"] = self.coverage
        exp_json["error"] = self.error
        exp_json["metric"] = self.metric.name
        exp_json["results"] = {}
        for solver, contig_set in self.results.items():
            exp_json["results"][solver] = contig_set.__dict__

        with open(os.path.join(path_to_results, json_filename), 'w') as outfile:
            json.dump(exp_json, outfile)


def int_value_from_filename(filename, value_prefix):
    return int(next(x for x in filename.split('.') if x.startswith(value_prefix)).replace(value_prefix, ''))


def error_value_from_filename(filename):
    filename_parts = filename.split('.')
    for i in range(len(filename_parts)):
        if filename_parts[i].startswith('er'):
            return float(filename_parts[i].replace('er', '')+'.'+filename_parts[i+1])


###############################################################################

if __name__ == "__main__":
    matrix_filename = sys.argv[1]
    metric = Metric.corr if matrix_filename.endswith(
        'corr') else Metric.matches
    distance_matrix = MatrixScallingUtil().read_and_scale_matrix(metric,
                                                                 os.path.join(
                                                                     path_to_data, matrix_filename))
    coverage = int_value_from_filename(matrix_filename, 'cov')
    reads_full_path = os.path.join(
        path_to_data, matrix_filename.replace(f'.{metric.name}', ''))
    overlaps_detectors = [GroundTruthOverlapDetector(distance_matrix,
                                                     reads_full_path,
                                                     int_value_from_filename(matrix_filename, 'L')),
                          CovBasedOverlapDetector(distance_matrix,
                                                  coverage)]
    solvers = [GoogleORToolsSolver(distance_matrix),
               DWaveLeapHybridSolver(distance_matrix, max_time=12, max_attempts=1)]
    solvers = [GoogleORToolsSolver(distance_matrix),
               DWaveVRPSolver(distance_matrix, solver_type="cpu")]
    exp = Experiment(matrix_filename=matrix_filename,
                     metric=metric,
                     coverage=coverage,
                     error=error_value_from_filename(matrix_filename),
                     solvers=solvers,
                     overlaps_detectors=overlaps_detectors,
                     sequence_name=matrix_filename.split('.', 1)[0],
                     sequence_length=int(matrix_filename[matrix_filename.index('.L')+2:].split('.', 1)[0]))
    print(f'\nCalculations for file {matrix_filename}')
    exp.run()
    exp.save()

