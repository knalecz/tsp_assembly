import sys
import os
import json
from enum import Enum


Clazz = Enum('Clazz', 'overlap contig')


if __name__ == '__main__':
    with open(sys.argv[1], 'r') as f:
        exp_result = json.load(f)

    print("\n" + os.path.basename(sys.argv[1]))
    for solver in ['GoogleORToolsSolver', 'DWaveVRPSolver']:
        result = exp_result['results'][solver]
        path = result['path']
        cost = result['cost']
        real_contigs = [tuple(c) for c in result['contigs']
                        ['GroundTruthOverlapDetector']]
        calculated_contigs = [
            tuple(c) for c in result['contigs']['CovBasedOverlapDetector']]
        TP, TN, FP, FN = (0, 0, 0, 0)
        for pair in zip(path[:-1], path[1:]):
            actual_class = Clazz.contig if pair in real_contigs else Clazz.overlap
            calculated_class = Clazz.contig if pair in calculated_contigs else Clazz.overlap

            if actual_class == calculated_class == Clazz.overlap:
                TP += 1
            elif actual_class == calculated_class == Clazz.contig:
                TN += 1
            elif (actual_class != calculated_class) and (calculated_class == Clazz.overlap):
                FP += 1
            elif (actual_class != calculated_class) and (calculated_class == Clazz.contig):
                FN += 1
        accuracy = 1.0 if (
            TP + TN + FP + FN) == 0 else (TP + TN)/(TP + TN + FP + FN)

        print(f'{solver} -> cost: {cost}, real: {len(real_contigs)}, calc: {len(calculated_contigs)}, TP: {TP}, FP: {FP}, accuracy: {accuracy:.2f}')


