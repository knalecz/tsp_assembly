import numpy as np
from sklearn import manifold


def distance_matrix2coordinates(distance_matrix):
    adist = np.array(distance_matrix)
    amax = np.amax(adist)
    adist = np.true_divide(adist, amax)
    mds = manifold.MDS(n_components=2, dissimilarity="precomputed", random_state=6)
    results = mds.fit(adist)
    return results.embedding_

