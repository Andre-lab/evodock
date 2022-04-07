import numpy as np
from numpy.linalg import norm


def fitness_distance_correlation(population, optimal_solution):
    # dataframe
    # |ind|fitness|rmsd|genotype|
    optimal_solution = np.array(optimal_solution)
    F = np.array([ind.score for ind in population])
    D = np.array([norm(ind.genotype - optimal_solution, 2) for ind in population])
    var_F = np.var(F)
    var_D = np.var(D)
    covariance = np.array(
        [(f - np.mean(F)) * (d - np.mean(D)) for f, d in zip(F, D)]
    ).sum()
    FDC = covariance / (var_F * var_D)
    return FDC
