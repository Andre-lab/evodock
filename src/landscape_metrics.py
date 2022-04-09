import numpy as np
from numpy.linalg import norm


def fitness_distance_correlation(population, optimal_solution, scfxn):
    # dataframe
    # |ind|fitness|rmsd|genotype|
    optimal_solution = np.array(optimal_solution)
    optimal_solution = np.array(scfxn.convert_positions_to_genotype(optimal_solution))

    F = np.array([ind.score for ind in population])
    D = np.array(
        [norm(np.array(ind.genotype) - optimal_solution, 2) for ind in population]
    )
    var_F = np.std(F)
    var_D = np.std(D)
    mean_F = np.average(F)
    mean_D = np.average(D)
    product_sum = np.array([(f - mean_F) * (d - mean_D) for f, d in zip(F, D)]).sum()

    if product_sum == 0:
        return 1
    else:
        return (product_sum / len(F)) / (var_F * var_D)
