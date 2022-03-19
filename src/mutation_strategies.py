import random


def ensure_bounds(vec, bounds):
    vec_new = []
    # cycle through each variable in vector
    for i, k in enumerate(vec):

        # variable exceedes the minimum boundary
        if vec[i] < bounds[i][0]:
            vec_new.append(bounds[i][0])

        # variable exceedes the maximum boundary
        if vec[i] > bounds[i][1]:
            vec_new.append(bounds[i][1])

        # the variable is fine
        if bounds[i][0] <= vec[i] <= bounds[i][1]:
            vec_new.append(vec[i])

    return vec_new


class MutationStrategyBuilder:
    def __init__(self, config):
        self.config = config
        self.scheme = config.scheme
        self.mutate = config.mutate
        self.recombination = config.recombination

    def build(self):
        if self.scheme == "BEST":
            return StrategyBest(self.config)
        elif self.scheme == "RANDOM":
            return StrategyRandom(self.config)
        elif self.scheme == "CURRENT":
            return StrategyCurrent(self.config)
        elif self.scheme == "pBEST":
            return StrategyPBest(self.config)
        else:
            return StrategyBest(self.config)


class Strategy:
    def __init__(self, config):
        self.config = config
        self.mutate = config.mutate
        self.recombination = config.recombination
        self.bounds = [(-1, 1)] * 6

    def select_parents(self, j, population, gen_scores):
        candidates = list(range(0, len(population)))
        candidates.remove(j)
        random_index = random.sample(candidates, 3)

        x_1 = population[random_index[0]].genotype
        x_2 = population[random_index[1]].genotype
        x_3 = population[random_index[2]].genotype

        return x_1, x_2, x_3

    def create_donor(self, j, population, gen_scores):
        x_1, x_2, x_3 = self.select_parents(j, population, gen_scores)

        # subtract x3 from x2, and create a new vector (x_diff)
        x_diff = [x_2_i - x_3_i for x_2_i, x_3_i in zip(x_2, x_3)]

        # multiply x_diff by the mutation factor (F) and add to x_1
        v_donor = [
            x_1_i + self.mutate * x_diff_i for x_1_i, x_diff_i in zip(x_1, x_diff)
        ]
        v_donor = ensure_bounds(v_donor, self.bounds)
        return v_donor


class StrategyRandom(Strategy):
    pass


class StrategyBest(Strategy):
    def select_parents(self, j, population, gen_scores):
        candidates = list(range(0, len(population)))
        candidates.remove(j)
        random_index = random.sample(candidates, 3)

        x_1 = population[gen_scores.index(min(gen_scores))].genotype
        x_2 = population[random_index[1]].genotype
        x_3 = population[random_index[2]].genotype

        return x_1, x_2, x_3


class StrategyPBest(Strategy):
    def select_parents(self, j, population, gen_scores):
        candidates = list(range(0, len(population)))
        candidates.remove(j)
        random_index = random.sample(candidates, 3)

        p = random.randint(1, len(population))
        top_individuals = random.sample(list(enumerate(gen_scores)), p)
        if p >= len(population) - 5:  # for computational efficiency
            best_top = sorted(top_individuals, key=lambda x: x[1])[0][0]
        else:
            best_top = gen_scores.index(min(gen_scores))
        x_1 = population[best_top].genotype
        x_2 = population[random_index[1]].genotype
        x_3 = population[random_index[2]].genotype

        return x_1, x_2, x_3


class StrategyCurrent(Strategy):
    def select_parents(self, j, population, gen_scores):
        candidates = list(range(0, len(population)))
        candidates.remove(j)
        random_index = random.sample(candidates, 3)

        x_1 = population[j].genotype
        x_2 = population[random_index[1]].genotype
        x_3 = population[random_index[2]].genotype

        return x_1, x_2, x_3
