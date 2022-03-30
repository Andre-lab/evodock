import random
from operator import itemgetter


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
        self.mutate = random.uniform(0.1, 0.7)
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
        self.mutate = random.uniform(0.1, 0.7)
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


class StrategyTriangular(Strategy):
    def update_best_and_worst(self, gen_scores):
        self.best = min(list(enumerate(gen_scores)), key=itemgetter(1))[1]
        self.worst = max(list(enumerate(gen_scores)), key=itemgetter(1))[1]
        print(self.best, self.worst)

    def select_parents(self, j, population, gen_scores):
        candidates = list(range(0, len(population)))
        candidates.remove(j)
        random_index = random.sample(candidates, 3)

        trial_1 = population[random_index[0]]
        trial_2 = population[random_index[1]]
        trial_3 = population[random_index[2]]
        score_1 = gen_scores[random_index[0]]
        score_2 = gen_scores[random_index[1]]
        score_3 = gen_scores[random_index[2]]
        data = [(trial_1, score_1), (trial_2, score_2), (trial_3, score_3)]
        data = sorted(data, key=lambda x: x[1])

        x_1 = data[0][0].genotype
        x_2 = data[1][0].genotype
        x_3 = data[2][0].genotype

        return x_1, x_2, x_3

    def update_weights(self, population):
        p1 = 1.0
        p2 = random.uniform(0.75, 1)
        p3 = random.uniform(0.5, p2)
        self.w1 = p1 / sum([p1, p2, p3])
        self.w2 = p2 / sum([p1, p2, p3])
        self.w3 = p3 / sum([p1, p2, p3])

    def create_donor(self, j, population, gen_scores):
        self.update_weights(population)
        self.mutate = random.uniform(0.1, 0.7)
        x_1, x_2, x_3 = self.select_parents(j, population, gen_scores)
        v_donor = []
        F = self.mutate
        for g in range(len(x_1)):
            conv = self.w1 * x_1[g] + self.w2 * x_2[g] + self.w3 * x_3[g]
            mut = (
                conv
                + F * (x_1[g] - x_2[g])
                + F * (x_2[g] - x_3[g])
                + F * (x_1[g] - x_3[g])
            )
            v_donor.append(mut)

        v_donor = ensure_bounds(v_donor, self.bounds)
        return v_donor
