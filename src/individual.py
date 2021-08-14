class Individual:
    def __init__(self, genotype, score, rmsd=1000, i_sc=1000, irms=1000):
        self.genotype = genotype
        self.score = score
        self.rmsd = rmsd
        self.i_sc = i_sc
        self.irms = irms
