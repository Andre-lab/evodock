import random


class Individual:
    def __init__(
        self,
        genotype,
        score,
        idx_ligand=1,
        idx_receptor=1,
        rmsd=1000,
        i_sc=1000,
        irms=1000,
    ):
        if idx_ligand == -1:
            idx_ligand = random.randit(1, 99)
        else:
            if idx_ligand > 100:
                print("error at ligand idx")

        if idx_receptor == -1:
            idx_receptor = random.randit(1, 99)
        else:
            if idx_receptor > 100:
                print("error at receptor idx")

        self.idx_receptor = idx_receptor
        self.idx_ligand = idx_ligand
        self.genotype = genotype
        self.score = score
        self.rmsd = rmsd
        self.i_sc = i_sc
        self.irms = irms
