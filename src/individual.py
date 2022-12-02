import random


class Individual:
    def __init__(
        self,
        idx,
        genotype,
        score,
        idx_ligand=1,
        idx_receptor=1,
        idx_subunit=1,
        rmsd=1000,
        i_sc=1000,
        irms=1000,
        ligand_name="",
        receptor_name="",
        subunit_name="",
    ):

        # FIXME: as you can have larger ensembles than 100 this does not make sense
        # if idx_ligand == -1:
        #     idx_ligand = random.randit(1, 99)
        # else:
        #     if idx_ligand > 100:
        #         print("error at ligand idx")
        #
        # if idx_receptor == -1:
        #     idx_receptor = random.randit(1, 99)
        # else:
        #     if idx_receptor > 100:
        #         print("error at receptor idx")
        #
        # if idx_subunit == -1:
        #     idx_subunit = random.randit(1, 99)
        # else:
        #     if idx_subunit > 100:
        #         print("error at receptor idx")

        self.idx = idx
        self.idx_receptor = idx_receptor
        self.idx_ligand = idx_ligand
        self.idx_subunit = idx_subunit
        self.genotype = genotype
        self.score = score
        self.rmsd = rmsd
        self.i_sc = i_sc
        self.irms = irms
        self.receptor_name = receptor_name
        self.ligand_name = ligand_name
        self.subunit_name = subunit_name
