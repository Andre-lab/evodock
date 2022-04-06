#!/usr/bin/env python
# coding: utf-8


from pyrosetta import Pose
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.core.pose import remove_virtual_residues

# from pyrosetta.rosetta.protocols.moves import PyMOLMover
# from src.utils import IP_ADDRESS


class RelaxOperator:
    def __init__(self, config, scfxn, local_search):
        # self.pymover = PyMOLMover(address=IP_ADDRESS, port=65000, max_packet_size=1400)
        self.relax = FastRelax(1)
        self.scfxn = scfxn
        self.relax.set_scorefxn(self.scfxn.scfxn_rosetta)
        self.scfxn_rosetta = self.scfxn.scfxn_rosetta
        self.config = config
        self.local_search = local_search.local_search_strategy

        self.jobid = self.config.jobid
        self.log_best = self.jobid + "/relax_best.log"
        if (
            self.config.docking_type_option == "Unbound"
            and self.config.bb_strategy == "relax_best"
        ):
            with open(self.log_best, "w") as file_object:
                file_object.write("#{}\n".format(self.jobid))

    def apply(self, input_pose):
        init_score = self.scfxn_rosetta.score(input_pose)
        # input_pose.pdb_info().name("rel_input_pose")
        # self.pymover.apply(input_pose)

        input_pose.dump_scored_pdb("best_to_relax.pdb", self.scfxn_rosetta)

        relax_pose = Pose()
        relax_pose.assign(input_pose)

        self.relax.apply(relax_pose)
        remove_virtual_residues(relax_pose)
        # relax_pose.pdb_info().name("rel_final_pose")
        # self.pymover.apply(relax_pose)
        final_score = self.scfxn_rosetta.score(relax_pose)

        if final_score < init_score:
            print("improved energy_score {} {}".format(init_score, final_score))
            self.scfxn.dock_pose.assign(relax_pose)
        else:
            print("NOT improved energy_score {} {}".format(init_score, final_score))
            self.scfxn.dock_pose.assign(input_pose)
        with open(self.log_best, "a") as file_object:
            file_object.write(
                "{} {} {}\n".format(init_score, final_score, final_score < init_score)
            )

        return self.scfxn.dock_pose, init_score, final_score
