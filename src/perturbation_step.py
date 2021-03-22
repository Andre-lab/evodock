from pyrosetta import MoveMap, SwitchResidueTypeSetMover, create_score_function
from pyrosetta.rosetta import protocols
from pyrosetta.rosetta.protocols.rigid import (RigidBodyPerturbMover,
                                               RigidBodyRandomizeMover,
                                               RigidBodySpinMover,
                                               partner_downstream,
                                               partner_upstream)


def conversion_movers(pose):
    to_centroid = SwitchResidueTypeSetMover("centroid")
    to_fullatom = SwitchResidueTypeSetMover("fa_standard")
    recover_sidechains = protocols.simple_moves.ReturnSidechainMover(pose)
    return to_centroid, to_fullatom, recover_sidechains


class PerturbationStepMover:
    def initialize(self, **args):
        full_atom_docking = (
            args["full_atom_docking"] if "full_atom_docking" in args else False
        )
        pose = args["pose"]
        dock_jump = args["dock_jump"]
        translation = args["translation"]
        rotation = args["rotation"]
        to_fullatom = args["to_fullatom"]
        recover_sidechains = args["recover_sidechains"]
        # to_centroid = args["to_centroid"]

        randomize_upstream = RigidBodyRandomizeMover(pose, dock_jump, partner_upstream)
        randomize_downstream = RigidBodyRandomizeMover(
            pose, dock_jump, partner_downstream
        )
        dock_pert = RigidBodyPerturbMover(dock_jump, translation, rotation)
        # this Mover randomizes a pose's partners (rotation)
        spin = RigidBodySpinMover(dock_jump)
        # this Mover uses the axis defined by the inter-body jump (jump 1) to move
        #    the docking partners close together
        slide_into_contact = protocols.docking.DockingSlideIntoContact(dock_jump)

        # 8. setup the MinMover
        # the MoveMap can set jumps (by jump number) as degrees of freedom
        movemap = MoveMap()
        movemap.set_jump(dock_jump, True)
        # the MinMover can minimize score based on a jump degree of freedom, this
        #    will find the distance between the docking partners which minimizes
        #    the score
        if full_atom_docking:
            minmover = protocols.minimization_packing.MinMover()
            minmover.movemap(movemap)
            scorefxn_high_min = create_score_function("docking", "docking_min")
            minmover.score_function(scorefxn_high_min)

        # 9. create a SequenceMover for the perturbation step
        perturb = protocols.moves.SequenceMover()
        perturb.add_mover(randomize_upstream)
        perturb.add_mover(randomize_downstream)
        perturb.add_mover(dock_pert)
        perturb.add_mover(spin)
        perturb.add_mover(slide_into_contact)

        if full_atom_docking:
            ## looks like fullatom perturbation not works... I think it is not a problem
            perturb.add_mover(to_fullatom)
            perturb.add_mover(recover_sidechains)
            perturb.add_mover(minmover)

        return perturb
