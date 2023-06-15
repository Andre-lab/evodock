import random
from pyrosetta.rosetta.core.pose.symmetry import is_symmetric
import glob
import time
from src.genotype_converter import RefineCluspro
from src.utils import make_trial
from src.symmetry import individual_is_within_bounds
from cubicsym.actors.cubicsymmetryslider import InitCubicSymmetrySlider
from src.utils import get_position_info
from symmetryhandler.reference_kinematics import set_jumpdof_str_str
from pathlib import Path
from cubicsym.utilities import get_jumpidentifier, add_id_to_pose_w_base

class InitializePopulation:
    # Default initizalize Popualtion for Global docking
    def __init__(self, config, logger, popul_calculator, scfxn):
        self.config = config
        self.logger = logger
        self.popul_calculator = popul_calculator
        self.popsize = config.popsize
        # if config.syminfo:
        #     self.ind_size = config.syminfo.genotype_size
        # else:
        #     self.ind_size = 6
        # self.bounds = [(-1, 1)] * self.ind_size
        self.scfxn = scfxn
        self.init_slider = None
        if self.config.symmetric:
            if self.config.docking_type_option == "Global":
                raise NotImplementedError
            elif self.config.docking_type_option == "GlobalFromMultimer":
                self.init_slider = InitCubicSymmetrySlider(scfxn.dock_pose, self.config.syminfo.input_symdef, self.config.syminfo.ccsc,
                                                            pymolmover=self.config.pmm)

    def popul_is_within_bounds(self, popul):
        """Checks that the individuals, if symmetric, are wihtin bounds."""
        if self.config.symmetric:
            for ind in popul:
                individual_is_within_bounds(self.config, self.scfxn, ind)

    def attach_names_to_popul(self, popul):
        """Add the current names of the subunits/ligands/receptor to each individual"""
        if self.config.flexbb:
            for ind in popul:
                if self.config.syminfo:
                    ind.subunit_name = self.popul_calculator.local_search.local_search_strategy.swap_operator.list_subunits[ind.idx_subunit]
                else:
                    #FIXME: Not implemented
                    pass


    def set_genotype_from_init_bounds(self):
        starting_dofs = self.config.syminfo.get_position_info(self.scfxn.dock_pose)
        init_bounds = self.config.syminfo.init_bounds
        fixed = False
        # fix x trans and the center of mass rotations (4 dofs in total) if self.config.init_input_fix_percent
        if self.config.init_input_fix_percent is not None and random.uniform(0, 1) < self.config.init_input_fix_percent:
            fixed = True
            jid = get_jumpidentifier(self.scfxn.dock_pose)
            jumps_to_fix = [f'JUMP{jid}fold111', f'JUMP{jid}fold111_x', f'JUMP{jid}fold111_y', f'JUMP{jid}fold111_z']
            jump_order = self.config.syminfo.dof_spec.jump_str
            positions = [random.uniform(*pertubation) + dof if jumpname not in jumps_to_fix else dof for jumpname, pertubation, dof in zip(jump_order, init_bounds, starting_dofs)]
        else:
            positions = [random.uniform(*pertubation) + dof for pertubation, dof in zip(init_bounds, starting_dofs)]
        genotype = self.scfxn.convert_positions_to_genotype(positions)
        assert all([i <= 1 or i >= -1 for i in genotype])
        return genotype, fixed

    def set_genotype_from_bounds(self):
        bounds = [(-1, 1)] * self.scfxn.converter.size
        genotype = []
        for j in range(len(bounds)):
            genotype.append(random.uniform(bounds[j][0], bounds[j][1]))
        return genotype

    def flip_genotype_around_angle_x(self, genotype):
        index = self.config.syminfo.dof_spec.dof_str.index('angle_x')
        current_dofs = self.config.syminfo.get_position_info(self.scfxn.apply_genotype_to_pose(genotype))
        current_dofs[index] += 180
        genotype = self.scfxn.convert_positions_to_genotype(current_dofs)
        assert all([i <= 1 or i >= -1 for i in genotype])
        return genotype

    def make_symmetric_genotype(self, idx_subunit=None):
        """Makes a random uniformly distributed genotype around its bounds or init_bounds if set.
        If the docking_type is set to Global or GlobalFromMultimer, a sliding step will be carried out."""
        flipped, fixed = False, False
        if is_symmetric(self.scfxn.dock_pose):
            # if init bounds are set we need to construct it from those bounds, else we use the regular bounds.
            if self.config.syminfo.init_bounds is not None:
                genotype, fixed = self.set_genotype_from_init_bounds()
            else:
                genotype = self.set_genotype_from_bounds()
            # if GlobalFromMultimer and allow_flip is set, flip approx 50% of the individuals around their x-axis
            if self.config.docking_type_option == "GlobalFromMultimer":
                if self.config.allow_flip:
                    if random.uniform(0, 1) < 0.5:
                        genotype = self.flip_genotype_around_angle_x(genotype)
                        flipped = True
            # Apply slide if either Global or GlobalFromMultimer is set. It will change the genotype accordingly.
            if self.init_slider:
                assert idx_subunit is not None
                assert self.config.docking_type_option in ("GlobalFromMultimer", "Global")
                genotype = self.apply_initial_slide(genotype, idx_subunit)
        else:
            genotype = self.set_genotype_from_bounds()
        return genotype, flipped, fixed

    def make_heterodimeric_genotype(self):
        """Makes a random uniformly distributed genotype."""""
        return self.set_genotype_from_bounds()

    def apply_initial_slide(self, genotype, idx_subunit):
        """Applies slide if either Global or GlobalFromMultimer is set. It will change the genotype accordingly."""
        pose = self.scfxn.apply_genotype_to_pose(genotype)
        add_id_to_pose_w_base(pose, idx_subunit)
        max_slide_hit = self.init_slider.apply(pose)
        # merge the genotype with the new z, x values
        positions = get_position_info(pose, self.config.syminfo)
        self.logger.info(f" docked staring position: {', '.join(['{:.2f}'.format(p) for p in positions])}{' (MAXIMUM AMOUNT OF SLIDES HIT!)' if max_slide_hit else ''}")
        slided_genotype = self.scfxn.convert_positions_to_genotype(positions)
        return slided_genotype


    def init_population(self):
        """Initialize a population of individuals."""
        self.logger.info(" init population")
        popul = []
        for i in range(0, self.popsize):
            idx_receptor, idx_ligand, idx_subunit, flipped, fixed = None, None, None, None, None
            # symmetric genotype creation
            if is_symmetric(self.scfxn.dock_pose):
                if self.config.flexbb:
                    idx_subunit = random.randint(0, self.config.subunit_library_size - 1)
                genotype, flipped, fixed = self.make_symmetric_genotype(idx_subunit)
            # heterodimeric genotype creation
            else:
                if self.config.flexbb:
                    idx_receptor = random.randint(0, self.config.receptor_library_size - 1)
                    idx_ligand = random.randint(0, self.config.ligand_library_size - 1)
                genotype = self.make_heterodimeric_genotype()
            popul.append(make_trial(i, genotype, ligand=idx_ligand, receptor=idx_receptor, subunit=idx_subunit, flipped=flipped, fixed=fixed))
        start = time.time()
        self.popul_is_within_bounds(popul)
        popul = self.popul_calculator.run(popul, init_population=True)
        self.popul_is_within_bounds(popul)
        self.attach_names_to_popul(popul)
        end = time.time()
        self.logger.info(f" population init in {end - start:.2f} seconds")
        return popul

# class InitializePopulationLocal(InitializePopulation):
#     def init_population(self):
#         population_calculator = self.popul_calculator
#         popsize = self.popsize
#         if popsize is None:
#             popsize = self.popsize
#         self.logger.info(" init population")
#         # default values
#         popul = []
#         for i in range(0, popsize):
#             if is_symmetric(self.scfxn.dock_pose):
#                 # todo:
#                 #  generate_genotype just sets the local bounds (8 degrees rotation + 3 Å translation) for the dimeric case.
#                 #  for symmetry i still choose to use the symbounds set in the config file. I dont have default values as of yet: todo!
#                 indv = []
#                 if self.config.syminfo.init_bounds:
#                     bounds = self.config.syminfo.init_bounds
#                 else:
#                     bounds = self.bounds
#                 for j in range(len(bounds)):
#                     indv.append(random.uniform(bounds[j][0], bounds[j][1]))
#             else:
#                 indv = generate_genotype(self.scfxn.dock_pose, self.max_translation)
#             popul.append(make_trial(i, indv))
#         init_population = True
#         start = time.time()
#         population = population_calculator.run(popul, init_population)
#         end = time.time()
#         self.logger.info(f" population init in {end - start:.2f} seconds")
#         return population

class InitializePopulationRefine(InitializePopulation):
    def init_population(self):
        population_calculator = self.popul_calculator
        popsize = self.popsize
        refCluspro = RefineCluspro(self.config, self.config.get_max_translation(None))

        if popsize is None:
            popsize = self.popsize
        self.logger.info(" init population")
        popul = []
        for i in range(0, popsize):
            indv = refCluspro.refine_cluspro(i)
            popul.append(make_trial(i, indv))

        init_population = True
        start = time.time()
        population = population_calculator.run(popul, init_population)
        end = time.time()
        self.logger.info(f" population init in {end - start:.2f} seconds")
        return population