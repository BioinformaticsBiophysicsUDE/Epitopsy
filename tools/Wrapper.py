'''
Created on Feb 24, 2012

@author: chris
'''

import os
import sys
import csv
import shutil
import random
import subprocess
import numpy as np

from epitopsy.Structure import PDBFile



class Gromacs(object):

    def __init__(self, pdb_path, production_time_ns, temp, md_data_dir,
                 force_field = None,
                 water_model = None,
                 time_step = 0.002,
                 em_minimization_steps = 1000,
                 nvt_equilibration_time_ps = 200,
                 npt_equilibration_time_ps = 200):
        '''
        This class performs a md simulation.
        '''
        if force_field is None:
            # Are Protein Force Fields Getting Better:
            # A Systematic Benchmark on 524 Diverse NMR Measurements
            self.force_field = 'amber99sb-ildn'
            self.possible_modified_ff = False
        else:
            self.force_field = force_field
            self.possible_modified_ff = True

        if water_model is None:
            self.water_model = 'tip3p'
        else:
            self.water_model = water_model

        # box parameters
        self.box_type = 'dodecahedron'
        self.box_dimension = '1.0'

        # solvation parameters
        self.water_coordinates = 'spc216.gro'

        # ion parameters
        self.counter_neutral_conc = '0.15'

        # Temperature
        self.temp = temp
        # path to the structure
        self.pdb_path = pdb_path

        # timing
        self.time_step = time_step
        # time for nvt equilibration
        self.nvt_equilibration_time_ps = float(nvt_equilibration_time_ps)
        self.nvt_nsteps = int(self.nvt_equilibration_time_ps * 10 ** (-12) / (10 ** (-12)) / self.time_step)
        # time for npt equilibration
        self.npt_equilibration_time_ps = float(npt_equilibration_time_ps)
        self.npt_nsteps = int(self.npt_equilibration_time_ps * 10 ** (-12) / (10 ** (-12)) / self.time_step)
        self.em_minimization_steps = em_minimization_steps

        # time for the production
        self.production_time_ns = float(production_time_ns)
        self.md_nsteps = int(self.production_time_ns * 10 ** (-9) / (10 ** (-12)) / self.time_step)

        #starting sequence
        self.input_structure = self.pdb_path
        self.pdb_name = self.pdb_path.replace('.pdb', '')
        self.md_data_dir = md_data_dir


        # mdp files
        self.em_mdp = 'em.mdp'
        self.nvt_mdp = 'nvt.mdp'
        self.npt_mdp = 'npt.mdp'
        self.md_mdp = 'md_prod.mdp'

        # gromacs files
        # pdb2gmx
        self.gro_processed = '{0}_processed.gro'.format(self.pdb_name)
        # editconf -> box
        self.gro_new_box = '{0}_new_box.gro'.format(self.pdb_name)
        # genbox -> solvate
        self.gro_solvated = '{0}_solvated.gro'.format(self.pdb_name)
        # genion -> add ions
        self.gro_solvated_ions = '{0}_solvated_ions.gro'.format(self.pdb_name)
        # grompp -> em
        self.gro_em = 'em.gro'
        # grompp -> nvt
        self.gro_nvt = 'nvt.gro'
        # grompp -> npt
        self.gro_npt = 'npt.gro'

        self.topol_top = 'topol.top'
        self.ions_tpr = 'ions.tpr'
        self.em_tpr = 'em.tpr'
        self.nvt_tpr = 'nvt.tpr'
        self.npt_tpr = 'npt.tpr'
        self.md_prod_tpr = 'md_prod.tpr'

        self.md_xtc = 'md_prod.xtc'
        self.md_no_pbc_xtc = 'md_prod_no_pbc.xtc'

        # required files
        self.req_file_list = [self.input_structure, self.em_mdp,
                              self.nvt_mdp, self.npt_mdp, self.md_mdp]


    def run(self):
        '''Set up the md, run equilibration and run the md.
        '''
        self._run_em_minimization()
        self._run_production()

    def _run_em_minimization(self):
        if os.path.exists(self.md_data_dir):
            raise AttributeError("'{0}' path already exists!".format(self.md_data_dir))
        else:
            os.mkdir(self.md_data_dir)
            shutil.copy(self.input_structure,
                        os.path.join(self.md_data_dir, self.input_structure))

            # check if there might be a modified forcefield
            if self.possible_modified_ff is True:
                ff_dir = "{0}.ff".format(self.force_field)
                if os.path.exists(ff_dir):
                    new_ff_dir = os.path.join(self.md_data_dir, ff_dir)
                    os.mkdir(new_ff_dir)
                    for item in os.listdir(ff_dir):
                        shutil.copy(os.path.join(ff_dir, item),
                                os.path.join(new_ff_dir, item))


            self._current_dir = os.getcwd()
            os.chdir(self.md_data_dir)

        # write mdp files to disk
        self._write_mdp_files()

        #### Set Up Topology
        self._set_up_topology()

        #### Generate Box
        # Notice that the created box can not be deduced from the positions of
        # the waters in pymol after you have converted the *.gro file to a
        # *.pdb file. pymol visualizes a square box instead of a dodecahedron.
        self._generate_box()

        #### Solvate Molecules and Box
        self._solvate_box()

        #### Generate Genion Input File
        self._generate_tpr_input(self.em_mdp, self.gro_solvated,
                self.topol_top, self.ions_tpr)

        #### Counter Ions
        self._add_counter_ions()

        #### Generate MD Input File
        self._generate_tpr_input(self.em_mdp, self.gro_solvated_ions,
                self.topol_top, self.em_tpr)

        #### MD for Energy Minimisation
        self._run_md('em')

        #### analyze minimization
        process = subprocess.Popen(['g_energy_d', '-f', 'em.edr', '-o',
                                    'potential_data.xvg'],
                                   stdin = subprocess.PIPE)
        process.communicate('Potential')

    def _run_production(self):
        #### Bring system to temperature
        self._generate_tpr_input(self.nvt_mdp, self.gro_em, self.topol_top,
                self.nvt_tpr)

        #### MD for temperature equilibration
        self._run_md('nvt')

        #### analyze temperature
        process = subprocess.Popen(['g_energy_d', '-f', 'nvt.edr', '-o',
                                    'temperature_data.xvg'],
                                   stdin = subprocess.PIPE)
        process.communicate('Temperature')

        #### Bring system to pressure
        self._generate_tpr_input(self.npt_mdp, self.gro_nvt, self.topol_top,
                self.npt_tpr)

        #### MD for pressure equilibration
        self._run_md('npt')

        #### analyze pressure
        process = subprocess.Popen(['g_energy_d', '-f', 'npt.edr', '-o',
                                    'pressure_data.xvg'],
                                   stdin = subprocess.PIPE)
        process.communicate('Pressure')

        #### analyze density
        process = subprocess.Popen(['g_energy_d', '-f', 'npt.edr', '-o',
                                    'density_data.xvg'],
                                   stdin = subprocess.PIPE)
        process.communicate('Density')

        #### Generate MD Input File
        self._generate_tpr_input(self.md_mdp, self.gro_npt, self.topol_top,
                self.md_prod_tpr)

        #### MD Simulation
        self._run_md('md_prod')

        #### correct for pbc
        process = subprocess.Popen(['trjconv_d', '-s', self.md_prod_tpr, '-f', self.md_xtc,
                                    '-o', self.md_no_pbc_xtc, '-pbc', 'mol', '-ur',
                                    'compact'], stdin = subprocess.PIPE)
        process.communicate('System')

        #### Analysis -> rmsd
        process = subprocess.Popen(['g_rms_d', '-s', self.md_prod_tpr, '-f',
                                    self.md_no_pbc_xtc, '-tu', 'ns', '-o',
                                    'rmsd_backbone.xvg'],
                                   stdin = subprocess.PIPE)
        process.stdin.write('Backbone\n')
        process.communicate('Backbone')

        #### Analysis -> rmsf res
        process = subprocess.Popen(['g_rmsf_d', '-s', self.md_prod_tpr, '-f',
                                    self.md_no_pbc_xtc, '-res', '-o',
                                    'rmsf_residue.xvg'],
                                   stdin = subprocess.PIPE)
        process.communicate('Protein')
        #### Analysis -> rmsf calpha
        process = subprocess.Popen(['g_rmsf_d', '-s', self.md_prod_tpr, '-f',
                                    self.md_no_pbc_xtc, '-o', 'rmsf_calpha.xvg'],
                                   stdin = subprocess.PIPE)
        process.communicate('C-alpha')
        #### Analysis -> radius of gyration
        process = subprocess.Popen(['g_gyrate_d', '-s', self.md_prod_tpr, '-f',
                                    self.md_no_pbc_xtc, '-o', 'gyrate.xvg'],
                                   stdin = subprocess.PIPE)
        process.communicate('Protein')


        os.chdir(self._current_dir)

    def _set_up_topology(self):
        subprocess.call(['pdb2gmx_d', '-ff', self.force_field, '-f',
                         self.input_structure, '-o', self.gro_processed, '-water',
                         self.water_model, '-ignh'])

    def _generate_box(self):
        # Notice that the created box can not be deduced from the positions of
        # the waters in pymol after you have converted the *.gro file to a
        # *.pdb file. pymol visualizes a square box instead of a dodecahedron.
        subprocess.call(['editconf_d', '-bt', self.box_type, '-f',
                         self.gro_processed, '-o', self.gro_new_box, '-d',
                         self.box_dimension, '-c'])

    def _solvate_box(self):
        subprocess.call(['genbox_d', '-cp', self.gro_new_box, '-cs',
                         self.water_coordinates, '-o', self.gro_solvated, '-p',
                         self.topol_top])

    def _add_counter_ions(self):
        counter_log = 'ion.log'
        process = subprocess.Popen(['genion_d', '-s', self.ions_tpr, '-o',
                                    self.gro_solvated_ions, '-neutral', '-conc',
                                    self.counter_neutral_conc, '-p',
                                    self.topol_top, '-g', counter_log],
                                   stdin = subprocess.PIPE)
        process.communicate('SOL')
        #
        # or:
        # process.stdin.write('SOL')
        # process.communicate()

    def _generate_tpr_input(self, mdp_file, gro_file, topol_file,
            tpr_output):
        subprocess.call(['grompp_d', '-f', mdp_file, '-c', gro_file, '-p',
                         topol_file, '-o', tpr_output])

    def _run_md(self, id_name):
        subprocess.call(['mdrun_d', '-deffnm', id_name])

    def _write_mdp_files(self):
        with open(self.em_mdp, 'w') as f:
            f.write("""; em.mdp - used as input into grompp to generate em.tpr
; Parameters describing what to do, when to stop and what to save
integrator    = steep        ; Algorithm (steep = steepest descent minimization)
emtol        = 500.0      ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
emstep          = 0.01          ; Energy step size
nsteps        = {0}          ; Maximum number of (minimization) steps to perform

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist        = 1        ; Frequency to update the neighbor list and long range forces
ns_type        = grid        ; Method to determine neighbor list (simple, grid)
rlist        = 1.0        ; Cut-off for making neighbor list (short range forces)
coulombtype    = PME        ; Treatment of long range electrostatic interactions
rcoulomb    = 1.0        ; Short-range electrostatic cut-off
rvdw        = 1.0        ; Short-range Van der Waals cut-off
pbc        = xyz         ; Periodic Boundary Conditions (yes/no)""".format(
                self.em_minimization_steps))

        with open(self.nvt_mdp, 'w') as f:
            f.write("""title        = NVT equilibration
define        = -DPOSRES    ; position restrain the protein
; Run parameters
integrator    = md        ; leap-frog integrator
nsteps        = {0}    ; {1} * {0} fs = {2} ps
dt        = {3}        ; {1} fs
; Output control
nstxout        = 100        ; save coordinates every 0.2 ps
nstvout        = 100        ; save velocities every 0.2 ps
nstenergy    = 100        ; save energies every 0.2 ps
nstlog        = 100        ; update log file every 0.2 ps
; Bond parameters
continuation    = no        ; first dynamics run
constraint_algorithm = lincs    ; holonomic constraints
constraints    = all-bonds    ; all bonds (even heavy atom-H bonds) constrained
lincs_iter    = 1        ; accuracy of LINCS
lincs_order    = 4        ; also related to accuracy
; Neighborsearching
ns_type        = grid        ; search neighboring grid cells
nstlist        = 5        ; 10 fs
rlist        = 1.0        ; short-range neighborlist cutoff (in nm)
rcoulomb    = 1.0        ; short-range electrostatic cutoff (in nm)
rvdw        = 1.0        ; short-range van der Waals cutoff (in nm)
; Electrostatics
coulombtype    = PME        ; Particle Mesh Ewald for long-range electrostatics
pme_order    = 4        ; cubic interpolation
fourierspacing    = 0.16        ; grid spacing for FFT
; Temperature coupling is on
tcoupl        = V-rescale    ; modified Berendsen thermostat
tc-grps        = Protein Non-Protein    ; two coupling groups - more accurate
tau_t        = 0.1    0.1    ; time constant, in ps
ref_t        = {4}     {4}    ; reference temperature, one for each group, in K
; Pressure coupling is off
pcoupl        = no         ; no pressure coupling in NVT
; Periodic boundary conditions
pbc        = xyz        ; 3-D PBC
; Dispersion correction
DispCorr    = EnerPres    ; account for cut-off vdW scheme
; Velocity generation
gen_vel        = yes        ; assign velocities from Maxwell distribution
gen_temp    = {4}        ; temperature for Maxwell distribution
gen_seed    = -1        ; generate a random seed""".format(self.nvt_nsteps,
                                                           self.time_step * 10 ** 3,
                                                           self.nvt_equilibration_time_ps,
                                                           self.time_step,
                                                           self.temp))

        with open(self.npt_mdp, 'w') as f:
            f.write("""title        = NPT equilibration
define        = -DPOSRES    ; position restrain the protein
; Run parameters
integrator    = md        ; leap-frog integrator
nsteps        = {0}    ; {1} * {0} fs = {2} ps
dt        = {3}        ; {1} fs
; Output control
nstxout        = 100        ; save coordinates every 0.2 ps
nstvout        = 100        ; save velocities every 0.2 ps
nstenergy    = 100        ; save energies every 0.2 ps
nstlog        = 100        ; update log file every 0.2 ps
; Bond parameters
continuation    = yes        ; Restarting after NVT
constraint_algorithm = lincs    ; holonomic constraints
constraints    = all-bonds    ; all bonds (even heavy atom-H bonds) constrained
lincs_iter    = 1        ; accuracy of LINCS
lincs_order    = 4        ; also related to accuracy
; Neighborsearching
ns_type        = grid        ; search neighboring grid cells
nstlist        = 5        ; 10 fs
rlist        = 1.0        ; short-range neighborlist cutoff (in nm)
rcoulomb    = 1.0        ; short-range electrostatic cutoff (in nm)
rvdw        = 1.0        ; short-range van der Waals cutoff (in nm)
; Electrostatics
coulombtype    = PME        ; Particle Mesh Ewald for long-range electrostatics
pme_order    = 4        ; cubic interpolation
fourierspacing    = 0.16        ; grid spacing for FFT
; Temperature coupling is on
tcoupl        = V-rescale    ; modified Berendsen thermostat
tc-grps        = Protein Non-Protein    ; two coupling groups - more accurate
tau_t        = 0.1    0.1    ; time constant, in ps
ref_t        = {4}     {4}    ; reference temperature, one for each group, in K
; Pressure coupling is on
pcoupl        = Parrinello-Rahman    ; Pressure coupling on in NPT
pcoupltype    = isotropic    ; uniform scaling of box vectors
tau_p        = 2.0        ; time constant, in ps
ref_p        = 1.0        ; reference pressure, in bar
compressibility = 4.5e-5    ; isothermal compressibility of water, bar^-1
refcoord_scaling = com
; Periodic boundary conditions
pbc        = xyz        ; 3-D PBC
; Dispersion correction
DispCorr    = EnerPres    ; account for cut-off vdW scheme
; Velocity generation
gen_vel        = no        ; Velocity generation is off""".format(self.npt_nsteps,
                                                                   self.time_step * 10 ** 3,
                                                                   self.npt_equilibration_time_ps,
                                                                   self.time_step,
                                                                   self.temp))

        with open(self.md_mdp, 'w') as f:
            f.write("""title        = MD
; Run parameters
integrator    = md        ; leap-frog integrator
nsteps        = {0}    ; {1} * {0} fs = {2} ns
dt        = {3}        ; {1} fs
; Output control
nstxout        = 500        ; save coordinates every 1 ps
nstvout        = 500        ; save velocities every 2 ps
nstxtcout    = 500        ; xtc compressed trajectory output every 1 ps
nstenergy    = 500        ; save energies every 1 ps
nstlog        = 1000        ; update log file every 2 ps
; Bond parameters
continuation    = yes        ; Restarting after NPT
constraint_algorithm = lincs    ; holonomic constraints
constraints    = all-bonds    ; all bonds (even heavy atom-H bonds) constrained
lincs_iter    = 1        ; accuracy of LINCS
lincs_order    = 4        ; also related to accuracy
; Neighborsearching
ns_type        = grid        ; search neighboring grid cells
nstlist        = 5        ; 10 fs
rlist        = 1.0        ; short-range neighborlist cutoff (in nm)
rcoulomb    = 1.0        ; short-range electrostatic cutoff (in nm)
rvdw        = 1.0        ; short-range van der Waals cutoff (in nm)
; Electrostatics
coulombtype    = PME        ; Particle Mesh Ewald for long-range electrostatics
pme_order    = 4        ; cubic interpolation
fourierspacing    = 0.16        ; grid spacing for FFT
; Temperature coupling is on
tcoupl        = V-rescale    ; modified Berendsen thermostat
tc-grps        = Protein Non-Protein    ; two coupling groups - more accurate
tau_t        = 0.1    0.1    ; time constant, in ps
ref_t        = {4}     {4}    ; reference temperature, one for each group, in K
; Pressure coupling is on
pcoupl        = Parrinello-Rahman    ; Pressure coupling on in NPT
pcoupltype    = isotropic    ; uniform scaling of box vectors
tau_p        = 2.0        ; time constant, in ps
ref_p        = 1.0        ; reference pressure, in bar
compressibility = 4.5e-5    ; isothermal compressibility of water, bar^-1
; Periodic boundary conditions
pbc        = xyz        ; 3-D PBC
; Dispersion correction
DispCorr    = EnerPres    ; account for cut-off vdW scheme
; Velocity generation
gen_vel        = no        ; Velocity generation is off """.format(self.md_nsteps,
                                                                   self.time_step * 10 ** 3,
                                                                   self.production_time_ns,
                                                                   self.time_step,
                                                                   self.temp))



class LatPack:

    def __init__(self, seq, fold, outFile, kT, step_list,
                 moveSet_list, runs, energy = None, seed = None, lat = 'FCC',
                 energyFile = None, out = 'S', lat_pack_path = None,
                 is_trajectory = True, traj_path = 'traj.data', conv_mode = 'a2x'):
        '''
        If a trajectory is calculated an already existing file of the same
        name is deleted.
        '''
        # check if step_list and moveSet_list match
        if len(step_list) != len(moveSet_list):
            raise AttributeError("Missmatch between length of step_list '{0}' and moveSet_list '{1}'!".format(len(step_list), len(moveSet_list)))

        # check which kind of moves
        for item in moveSet_list:
            if item == 'PivotM' or item == 'PullM':
                continue
            else:
                raise ValueError('Unknown move set: {0}'.format(item))

        if lat_pack_path is None:
            self.lat_fold_path = os.path.join(os.getenv('HOME'), 'Programs',
                                         'tmpLatpackminE_minFold', 'bin',
                                         'latFold')
            self.lat_conv_path = os.path.join(os.getenv('HOME'), 'Programs',
                                         'tmpLatpackminE_minFold', 'bin',
                                         'latConv')
        else:
            self.lat_fold_path = os.path.join(lat_pack_path, 'bin', 'latFold')
            self.lat_conv_path = os.path.join(lat_pack_path, 'bin', 'latConv')

        if energyFile is None:
            energyFile = os.path.join(os.getenv('HOME'), 'Programs',
                                      'tmpLatpackminE_minFold', 'energy_files',
                                      'MJ')


        if not os.path.exists(self.lat_fold_path):
            raise AttributeError("Could not find program latFold at '{0}'!".format(self.lat_fold_path))

        # sequence
        self.seq = seq
        # absolute move string == fold
        self.fold = fold
        # energyFile (MJ)
        self.energyFile = energyFile
        # energy of the current fold
        self.energy = energy
        # folding temperature
        self.kT = kT
        # list of maximal number of folding steps
        self.step_list = step_list
        # number of seeds / folding simulations
        self.runs = runs
        # lattice used
        self.lat = lat
        # list of applied move sets
        self.moveSet_list = moveSet_list
        # output mode
        self.out = out
        # output file
        self.outFile = outFile
        # is there a trajectory?
        self.is_trajectory = is_trajectory
        # path for the trajectory
        self.traj_path = traj_path
        if self.is_trajectory and self.traj_path is not None:
            if os.path.exists(self.traj_path):
                os.remove(self.traj_path)

        # seed option
        self.seed = seed

        # conversion mode
        self.conv_mode = conv_mode

        if self.energy is None:
            # calculate the energy for the fold by performing 0 steps
            cmd = self._get_fold_command(self.fold, 0, 'PivotM')
            self.energy, fold = self._run_fold(cmd)
            # this should not happen!
            # except ones starts from an extended chain
            if self.fold is not None and fold != self.fold:
                raise AttributeError('Fold mismatch after scoring the inital fold!\n{0}\n!=\n{1}'.format(self.fold, fold))


    def run_fold(self):
        '''
        Run the program.
        '''
        E_best = 0
        fold_best = self.fold
        new_fold = self.fold
        for step_size, moveSet in zip(self.step_list, self.moveSet_list):
            # get the commant to run the program
            command = self._get_fold_command(new_fold, step_size, moveSet)
            # get energy and fold
            energy, new_fold = self._run_fold(command)

#            # check if we have found a better fold
#            if energy < E_best:
#                E_best = energy
#                fold_best = new_fold

        return [energy, new_fold]
#        return [E_best, fold_best]

    def _run_fold(self, command):
        # run the command
        p = subprocess.call(command, stdout = subprocess.PIPE,
                            stderr = subprocess.PIPE)
        if p != 0:
            raise AttributeError("latFold error for command:\n{0}".format(' '.join(command)))

        if self.is_trajectory is True:
            with open(self.outFile) as f:
                content = f.readlines()
            with open(self.traj_path, 'a') as f:
                for line in content:
                    f.write(line)

            energy = float(content[-1].split(' ')[0])
            new_fold = content[-1].split(' ')[1].rstrip('\n')
        else:
            print('not okay')
            # read calculated energy and fold
            with open(self.outFile) as f:
                for i, line in enumerate(f):
                    if i == 0:
                        energy = float(line.split(' ')[0])
                        new_fold = line.split(' ')[1].rstrip('\n')

        # remove the out file
        os.remove(self.outFile)

        return energy, new_fold

    def get_rmsd(self, fold_1, seq_1, fold_2, seq_2, pdb_path):
        '''
        This method converts the first structure and saves it as pdb_path,
        then it reads the data and converts the second and overrides
        pdb_path. After the calculation the rmsd is returned and the pdb_path
        will be removed.
        The second structure will be rotated
        '''
        # convert first structure
        self.convert_fold(fold_1, seq_1, pdb_path)
        pdb_1 = PDBFile(pdb_path)
        # convert second structure
        self.convert_fold(fold_2, seq_2, pdb_path)
        pdb_2 = PDBFile(pdb_path)
        # remove the pdb file
        os.remove(pdb_path)
        rmsd = pdb_1.get_rmsd_rotation_translation_from_superposition(pdb_2)['rmsd']
        return rmsd

    def convert_fold(self, fold, seq, pdb_file):
        cmd = self._get_conv_command(fold, seq)
        with open(pdb_file, 'w') as f:
            subprocess.call(cmd, stdout = f, stderr = subprocess.PIPE)

    def get_energy_and_fold_from_trajectory(self):
        '''
        Returns:
            A list with the energy as the first and the most populated
            fold as the second item.
        '''
        with open(self.traj_path) as f:
            content = f.readlines()

        traj_start = int(round(len(content) / 2))
        fold_dict = {}

        # set starting fold and energy
        pop_fold = content[traj_start].split(' ')[1].strip()
        pop_energy = float(content[traj_start].split(' ')[0])

        # count mc steps for the mean energy
        mc_counter = 0.
        # add all energies
        mc_energy = 0

        for line in content[traj_start:]:
            mc_counter += 1.
            energy = float(line.split(' ')[0])
            mc_energy = mc_energy + energy

            fold = line.split(' ')[1].strip()

            # check if this is a new fold or not
            if fold not in fold_dict:
                fold_dict[fold] = 1
            else:
                fold_dict[fold] += 1

            if fold_dict[fold] > fold_dict[pop_fold]:
                # the new fold has a higher population
                pop_energy = energy
                pop_fold = fold
            elif fold_dict[fold] == fold_dict[pop_fold]:
                # the new fold is equal to the old best one ...
                if energy < pop_energy:
                    # ... but has a lower energy
                    pop_energy = energy
                    pop_fold = fold

        mean_energy = mc_energy / mc_counter
        return [mean_energy, pop_fold]


    def get_min_energy_and_fold_from_trajectory(self):
        '''
        Returns:
            A list with the minimum energy in the trajectory as the first and the corresponding
            fold as the second item.
        '''

        energy_dict = {}

        with open(self.traj_path) as f:
            content = f.readlines()

        min_energy = 0
        energy_dict[min_energy] = []
        for line in content:
            energy = float(line.split(' ')[0])
            fold = line.split(' ')[1].strip()

            if energy < min_energy:
                min_energy = energy
                energy_dict[min_energy] = [fold]
            elif energy == min_energy and fold not in energy_dict[min_energy]:
                energy_dict[min_energy].append(fold)

        return [min_energy, energy_dict[min_energy]]


    def _get_fold_command(self, fold, step_size, moveSet):
        '''
        generate latFold run command
        '''

        # command is a list containing all the parameters (obligatory for subprocess in Evaluator)
        command = []

        command.append(self.lat_fold_path)

        # append sequence to fold
        command.append('-seq={0}'.format(self.seq))
        if fold is not None:
            # append absolute move path
            command.append('-abs={0}'.format(fold))
        # append energy file path
        command.append('-energyFile={0}'.format(self.energyFile))
        # append kT value
        command.append('-kT={0:.2f}'.format(self.kT))
        # append maximum number of steps to perform
        command.append('-maxSteps={0:d}'.format(step_size))
        # append number of runs to perform
        command.append('-runs={0:d}'.format(self.runs))
        # append lattice
        command.append('-lat={0}'.format(self.lat))
        # append moveSet
        command.append('-moveSet={0}'.format(moveSet))
        # append output mode
        command.append('-out={0}'.format(self.out))
        # append output file name
        command.append('-outFile={0}'.format(self.outFile))
        # append seed option
        if self.seed is None:
            seed = random.randint(1, 2147483647)
        else:
            seed = self.seed
        command.append('-seed={0:d}'.format(seed))
        # append verbose option
        #self.command.append('-v')

        return command

    def _get_conv_command(self, fold, seq):
        command = []
        command.append(self.lat_conv_path)
        command.append('-m={0}'.format(self.conv_mode))
        command.append('-str={0}'.format(fold))
        command.append('-seq={0}'.format(seq))
        command.append('-lat={0}'.format(self.lat))
        return command



class Rosetta(object):
    '''
    classdocs
    '''


    def __init__(self, rosetta_executables_path = None,
                 rosetta_database_path = None):
        '''

        '''
        if rosetta_executables_path is None:
            self.rosetta_executables_path = '/home/chris/tmp/rosetta3.4/rosetta_source/bin/'
        else:
            self.rosetta_executables_path = rosetta_executables_path

        if rosetta_database_path is None:
            self.rosetta_database_path = '/home/chris/tmp/rosetta3.4/rosetta_database'
        else:
            self.rosetta_database_path = rosetta_database_path

    def score_pdb(self, pdb_path_list, energy_func = 'standard'):
        '''
        This functions scores the given pdbs and returns their score as a
        list. The returned list should contain the energies in the same order
        as the given pdb_path_list
        '''
        if not isinstance(pdb_path_list, list):
            raise ValueError("Given object '{0}' is not a list!".format(pdb_path_list))

        list_path = 'list_to_score.txt'
        with open(list_path, 'w') as f:
            for item in pdb_path_list:
                f.write('{0}\n'.format(item))

        executable = 'score.linuxgccrelease'
        cmd = [os.path.join(self.rosetta_executables_path, executable),
               '-database {0}'.format(self.rosetta_database_path),
               '-in:file:l {0}'.format(list_path),
               '-score:weights {0}'.format(energy_func),
               '-rescore:verbose']
        p = subprocess.Popen(cmd, stdout = subprocess.PIPE,
                             stderr = subprocess.PIPE)
        rosetta_output, rosetta_error = p.communicate()

        rosetta_output = rosetta_output.split('\n')

        energy_list = []
        for line in rosetta_output:
            if line.startswith(' Total weighted score:'):
                line_content = line.split(':')
                energy_list.append(float(line_content[-1]))


        # remove pointless energy list -> it is without any names, so the
        # calculated energies can not be matched to their structures ...
        os.remove('default.sc')

        return energy_list


    def refine_pdb(self, pdb_path, new_path):
        '''
        This function calls rosetta to minimize the given pdb. The input and
        the output are superimposed afterwards.
        '''
        executable = 'relax.linuxgccrelease'
        wt_pdb = PDBFile(pdb_path)
        cmd = [os.path.join(self.rosetta_executables_path, executable),
               '-database {0}'.format(self.rosetta_database_path),
               '-in:file:s {0}'.format(pdb_path),
               '-relax:thorough',
               '-mute all']
        p = subprocess.Popen(cmd, stdout = subprocess.PIPE,
                             stderr = subprocess.PIPE)
        rosetta_output, rosetta_error = p.communicate()

        rosetta_output = rosetta_output.split('\n')

        # rename minimized structure
        shutil.move('{0}_0001.pdb'.format(pdb_path.replace('.pdb', '')), new_path)
        # superimpose structures
        refined_pdb = PDBFile(new_path)
        refined_pdb.superimpose_self_onto_given_pdb(wt_pdb, atom_types = 'CA')
        refined_pdb.save_to_file(new_path)
        os.remove('score.sc')

    def calculate_ddG_monomer(self, refined_wt_path, new_sequence,
                              new_path_prefix = None, iterations = 50,
                              energy_func = 'standard'):
        '''
        This function calls rosetta to mutate a given pdb. It is recommended to
        refine the structure before mutating it. The wild type should be copied
        to the 'new_dir' previously. The calculated ddG will be returned and is
        calculated like this:
            ddG = dG_mut - dG_wt
        and is < 0 for stabilizing mutations and > 0 for destabilizing ones.
        If no 'new_path_prefix' is given, it will delete the mutant pdbs,
        otherwise they will be named like this:
            '{0}_{1}.pdb'.format(new_path_prefix, i) with i in range(iterations).
        '''
        # check for different chain lengths
        wt_pdb = PDBFile(refined_wt_path)
        old_seq = wt_pdb.get_sequence().values()[0]

        if len(old_seq) != len(new_sequence):
            raise ValueError('Sequence lengths do not match! old {0} vs new {1}!'.format(len(old_seq), len(new_sequence)))

        # get mut list from the old and the new sequence
        # -> mutation instructions
        mut_list = []
        pos_counter = 0
        for i in range(len(old_seq)):
            if old_seq[pos_counter] != new_sequence[pos_counter]:
                mut_list.append('{0}{1}{2}'.format(old_seq[pos_counter],
                                                   pos_counter + 1, # rosetta starts at 1
                                                   new_sequence[pos_counter]))

            pos_counter = pos_counter + 1

        if len(mut_list) is 0:
            print('Given sequence is the same as the old one! Returning 0.')
            return 0.0

        # rosetta name
        rosetta_name = ''.join(mut_list)

        # file paths
        resfile = 'resfile.txt'
        ddg_file_path = 'ddg_predictions.out'

        # load wt_pdb for pair fitting
        wt_pdb = PDBFile(refined_wt_path)

        # prepare resfile:
        num_mutation = len(mut_list)
        with open(resfile, 'w') as f:
            f.write('total {0}\n'.format(num_mutation))
            f.write('{0}\n'.format(num_mutation))
            for item in mut_list:
                f.write('{0} {1} {2}\n'.format(item[0], int(item[1:-1]),
                                               item[-1]))

        executable = 'ddg_monomer.linuxgccrelease'
        cmd = [os.path.join(self.rosetta_executables_path,
                            executable),
                '-database {0}'.format(self.rosetta_database_path),
                '-in:file:s {0}'.format(refined_wt_path),
                '-in:file:fullatom',
                '-ignore_unrecognized_res',
                '-ddg:mut_file {0}'.format(resfile),
                '-ddg:iterations {0}'.format(iterations),
                '-ddg:dump_pdbs true',
                '-ddg:local_opt_only false',
                '-ddg:min_cst false',
                '-ddg:suppress_checkpointing true',
                '-ddg:mean false',
                '-ddg:min true',
                '-ddg:sc_min_only false',
                '-ddg:ramp_repulsive true',
                '-ddg:weight_file soft_rep_design',
                #'-mute all',
                '-ddg:minimization_scorefunction score12']
        p = subprocess.Popen(cmd, stdout = subprocess.PIPE,
                             stderr = subprocess.PIPE)

        rosetta_output, rosetta_error = p.communicate()

        rosetta_output = rosetta_output.split('\n')

        # iterate to find all scores
        wt_score = []
        wt_counter = 1 # rosetta starts with 1
        mut_score = []
        mut_counter = 1 # rosetta starts with 1
        for line in rosetta_output:
            wt_start_string = 'protocols.moves.ddGMover: {0} score before mutation: residue'.format(wt_counter)
            mut_start_string = 'protocols.moves.ddGMover: round {0} mutate {1}'.format(mut_counter, rosetta_name)
            if line.startswith(wt_start_string):
                score = float(line.replace(wt_start_string, '').strip().split(' ')[0])
                wt_score.append(score)
                wt_counter = wt_counter + 1
            elif line.startswith(mut_start_string):
                score = float(line.replace(mut_start_string, ''))
                mut_score.append(score)
                mut_counter = mut_counter + 1

        wt_score = np.array(wt_score)
        mut_score = np.array(mut_score)

        self.wt_score = wt_score
        self.mut_score = mut_score

        # read ddG, ddG_rosetta = wt - mut -> ddG = -ddG_rosetta or is it not?
        with open(ddg_file_path, 'r') as f:
            for line in f:
                if line.startswith('ddG: {0}'.format(rosetta_name)):
                    line_content = line.split('   ')
                    ddG_rosetta = float(line_content[1])

        # get list of pdb paths
        # could be done in one step, but it is messy ...
        my_wt_list = []
        my_mut_list = []
        for i in range(1, iterations + 1):
            my_wt_list.append('repacked_wt_round_{0}.pdb'.format(i))
            my_mut_list.append('mut_{0}_round_{1}.pdb'.format(rosetta_name, i))

        my_wt_score = np.array(self.score_pdb(my_wt_list, energy_func))
        my_mut_score = np.array(self.score_pdb(my_mut_list, energy_func))

        # rename/delete mutant pdbs and clean up
        os.remove('wt_traj')
        os.remove('mutant_traj{0}'.format(rosetta_name))
        os.remove(ddg_file_path)
        os.remove(resfile)
        for i in range(1, iterations + 1):
            os.remove('repacked_wt_round_{0}.pdb'.format(i))
            if new_path_prefix is None:
                os.remove('mut_{0}_round_{1}.pdb'.format(rosetta_name, i))
            else:
                new_path = '{0}_{1}.pdb'.format(new_path_prefix, i)
                shutil.move('mut_{0}_round_{1}.pdb'.format(rosetta_name, i),
                            new_path)
                mut_pdb = PDBFile(new_path)
                mut_pdb.superimpose_self_onto_given_pdb(wt_pdb, atom_types = 'CA')
                mut_pdb.save_to_file(new_path)

        # check if everything is alright
        if len(mut_score) != len(wt_score):
            raise ValueError('There has been an error with the identification of the scored structures!')

        if len(wt_score) != iterations:
            raise ValueError('Not all scores of the created structures could be retrieved!')

        ddG_min = np.amin(mut_score - wt_score)
        ddG_max = np.amax(mut_score - wt_score)
        ddG_max_diff = np.amin(mut_score) - np.amin(wt_score)
        ddG_mean = np.mean(mut_score - wt_score)
        my_ddG = np.amin(my_mut_score) - np.amin(my_wt_score)

        result = {'rosetta':ddG_rosetta,
                  'mean':ddG_mean,
                  'max_diff':ddG_max_diff,
                  'minimum':ddG_min,
                  'maximum':ddG_max,
                  'wt_score':wt_score,
                  'mut_score':mut_score,
                  'my_wt_score':my_wt_score,
                  'my_mut_score':my_mut_score,
                  'my_ddg':my_ddG}

        return result

class FoldX(object):
    '''
    This class implements a method to run a FoldX simulation.
    '''


    def __init__(self, pdb_name, sequence_parent = None, sequence_child = None,
                 repair_flag = True, new_name = None, clean_up = False,
                 Temp = 310, pH = 7, number_of_seeds = 1,
                 fixed_during_repair = None):
        '''
        * pdb_name refers to a structure that should be mutated and analyzed#
            without '.pdb'!
        * sequence_parent is the sequence of the pdb (pdb_name)
        * sequence_child is a mutated version of sequence_parent
        * repair_flag indicates, if the pdb should be repaired at first
        * new_name is the name of the scored structure
        * Temp is the temperature of the experiment in K
        * pH is the pH of the experiment
        * number_of_seeds is the number of runs with different starting
            conditions (in the end the mean will be returned)
        * fixed_during_repair are the residues which should be kept fixed
            Format: (aminoacid)(chain)(number), e.g. LC7, RA17
        '''
        # repair or don't
        self.repair_flag = repair_flag
        # specify residues which should stay fixed during the repair process
        # does not work ...
        self.fixed_during_repair = fixed_during_repair

        # all names miss the *.pdb at the end
        if pdb_name.endswith('.pdb'):
            self.pdb_name = pdb_name.replace('.pdb', '')
        else:
            self.pdb_name = pdb_name

        # name of a repaired pdb
        self.repaired_pdb_name = "RepairPDB_{0}".format(self.pdb_name)

        # name for build model
        if self.repair_flag is True:
            self.build_model_pdb = self.repaired_pdb_name
        else:
            self.build_model_pdb = self.pdb_name

        # name after build model
        if new_name is not None:
            if new_name.endswith('.pdb'):
                self.new_pdb_name = new_name.replace('.pdb', '')
            else:
                self.new_pdb_name = new_name
        else:
            self.new_pdb_name = '{0}_1'.format(self.build_model_pdb)

        self.sequence_parent = sequence_parent
        self.sequence_child = sequence_child

        self.clean_up = clean_up

        self.Temp = Temp
        self.pH = pH

        # number of runs to perform, in the end the mean will be returned
        self.number_of_seeds = number_of_seeds

        # filenames
        # list.txt
        self.pdb_list = 'list_{0}.txt'.format(self.pdb_name)
        # mutant_file.txt
        self.mutant_file = 'mutant_file_{0}.txt'.format(self.pdb_name)
        # run-repair.txt
        self.run_repair_file = 'run_repair_{0}.txt'.format(self.pdb_name)
        # run-buildModel
        self.run_build_file = 'run_buildModel_{0}.txt'.format(self.pdb_name)


        # FoldX path
        self.foldx_path = './FoldX_3_5_1'

        #ddG
        # stable < less stable < totally unstable
        self.ddG = None

    def _make_FoldX_run_Repair(self):
        with open(self.run_repair_file, 'w') as outdata:
            # first part
            outdata.write('<TITLE>FOLDX_runscript;\n<JOBSTART>#;\n<PDBS>#;\n<BATCH>{0};\n<COMMANDS>FOLDX_commandfile;\n'.format(self.pdb_list))
            if self.fixed_during_repair is None:
                outdata.write('<RepairPDB>#;\n')
            else:
                outdata.write('<RepairPDB>#, Fixed:{0};\n'.format(self.fixed_during_repair))
            outdata.write('<END>#;\n<OPTIONS>FOLDX_optionfile;\n<Temperature>')
            # write Temperature
            outdata.write(str(self.Temp) + ';\n')
            # write more stuff
            outdata.write('<R>#;\n<pH>')
            # write PH
            outdata.write(str(self.pH) + ';\n')
            # write the rest
            outdata.write('<IonStrength>0.050;\n<water>-CRYSTAL;\n<metal>-CRYSTAL;\n<VdWDesign>2;\n<OutPDB>true;\n<pdb_hydrogens>false;\n<END>#;\n<JOBEND>#;\n<ENDFILE>#;')


    def _make_FoldX_run_Build(self):
        with open(self.run_build_file, 'w') as outdata:
            # first part
            outdata.write('<TITLE>FOLDX_runscript;\n<JOBSTART>#;\n<PDBS>#;\n<BATCH>{0};\n<COMMANDS>FOLDX_commandfile;\n<BuildModel>#,{1};\n<END>#;\n<OPTIONS>FOLDX_optionfile;\n<Temperature>'.format(self.pdb_list, self.mutant_file))
            # write Temperature
            outdata.write(str(self.Temp) + ';\n')
            # write more stuff, write PH
            outdata.write('<R>#;\n<pH>{0};\n'.format(self.pH))
            # write number of seeds
            outdata.write('<numberOfRuns>{0};\n'.format(self.number_of_seeds))
            # write the rest
            outdata.write('<IonStrength>0.050;\n<water>-CRYSTAL;\n<metal>-CRYSTAL;\n<VdWDesign>2;\n<OutPDB>true;\n<pdb_hydrogens>false;\n<complex_with_DNA>false;\n<END>#;\n<JOBEND>#;\n<ENDFILE>#;')


    # repair!
    def _make_FoldX_mut_List(self):
        with open(self.mutant_file, 'w') as individualList:
            individualList.write(self.sequence_parent + '\n')
            individualList.write(self.sequence_child)


    def _make_FoldX_pdb_List(self, pdb_name):
        with open(self.pdb_list, 'w') as pdblist:
            pdblist.write("{0}.pdb".format(pdb_name))

    def _extract_FoldX_ddG(self):
        ddGfilename = 'Dif_BuildModel_{0}.fxout'.format(self.build_model_pdb)
        ddGcontent = csv.reader(open(ddGfilename), delimiter = '\t')
        content = []
        for line in ddGcontent:
            content.append(line)

        mean_ddG = []
        for item in range(9, 9 + self.number_of_seeds):
            mean_ddG.append(float(content[item][1]))

        self.ddG = np.mean(mean_ddG)

    def extract_FoldX_ddG(self):
        return self.ddG

    # repair
    def run_FoldX(self):
        # check if rotabase exists
        if not os.path.exists('rotabase.txt'):
            raise AttributeError('Missing rotabase.txt in {0}!'.format(os.getcwd()))
        else:
            # trash goes nowhere
            foldX_trash_file = '/dev/null'

            # repair PDB
            if self.repair_flag is True:
                self._make_FoldX_pdb_List(self.pdb_name)
                self._make_FoldX_run_Repair()
                f = open(foldX_trash_file, 'w')
                subprocess.call([self.foldx_path, '-runfile', self.run_repair_file], stdout = f)
                f.close()



            self._make_FoldX_pdb_List(self.build_model_pdb)
            self._make_FoldX_mut_List()
            self._make_FoldX_run_Build()
            f = open(foldX_trash_file, 'w')
            subprocess.call([self.foldx_path, '-runfile', self.run_build_file], stdout = f)
            f.close()
            # remove FoldX header ...
            if self.number_of_seeds == 1:
                f = open('{0}_1.pdb'.format(self.build_model_pdb), 'r')
                pdb_data = f.readlines()
                f.close()
                with open('{0}.pdb'.format(self.new_pdb_name), 'w') as f:
                    for line in pdb_data:
                        if (line[0:4] == 'ATOM' or line[0:3] == 'TER'
                            or line[0:6] == 'HETATM' or line[0:3] == 'END'):
                            f.write(line)
                # test for sequence missmatch
                new_pdb = PDBFile('{0}_1.pdb'.format(self.build_model_pdb))
                new_sequence = new_pdb.get_sequence().values()[0]
                if new_sequence != self.sequence_child:
                    raise ValueError("The sequence of the mutant is not as it should be!\n'{0}' != '{1}'".format(new_sequence, self.sequence_child))
            else:
                for i in range(self.number_of_seeds):
                    f = open('{0}_1_{1}.pdb'.format(self.build_model_pdb, i), 'r')
                    pdb_data = f.readlines()
                    f.close()
                    with open('{0}_{1}.pdb'.format(self.new_pdb_name, i), 'w') as f:
                        for line in pdb_data:
                            if (line[0:4] == 'ATOM' or line[0:3] == 'TER'
                                or line[0:6] == 'HETATM' or line[0:3] == 'END'):
                                f.write(line)

                    # test for sequence missmatch
                    new_pdb = PDBFile('{0}_{1}.pdb'.format(self.new_pdb_name, i))
                    new_sequence = new_pdb.get_sequence().values()[0]
                    if new_sequence != self.sequence_child:
                        raise ValueError("The sequence of the mutant is not as it should be!\n'{0}' != '{1}'".format(new_sequence, self.sequence_child))

            # get ddG
            self._extract_FoldX_ddG()

            # clean up
            if self.clean_up is True:
                os.remove(self.pdb_list)
                os.remove(self.mutant_file)
                os.remove(self.run_build_file)
                if self.repair_flag is True:
                    # names
                    repaired_name = self.repaired_pdb_name
                    # remove
                    os.remove(self.run_repair_file)
                    os.remove('RepairPDB_{0}.fxout'.format(self.pdb_name))
                    os.remove('{0}.pdb'.format(repaired_name))

                os.remove('BuildModel_{0}.fxout'.format(self.build_model_pdb))
                os.remove('PdbList_BuildModel_{0}.fxout'.format(self.build_model_pdb))
                os.remove('Dif_BuildModel_{0}.fxout'.format(self.build_model_pdb))
                os.remove('Average_BuildModel_{0}.fxout'.format(self.build_model_pdb))
                os.remove('Raw_BuildModel_{0}.fxout'.format(self.build_model_pdb))
                if self.number_of_seeds == 1:
                    os.remove('WT_{0}_1.pdb'.format(self.build_model_pdb))
                    os.remove('{0}_1.pdb'.format(self.build_model_pdb))
                else:
                    for i in range(self.number_of_seeds):
                        os.remove('WT_{0}_1_{1}.pdb'.format(self.build_model_pdb, i))
                        os.remove('{0}_1_{1}.pdb'.format(self.build_model_pdb, i))


class Eris(object):

    def __init__(self, unique_id, storage_path, reference_pdb_path, new_seq,
                 pre_relaxation = True, flexibile_backbone = True,
                 eris_path = None, remove_unfitted_pdbs = True):
        '''
        This calls eris to perform the mutation and the result can be obtained
        by 'calculate_ddG'.



        '''
        # unique id
        self.unique_id = unique_id

        # storage path
        self.storage_path = storage_path

        self.ref_pdb_path = reference_pdb_path
        self.ref_pdb = PDBFile(self.ref_pdb_path)
        self.chain = self.ref_pdb.structure.child_list[0].id
        self.old_seq = self.ref_pdb.get_sequence()[self.chain]
        self.new_seq = new_seq

        # check sequence length
        if len(self.old_seq) != len(new_seq):
            raise AttributeError('Sequencelength missmatch! Reference sequence {0} vs {1} new sequence.'.format(len(self.old_seq), len(new_seq)))

        # get mutation instructions
        self.mut_instructions = self._get_mut_instructions(self.old_seq, new_seq)

        # get path of eris
        if eris_path is not None:
            self.eris_path = eris_path
        else:
            self.eris_path = os.path.join(os.getenv('HOME'), 'Programs',
                                          'eris_standalone', 'bin', 'eris.pl')

        # pre relaxation
        if pre_relaxation is True:
            self.pre_relaxation = '1'
        else:
            self.pre_relaxation = '0'

        # flexible backbone
        if flexibile_backbone is True:
            self.flexible_backbone = '1'
        else:
            self.flexible_backbone = '0'

        # mutant path
        self.mutant_pdb_path = ['Dir_{0}'.format(self.storage_path)]
        self.mutant_pdb_path.append('FLEX')
        self.mutant_pdb_path.append('pdb-{0}'.format(self.unique_id))
        self.mutant_pdb_path = '{0}'.format(os.path.sep).join(self.mutant_pdb_path)


        self.remove_unfitted_pdbs = remove_unfitted_pdbs

    def calculate_ddG(self):
        '''
        This method runs eris and returns the score.
        '''
        # check if this unique id has already been calculated
        # this could have happend when the GA crashed
        result_path = ['Dir_{0}'.format(self.storage_path)]
        result_path.append('DATA')
        result_path.append('mutations')
        result_path = '{0}'.format(os.path.sep).join(result_path)
        if os.path.exists(result_path):
            with open(result_path) as f:
                content = f.readlines()

            for line in content:
                if line.split(' ')[0] == self.unique_id:
                    raise AttributeError("This unique_id '{0}' has already been calculated! {1}".format(self.unique_id, self.new_seq))

        cmd = [self.eris_path]
        cmd.append(self.unique_id)
        cmd.append(self.storage_path)
        cmd.append(self.ref_pdb_path)
        cmd.append(self.pre_relaxation)
        cmd.append(self.flexible_backbone)
        cmd.append(self.mut_instructions)

        trash_file = '/dev/null'
        f = open(trash_file, 'w')
        subprocess.call(cmd, stdout = f, stderr = f)
        f.close()

        pdb_list = os.listdir(self.mutant_pdb_path)
        for i, pdb_name in enumerate(pdb_list):
            # what a shitty program ... it can happen, that a pdb file exists
            # with a lot of numbers but there are no coordinates in it? wtf???
            # so here we go and read the file and check if there are
            # coordinates
            pdb_with_atom_coord = False
            with open(os.path.join(self.mutant_pdb_path, pdb_name)) as f:
                for line in f:
                    if line.startswith('ATOM'):
                        pdb_with_atom_coord = True
                        break

            if pdb_with_atom_coord is False:
                os.remove(os.path.join(self.mutant_pdb_path, pdb_name))

        pdb_list = os.listdir(self.mutant_pdb_path)

        if len(pdb_list) == 0:
            # no valid structures
            return float('NaN')

        # check for wrong sequences
        for i, pdb_name in enumerate(pdb_list):
            mutant_pdb = PDBFile(os.path.join(self.mutant_pdb_path, pdb_name))
            if mutant_pdb.get_sequence().values()[0] != self.new_seq:
                # eris is such a nice programm, it can happen, that a
                # mutant has not the new sequence ...
                #os.remove(mutant_pdb.structure_path)
                os.remove(os.path.join(self.mutant_pdb_path, pdb_name))

        pdb_list = os.listdir(self.mutant_pdb_path)

        if len(pdb_list) == 0:
            # no valid structures
            return float('NaN')

        for i, pdb_name in enumerate(pdb_list):
            ref_pdb = PDBFile(self.ref_pdb_path)
            mutant_pdb = PDBFile(os.path.join(self.mutant_pdb_path, pdb_name))

            ref_pdb.superimpose_given_pdb_onto_self(mutant_pdb, 'CA')

            # give them a new name
            mutant_pdb.save_to_file('{0}_{1}.pdb'.format(os.path.join(self.mutant_pdb_path, self.unique_id), i))

        # read data
        with open(result_path) as f:
            new_content = f.readlines()
        new_content.reverse()
        for line in new_content:
            if line.split(' ')[0] == self.unique_id:
                eris_ddg = float(line.strip().split(' ')[-1])

        return eris_ddg

    def _get_mut_instructions(self, old_seq, new_seq):
        '''
        Get mutation instructions:
            old_seq: AYTP
            index:   1234
            new_seq: ARTN

        return '\"Y2R P4N\"'
        '''
        # get first index
        index = self.ref_pdb.structure[self.chain].child_list[0].id[1]

        if old_seq == new_seq:
            # dummy calculation
            print('dummy calculation')
            return '\"{0}{1}{2}\"'.format(old_seq[0], index, old_seq[0])

        mut_instructions = []
        # chains are of the same length, already checked that
        for i in range(len(old_seq)):
            if old_seq[i] != new_seq[i]:
                mut_instructions.append('{0}{1}{2}'.format(old_seq[i],
                                                           index + i,
                                                           new_seq[i]))
        mut_instructions = ' '.join(mut_instructions)
        return '\"{0}\"'.format(mut_instructions)
