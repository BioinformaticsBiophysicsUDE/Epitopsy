'''
Created on Jan 10,  2012

@author: chris
'''

import os
import stat
import subprocess
import numpy as np

import epitopsy
from epitopsy.Structure import PDBFile, PQRFile
from epitopsy.DXFile import DXReader

class APBSWrapper:
    '''
    classdocs
    '''


    def __init__(self, APBSPath = "apbs", PDB2PQRPath = "pdb2pqr"):
        self.APBSPath = APBSPath
        self.PDB2PQRPath = PDB2PQRPath
        self.binding_energy = None
        self.debye_length = []

    def runPDB2PQR(self, pdb_path, pqr_path, force_field = 'amber', pdb2pqr_argv = None):
        '''Call pdb2pqr to replace the b-factor and the occupancy information
        in the pdb file with the charges and the vdw radii.

        Args:
            pdb_path -> Path to the pdb file.
            pqr_path -> Path for the new pqr file.
            ph -> ph at for which the charges should be calculated. Default is
                None, which is close to a ph of 7.
            force_field -> Forcefield from which charges and radii should be
                taken. Default is amber.
            pdb2pqr_argv -> Can contain additional arguments to pdb2pqr as a
                list (e.g. ['--assign-only'], oder ['--noopt']). If multiple
                additional arguments are given, they also have to be given as
                a list (e.g. ['--assign-only', '--noopt']).

        Returns:
            None.
        '''
        pdb_object = PDBFile(pdb_path)
        pqr_object = pdb_object.get_pqr_structure(pqr_path, force_field = force_field,
                pdb2pqr_argv=pdb2pqr_argv)


    def runAPBS(self, apbsInParameters, inFilename):
        '''
        Args:
            apbsInParameters -> Contains all neccessary parameters.
            inFilename -> Name of the inputfile for apbs.

        Returns:
            None.
        '''
        if not(os.path.exists(apbsInParameters.pqr_path)):
            raise NameError('Error running APBS: File not found: {0}'.format(apbsInParameters.getPQRFilePath()))
        apbsInParameters.write(inFilename)

        p = subprocess.Popen([self.APBSPath, inFilename], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        apbs_output, apbs_error = p.communicate()

        apbs_error = apbs_error.split('\n')
        apbs_output = apbs_output.split('\n')
        # there are some error messages (something with fscan) but they
        # are not relevant as they just say that the programm encountered
        # an end of file during the calculation
        for line in apbs_error:
            if line.startswith('Error while parsing input file.'):
                raise AttributeError('Error while parsing input file to apbs.')

        for line in apbs_output:
            # find debye length
            if line.startswith('  Debye length:'):
                self.debye_length.append(float(line.split(':')[1].rstrip('A')))
            # get binding energy
            if(apbsInParameters.calculation_type == apbsInParameters.binding_energy
               or apbsInParameters.calculation_type == apbsInParameters.binding_energy_long):
                if line.startswith('  Global net ELEC energy ='):
                    line = line.rstrip(' kJ/mol')
                    self.binding_energy = float(line.split('=')[1])


        #os.remove('io.mc')

    def get_dxbox(self, pdb_path, mesh_size, no_return = False, **kwds):
        '''
        This method starts pdb2pqr and apbs and returns the specified boxes.
        Be careful with absolute paths!

        Args:
            pdb_path -> path to the structure
            mesh_size -> grid dimensions
            --- below are optional
            pqr_path -> if None is supplied, it uses pdb2pqr to calculate one
            box_dim -> if box dimensions and a box center are supplied, it
                will not use the pdb for the construction of the
                box. Furthermore it is not possible to extend a
                box, if the box properties are supplied.
            box_center -> determine the center of the box. If None is
                supplied, it calculates the center of the pqr file.
            extend -> Extend grid by the given number of Angstroems
            box_type -> specify the returned objects: 'esp', 'vdw', 'smol'
            cubic_box -> if the box is generated from the pqr file, it can
                be made to be cubic
            close_boundaries -> if True it uses another algorithm for the
                boundary conditions of the box which
                yields better results but is slower
            temperatur -> well it is the Temperatur in Kelvin
            ph -> it's the pH

        Returns:
            A dictionary with DXBox objects, which keys are the requested
            box types.
        '''
        if 'box_dim' in kwds:
            box_dim = kwds['box_dim']
        else:
            box_dim = None
        if 'box_center' in kwds:
            box_center = kwds['box_center']
        else:
            box_center = None
        if 'extend' in kwds:
            extend = kwds['extend']
        else:
            extend = None
        if 'box_type' in kwds:
            box_type = kwds['box_type']
        else:
            box_type = ['esp', 'vdw']
        if 'cubic_box' in kwds:
            cubic_box = kwds['cubic_box']
        else:
            cubic_box = True
        if 'close_boundaries' in kwds:
            close_boundaries = kwds['close_boundaries']
        else:
            close_boundaries = True
        if 'temperature' in kwds:
            temperature = kwds['temperature']
        else:
            temperature = 310
        if 'ph' in kwds:
            ph = kwds['ph']
        else:
            ph = None

        pdb = PDBFile(pdb_path)
        if 'pqr_path' in kwds:
            pqr_path = kwds['pqr_path']

            if pqr_path is None:
                # this could somehow happen
                # if the given path is None, calculate it
                pqr_path = pdb_path.replace('.pdb', '.pqr')
                pdb2pqr_argv = []

                if 'pdb2pqr_argv' in kwds:
                    pdb2pqr_argv = kwds['pdb2pqr_argv']
                if ph is None:
                    pqr = pdb.get_pqr_structure(pqr_path, pdb2pqr_argv=pdb2pqr_argv)
                    pqr.save_to_file(pqr_path)
                else:
                    pdb2pqr_argv.append("--with-ph={0}".format(ph))
                    pqr = pdb.get_pqr_structure(pqr_path, pdb2pqr_argv=pdb2pqr_argv)
                    pqr.save_to_file(pqr_path)

            else:
                # the pqr path should exist
                if not os.path.exists(pqr_path):
                    raise ValueError("Could not find pqr path: {0}".format(pqr_path))

        else:
            pqr_path = pdb_path.replace('.pdb', '.pqr')
            pdb2pqr_argv = []
            if 'pdb2pqr_argv' in kwds:
                pdb2pqr_argv = kwds['pdb2pqr_argv']
            if ph is None:
                pqr = pdb.get_pqr_structure(pqr_path, pdb2pqr_argv=pdb2pqr_argv)
                pqr.save_to_file(pqr_path)
            else:
                pdb2pqr_argv.append("--with-ph={0}".format(ph))
                pqr = pdb.get_pqr_structure(pqr_path, pdb2pqr_argv=pdb2pqr_argv)
                pqr.save_to_file(pqr_path)

        template_in = InFile(pqr_path = pqr_path,
                             calculation_type = 'potential',
                             box_mesh_size = mesh_size,
                             box_dim = box_dim,
                             box_center = box_center,
                             extend = extend, cubic_box = cubic_box,
                             box_type = box_type)

        template_in.temp = temperature

        if close_boundaries:
            template_in.bcfl = "mdh"
        else:
            template_in.bcfl = 'sdh'

        self.runAPBS(template_in, '{0}.in'.format(pqr_path.replace('.pqr', '')))
        '''
        Read the grids as boxes, no need to use the grid class, because i
        need arrays:
        '''
        box_dict = {}

        for item in template_in.box_type:
            def fix_apbs_filename(orig_filename):
                '''
                APBS installed from the ubuntu repositories inserts "-PE0"
                into the generated dx files. This function checks which 
                output was generated
                '''
                filename = orig_filename.replace("-PE0.dx",".dx")
                variant_1 = filename
                variant_2 = filename.replace(".dx","-PE0.dx")
                if (os.path.exists(variant_1) and os.path.exists(variant_2)):
                    raise ValueError("Both files exist, '{0}' and '{1}'!".format(
                        variant_1,
                        variant_2))
                elif os.path.exists(variant_1):
                    return variant_1
                elif os.path.exists(variant_2):
                    return variant_2
                else:
                    raise ValueError("Could not find file: {0}".format(orig_filename))
                    
            filename = "{0}_{1}.dx".format(pqr_path.replace('.pqr', ''), item)
            if not no_return:
                box_dict[item] = DXReader().parse(fix_apbs_filename(filename), item)

        return box_dict

    def get_binding_energy(self, complex_pqr_path, ligand_pqr_path,
                           fixed_pqr_path, box_mesh_size, extend = None,
                           ** kwds):
        '''
        Binding energy is in kJ/mol.
        Obtain the change in total polar solvation energy by:
            dG_bind = comp_solv - fixed_solv - ligand_solv
        It gives the total polar solvation energy which includes both the
        reaction field and coulombic contributions.
        Positive values are favorable.

        **kwds is used to set the options of the infile for the apbs
        calculation.

        Args:
            complex_pqr_path -> Path of the complex of ligand and fixed
                protein.
            ligand_pqr_path -> Path of the isolated ligand.
            fixed_pqr_path -> Path of the isolated fixed protein.
            box_mesh_size -> Specify the mesh size of the grid in Angstroem.
            extend -> Increase the box dimensions.

        Returns:
            Binding energy as a float.
        '''
        inFilename = complex_pqr_path.replace('.pqr', '.in')
        pqr = PQRFile(complex_pqr_path)
        box_center = pqr.determine_geometric_center()
        box_dim = pqr.get_dxbox_dim(box_mesh_size, extend=extend)
        apbsInParameters = InFile(pqr_path = complex_pqr_path,
                                  calculation_type = 'binding_energy',
                                  box_mesh_size = box_mesh_size,
                                  box_dim=box_dim,
                                  box_center=box_center,
                                  ligand_pqr_path = ligand_pqr_path,
                                  fixed_pqr_path = fixed_pqr_path)

        apbsInParameters.set_options(kwds)

        self.runAPBS(apbsInParameters, inFilename)
        return self.binding_energy

    def get_binding_energy_long(self, complex_pqr_path, ligand_pqr_path,
                           fixed_pqr_path, box_mesh_size, extend = None,
                           ** kwds):
        '''
        Binding energy is in kJ/mol.

        solvation_energy = (complex_solv - complex_ref)
                            - (ligand_solv - ligand_ref)
                            - (fixed_solv - fixed_ref)

        coulomb_energy = (complex_coulomb - ligand_coulomb - fixed_coulomb) / pdie

        binding_energy = solvation_energy + coulomb_energy

        Positive values are favorable.

        **kwds is used to set the options of the infile for the apbs
        calculation.

        Args:
            complex_pqr_path -> Path of the complex of ligand and fixed
                protein.
            ligand_pqr_path -> Path of the isolated ligand.
            fixed_pqr_path -> Path of the isolated fixed protein.
            box_mesh_size -> Specify the mesh size of the grid in Angstroem.
            extend -> Increase the box dimensions.

        Returns:
            Binding energy as a float.
        '''
        inFilename = complex_pqr_path.replace('.pqr', '.in')
        pqr = PQRFile(complex_pqr_path)
        box_center = pqr.determine_geometric_center()
        box_dim = pqr.get_dxbox_dim(box_mesh_size, extend=extend)
        apbsInParameters = InFile(pqr_path = complex_pqr_path,
                                  calculation_type = 'binding_energy_long',
                                  box_mesh_size = box_mesh_size,
                                  box_dim=box_dim,
                                  box_center=box_center,
                                  ligand_pqr_path = ligand_pqr_path,
                                  fixed_pqr_path = fixed_pqr_path)

        apbsInParameters.set_options(kwds)

        self.runAPBS(apbsInParameters, inFilename)

        # call coulomb script
        coulomb_energy_complex = get_coulomb_energy(complex_pqr_path)
        coulomb_energy_ligand = get_coulomb_energy(ligand_pqr_path)
        coulomb_energy_fixed = get_coulomb_energy(fixed_pqr_path)
        coulomb_energy = (coulomb_energy_complex
                          - coulomb_energy_ligand
                          - coulomb_energy_fixed)

        # correct for dieelectric constant
        coulomb_energy = coulomb_energy / float(apbsInParameters.pdie)

        # get the solvation energy
        solvation_energy = self.binding_energy
        total_energy = solvation_energy + coulomb_energy
#        print("solv: {0}, coul: {1} -> bind: {2}".format(solvation_energy, coulomb_energy, total_energy))

        return total_energy

    def get_dissociation_energy(self, complex_pqr_path, ligand_pqr_path,
                                fixed_pqr_path, box_mesh_size, extend = None,
                                **kwds):
        '''
        Dissociation energy is in kJ/mol.
        Obtain the change in total polar solvation energy by:
            dG_diss = - (comp_solv - fixed_solv - ligand_solv)
        It gives the total polar solvation energy which includes both the
        reaction field and coulombic contributions.

        **kwds is used to set the options of the infile for the apbs
        calculation.

        Args:
            complex_pqr_path -> Path of the complex of ligand and fixed
                protein.
            ligand_pqr_path -> Path of the isolated ligand.
            fixed_pqr_path -> Path of the isolated fixed protein.
            box_mesh_size -> Specify the mesh size of the grid in Angstroem.
            extend -> Increase the box dimensions.

        Returns:
            Dissociation energy as a float.
        '''
        inFilename = complex_pqr_path.replace('.pqr', '.in')
        apbsInParameters = InFile(pqr_path = complex_pqr_path,
                                  calculation_type = 'binding_energy',
                                  box_mesh_size = box_mesh_size,
                                  extend = extend,
                                  cubic_box = True,
                                  ligand_pqr_path = ligand_pqr_path,
                                  fixed_pqr_path = fixed_pqr_path)

        apbsInParameters.set_options(kwds)

        self.runAPBS(apbsInParameters, inFilename)
        return -self.binding_energy

def get_coulomb_energy(protein_pqr_path):
    '''
    Args:
        protein_pqr_path -> Path to the pqr of the protein.

    Return:
        Energy in kJ/mol.
    '''
    coulomb_path = os.path.join(epitopsy.__path__[0], "external_scripts", "coulomb")
    # make sure everyone can execute the file
    os.chmod(coulomb_path, stat.S_IRWXU | stat.S_IRGRP | stat.S_IXGRP | stat.S_IROTH | stat.S_IXOTH)
    p = subprocess.Popen([coulomb_path, protein_pqr_path], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    coulomb_output, coulomb_error = p.communicate()

    coulomb_error = coulomb_error.split('\n')
    coulomb_output = coulomb_output.split('\n')

    energy = None
    for line in coulomb_output:
        if line.startswith("Total energy ="):
            content = line.replace("kJ/mol in vacuum.","")
            content = content.split("=")
            energy = float(content[-1])

    return energy


class InFile:
    '''
    This class holds parameters to run APBS calculations.
    Only the most basic parameters are covered here,  intended for the
    generation of basic grids holding information of electrostatic
    potential and/or the van-der-Waals surface
    '''
    OUTPUT_TYPE_VDW = 'vdw'
    OUTPUT_TYPE_ESP = 'esp'
    OUTPUT_TYPE_SMOL = 'smol'
    OUTPUT_TYPE_NDENS = 'ndens'
    OUTPUT_TYPE_QDENS = 'qdens'

    binding_energy = 'binding_energy'
    binding_energy_long = 'binding_energy_long'
    potential = 'potential'
    def __init__(self, pqr_path, calculation_type, box_mesh_size,
                 box_dim = None, box_center = None,
                 extend = None, cubic_box = True,
                 box_type = None, ligand_pqr_path = None,
                 fixed_pqr_path = None):
        '''
        Args:
            pqr_path -> pqr path, from which the potential energy should be
                calculated or the path of the complex for which the binding
                energy should be calculated
            calculation_type -> potential or binding_energy
            box_mesh_size -> mesh size
            box_dim -> special dimension can be supplied (although they might be
                fixed to work properly with apbs), if None is given it
                will be calculated from the size of the protein.
            box_center -> center of the box, if None is given if will be set to the
                geometric center of the protein
            extend -> can only be used when no box_dim is supplied, extends the box
                size by the given amount (in angstroem)
            cubic_box -> determine wheter it is a cubic box or not, in the case of
                a cubic box the largest dimension is used for all
                dimensions
            box_type ->
                esp : electrostatic potential
                vdw : van der waals based solvent accessibility
                smol : solvent accessibility
                ndens : total mobile ion number density in units of M
                qdens : total mobile charge density in units of e_c M
            ligand_pqr_path -> needs to be supplied if the calculation is a
                binding energy calculation
            fixed_pqr_path -> needs to be supplied if the calculation is a
                binding energy calculation

        Returns:
            None.
        '''
        # type of calculation
        if calculation_type == self.potential:
            self.calculation_type = self.potential
        elif calculation_type == self.binding_energy:
            self.calculation_type = self.binding_energy
            if ligand_pqr_path is None or fixed_pqr_path is None:
                raise ValueError('For a binding energy calculation a ligand and a fixed protein are necessary!')
        elif calculation_type == self.binding_energy_long:
            self.calculation_type = self.binding_energy_long
            if ligand_pqr_path is None or fixed_pqr_path is None:
                raise ValueError('For a binding energy calculation a ligand and a fixed protein are necessary!')
        else:
            raise ValueError("Unknown Argument for 'calculation_type': '{0}'!".format(calculation_type))

        # grid parameters
        self.gridSizeX = 0
        self.gridSizeY = 0
        self.gridSizeZ = 0
        self.meshSizeX = 0.0
        self.meshSizeY = 0.0
        self.meshSizeZ = 0.0
        # grid center
        self.centerX = 0.0
        self.centerY = 0.0
        self.centerZ = 0.0
        # file informations
        self.pqr_path = pqr_path
        self.ligand_pqr_path = ligand_pqr_path
        self.fixed_pqr_path = fixed_pqr_path

        # check for list
        if isinstance(box_type, str):
            box_type = [box_type]

        self.box_type = box_type

        # various parameters
        # elec statement
        #     mg-auto -> Automatically-configured sequential focusing
        #                multigrid Poisson-Boltzmann calculations.
        #     mg-para -> Automatically-configured parallel focusing
        #                multigrid Poisson-Boltzmann calculations.
        #     mg-manual -> Manually-configured multigrid Poisson-Boltzmann
        #                  calculations.
        #     fe-manual -> Manually-configured adaptive finite element
        #                  Poisson-Boltzmann calculations.
        #     mg-dummy -> Calculations of surface and charge distribution
        #                 properties which do not require solution of the PBE.
        #self.elec_type = 'mg-manual'
        self.elec_type = 'mg-auto'
        # type of Poisson Boltzmann function:
        #     npbe
        #     lpbe
        self.pbeType = 'npbe'
        # type of boundary conditions
        #     sdh -> fast, works best, if boundaries are far away
        #     mdh -> slow, works better then sdh, if boundaries are close
        self.bcfl = "mdh"
        # protein interior dielectric constant
        self.pdie = 2.0
        # solvent dielectric constant
        self.sdie = 79.0
        # depth of multilevel hierarchy? marked as deprecated
        self.nlev = 4
        # model to construct the dielectric and ion-accessibiity
        #     mol -> dielectric coefficient is defined based on molecular
        #            surface definition
        #     smol -> same as mol, but smoothed with a 9-point harmonic
        #            averaging to reduce sensitivity to grid setup
        #     spl2 -> dielectric and ion-accessibility are defined by a cubic
        #            spline surface, can generate unphysical results
        #     spl4 -> dielectric and ion-accessibility are defined a 7th order
        #            polynomial
        self.srfm = "smol"
        # method by which point charges (i.e. Dirac delta functions) are mapped
        # to the grid
        #    sp10 -> trilinear interpolation (linera splines), very sensitive
        #            to grid setup
        #    sp12 -> cubic B-spline interpolation, not so sensitive
        #    sp14 -> quintic B-spline interpolation (125 grid points receive
        #            charge density, next-next-nearest neighbors)
        self.chgm = "spl2"
        # number of grid points per square angstrom to use in discontinuos
        # surface constructions
        self.sdens = 10.0
        # radius of the solvent molecules, water = 1.4
        self.srad = 1.4
        # size of the support (i.e., rate of change) for spline based surface
        # definitions
        self.swin = 0.3
        # Temperature
        self.temp = 310.00
        # ion charge -> bulk system must be electroneutral!
        #     charge in e_c
        #     conc in M
        #     radius in A
        self.ion_pos_charge = 1
        self.ion_pos_conc = 0.15
        self.ion_pos_radius = 2.0 # from APBS Tutorial
        self.ion_neg_charge = -self.ion_pos_charge
        self.ion_neg_conc = self.ion_pos_conc
        self.ion_neg_radius = 2.0 # from APBS Tutorial
        # tolerance for iterations of the pmg partial differential equation
        # solver
        self.etol = 1.0e-6

        if(self.calculation_type == self.binding_energy
           or self.calculation_type == self.binding_energy_long):
            self.calcenergy = "total"
        else:
            self.calcenergy = 'no'
        self.calcforce = "no"

        if box_center is None:
            # if no dimension is available, determine ma
            self.setGridCenter(PQRFile(self.pqr_path).determine_geometric_center())
        else:
            self.setGridCenter(box_center)

        self.setMeshSize(box_mesh_size)

        if box_dim is not None and extend is not None:
            raise AttributeError('Input error! It does not make sense to provide the box dimension AND an extend parameter!')

        if box_dim is None:
            if cubic_box is True:
                max_dim = PQRFile(self.pqr_path).determine_max_diameter()
                box_dim = np.array([max_dim, max_dim, max_dim])
            else:
                extremes = PQRFile(self.pqr_path).determine_coordinate_extremes()
                box_dim = np.zeros(3)
                box_dim[0] = abs((extremes[0][1] - extremes[0][0]))
                box_dim[1] = abs((extremes[1][1] - extremes[1][0]))
                box_dim[2] = abs((extremes[2][1] - extremes[2][0]))

            if extend is not None:
                # mesh_size: x = y = z
                box_dim = np.ceil(box_dim / np.array(box_mesh_size)) + 2 * extend / box_mesh_size[0]
            else:
                box_dim = np.ceil(box_dim / np.array(box_mesh_size))

        self.setGridSize(box_dim)
        self.fixGridSize()

    def  setGridCenter(self, center):
        self.centerX = center[0]
        self.centerY = center[1]
        self.centerZ = center[2]


    def  setGridSize(self, size):

        self.gridSizeX = size[0]
        self.gridSizeY = size[1]
        self.gridSizeZ = size[2]


    def  setMeshSize(self, meshSize):
        self.meshSizeX = meshSize[0]
        self.meshSizeY = meshSize[1]
        self.meshSizeZ = meshSize[2]

    def  generateFromPDB(self, pdb, padding, cubicBox):
        self.generateFromPDB2(pdb, padding, 0.0, cubicBox)

    def  generateFromPDB2(self, pdb, padding, minDiameter, cubicBox):
        if(self.meshSizeX == 0.0 or self.meshSizeY == 0.0 or self.meshSizeX == 0.0):
            raise NameError("ERROR: Cannot calculate grid size mesh sizes < =  0.0!")
        ''''
         * TODO: use method to determine geometric center of protein-associated atoms only!
         * Same goes for determination of coordinate extremes. Only use protein-associated atoms!
         ???
         '''
        self.setGridCenter(pdb.determine_geometric_center())

        diameterX = 0.0
        diameterY = 0.0
        diameterZ = 0.0
        if cubicBox is True:
            maxDim = pdb.determine_max_diameter()
            diameterX = maxDim
            diameterY = maxDim
            diameterZ = maxDim
        else:
            extremes = pdb.determine_coordinate_extremes()
            diameterX = abs((extremes[0][1] - extremes[0][0]))
            diameterY = abs((extremes[1][1] - extremes[1][0]))
            diameterZ = abs((extremes[2][1] - extremes[2][0]))
        if(diameterX < minDiameter):
            diameterX = minDiameter
        if(diameterY < minDiameter):
            diameterY = minDiameter
        if(diameterZ < minDiameter):
            diameterZ = minDiameter
        self.gridSizeX = (np.ceil(diameterX / self.getMeshSizeX()) + padding * 2)
        self.gridSizeY = (np.ceil(diameterY / self.getMeshSizeY()) + padding * 2)
        self.gridSizeZ = (np.ceil(diameterZ / self.getMeshSizeZ()) + padding * 2)
        self.fixGridSize()

    def  fixGridSize(self):
        fixed_size = fix_grid_size([self.gridSizeX,
                self.gridSizeY, self.gridSizeZ])
        self.gridSizeX = fixed_size[0]
        self.gridSizeY = fixed_size[1]
        self.gridSizeZ = fixed_size[2]

    def write(self, file_path):
        if self.calculation_type == self.potential:
            self.write_potential(file_path)
        elif self.calculation_type == self.binding_energy:
            self.write_binding_energy(file_path)
        elif self.calculation_type == self.binding_energy_long:
            self.write_binding_energy_long(file_path)


    def write_potential(self, file_path):
        '''
        This is the function that writes an infile for the calculation of
        potential grids.

        Args:
            file_path -> Path to the new infile.

        Returns:
            None.
        '''
        # check if everything is okay:
        self.fixGridSize()
        if self.checkInFile() == False:
            raise NameError('Error in consistency check')
        instructions = ''

        #### read statement
        instructions = instructions + 'read\n'
        # complex
        instructions = instructions + '\tmol pqr {0}\n'.format(self.pqr_path)
        # ligand
        #instructions = instructions + '\tmol pqr {0}\n'.format(self.ligand_pqr_path)
        # fixed
        #instructions = instructions + '\tmol pqr {0}\n'.format(self.fixed_pqr_path)
        instructions = instructions + 'end\n'

        instructions = instructions + 'elec\n'
        instructions = instructions + '\t{0}\n'.format(self.elec_type)
        # grid size, needs to be: c*2^(l+1) + 1
        instructions = instructions + '\tdime {0} {1} {2}\n'.format(self.gridSizeX,
                                                                    self.gridSizeY,
                                                                    self.gridSizeZ)
        # starting grid length
        instructions = instructions + "\tcglen {0:.2f} {1:.2f} {2:.2f}\n".format(self.gridSizeX * self.meshSizeX,
                                                        self.gridSizeY * self.meshSizeY,
                                                        self.gridSizeZ * self.meshSizeZ)
        # fine mesh domain length in a multigrid focusing calculation
        shrink_fac = 1.
        instructions = instructions + "\tfglen {0:.2f} {1:.2f} {2:.2f}\n".format(
                self.gridSizeX * self.meshSizeX * shrink_fac,
                self.gridSizeY * self.meshSizeY * shrink_fac,
                self.gridSizeZ * self.meshSizeZ * shrink_fac)
        # center grid at position given by the center of the pdb
        instructions = instructions + '\tcgcent {x} {y} {z}\n'.format(x=self.centerX, y=self.centerY, z=self.centerZ)
        # center of the fine grid (in a focusing calculation)
        instructions = instructions + '\tfgcent {x} {y} {z}\n'.format(x=self.centerX, y=self.centerY, z=self.centerZ)
        # molecule for which the PBE is to be solved, definition according
        # to the list of molecules given after 'read' starting with 1
        instructions = instructions + '\tmol 1\n'
        # type of equation to be solved, lpbe, npbe, smpbe
        instructions = instructions + '\t{0}\n'.format(self.pbeType)
        instructions = instructions + '\tbcfl {0}\n'.format(self.bcfl)
        # dielectric constant of the protein
        instructions = instructions + '\tpdie {0}\n'.format(self.pdie)
        # dielectric constant of the solvent
        instructions = instructions + '\tsdie {0}\n'.format(self.sdie)
        # bulk concentration of mobile ions, total bulk system of ions has
        # to be electroneutral!
        instructions = instructions + '\tion charge {0} conc {1:.3f} radius {2:.1f}\n'.format(self.ion_pos_charge,
                                                                                              self.ion_pos_conc,
                                                                                              self.ion_pos_radius)
        instructions = instructions + '\tion charge {0} conc {1:.3f} radius {2:.1f}\n'.format(self.ion_neg_charge,
                                                                                              self.ion_neg_conc,
                                                                                              self.ion_neg_radius)
        # model to construct dielectric constant and ion accessibility
        instructions = instructions + '\tsrfm {0}\n'.format(self.srfm)
        # method by which point charges are mapped to the grid
        instructions = instructions + '\tchgm {0}\n'.format(self.chgm)
        # number of grid points per square angstrom, ignored when srad is
        # 0.0 or srfm is sp12
        instructions = instructions + '\tsdens {0}\n'.format(self.sdens)
        # radius of solvent molceules, 1.4 for water, ignored for sp12
        instructions = instructions + '\tsrad {0}\n'.format(self.srad)
        # size of support (i.e. the rate of change) for spline based
        # surface definitions
        instructions = instructions + '\tswin {0}\n'.format(self.swin)
        # Temperature
        instructions = instructions + '\ttemp {0}\n'.format(self.temp)
        # VDW
        if self.OUTPUT_TYPE_VDW in self.box_type:
            instructions = instructions + "\twrite vdw dx {0}_vdw\n".format(self.pqr_path.replace('.pqr', ''))
        # ESP
        if self.OUTPUT_TYPE_ESP in self.box_type:
            instructions = instructions + "\twrite pot dx {0}_esp\n".format(self.pqr_path.replace('.pqr', ''))
        # SMOL
        if self.OUTPUT_TYPE_SMOL in self.box_type:
            instructions = instructions + "\twrite smol dx {0}_smol\n".format(self.pqr_path.replace('.pqr', ''))
        # QDENS
        if self.OUTPUT_TYPE_QDENS in self.box_type:
            instructions = instructions + "\twrite qdens dx {0}_qdens\n".format(self.pqr_path.replace('.pqr', ''))
        # NDENS
        if self.OUTPUT_TYPE_NDENS in self.box_type:
            instructions = instructions + "\twrite ndens dx {0}_ndens\n".format(self.pqr_path.replace('.pqr', ''))
        instructions = instructions + 'end\n'
        instructions = instructions + 'quit\n'

        with open(file_path, 'w') as f:
            f.write(instructions)

    def write_binding_energy(self, file_path):
        '''
        This is the function that writes an infile for binding energy
        calculations.

        Args:
            file_path -> Path to the new infile.

        Returns:
            None.
        '''
        # check if everything is okay:
        self.fixGridSize()
        if self.checkInFile() == False:
            raise NameError('Error in consistency check')
        instructions = ''

        #### read statement
        instructions = instructions + 'read\n'
        # complex
        instructions = instructions + '\tmol pqr {0}\n'.format(self.pqr_path)
        # ligand
        instructions = instructions + '\tmol pqr {0}\n'.format(self.ligand_pqr_path)
        # fixed
        instructions = instructions + '\tmol pqr {0}\n'.format(self.fixed_pqr_path)
        instructions = instructions + 'end\n'

        for i in range(1, 4):
            if i == 1:
                elec_name = 'complex'
            elif i == 2:
                elec_name = 'ligand'
            elif i == 3:
                elec_name = 'fixed'
            #### elec complex
            instructions = instructions + 'elec name {0}\n'.format(elec_name)
            instructions = instructions + '\t{0}\n'.format(self.elec_type)
            # grid size, needs to be: c*2^(l+1) + 1
            instructions = instructions + '\tdime {0} {1} {2}\n'.format(self.gridSizeX,
                                                                        self.gridSizeY,
                                                                        self.gridSizeZ)
            # starting grid length
            instructions = instructions + "\tcglen {0:.2f} {1:.2f} {2:.2f}\n".format(self.gridSizeX * self.meshSizeX,
                                                            self.gridSizeY * self.meshSizeY,
                                                            self.gridSizeZ * self.meshSizeZ)
            # fine mesh domain length in a multigrid focusing calculation
            shrink_fac = 1.
            instructions = instructions + "\tfglen {0:.2f} {1:.2f} {2:.2f}\n".format(
                    self.gridSizeX * self.meshSizeX * shrink_fac,
                    self.gridSizeY * self.meshSizeY * shrink_fac,
                    self.gridSizeZ * self.meshSizeZ * shrink_fac)

            # center grid at position given by the center of the pdb
            instructions = instructions + '\tcgcent {x} {y} {z}\n'.format(
                    x=self.centerX,
                    y=self.centerY,
                    z=self.centerZ)
            # center of the fine grid (in a focusing calculation)
            instructions = instructions + '\tfgcent {x} {y} {z}\n'.format(
                    x=self.centerX,
                    y=self.centerY,
                    z=self.centerZ)

            # molecule for which the PBE is to be solved, definition according
            # to the list of molecules given after 'read' starting with 1
            instructions = instructions + '\tmol {0}\n'.format(i)
            # type of equation to be solved, lpbe, npbe, smpbe
            instructions = instructions + '\t{0}\n'.format(self.pbeType)
            instructions = instructions + '\tbcfl {0}\n'.format(self.bcfl)
            # dielectric constant of the protein
            instructions = instructions + '\tpdie {0}\n'.format(self.pdie)
            # dielectric constant of the solvent
            instructions = instructions + '\tsdie {0}\n'.format(self.sdie)
            # bulk concentration of mobile ions, total bulk system of ions has
            # to be electroneutral!
            instructions = instructions + '\tion charge {0} conc {1:.3f} radius {2:.1f}\n'.format(self.ion_pos_charge,
                                                                                                  self.ion_pos_conc,
                                                                                                  self.ion_pos_radius)
            instructions = instructions + '\tion charge {0} conc {1:.3f} radius {2:.1f}\n'.format(self.ion_neg_charge,
                                                                                                  self.ion_neg_conc,
                                                                                                  self.ion_neg_radius)
            # model to construct dielectric constant and ion accessibility
            instructions = instructions + '\tsrfm {0}\n'.format(self.srfm)
            # method by which point charges are mapped to the grid
            instructions = instructions + '\tchgm {0}\n'.format(self.chgm)
            # number of grid points per square angstrom, ignored when srad is
            # 0.0 or srfm is sp12
            instructions = instructions + '\tsdens {0}\n'.format(self.sdens)
            # radius of solvent molceules, 1.4 for water, ignored for sp12
            instructions = instructions + '\tsrad {0}\n'.format(self.srad)
            # size of support (i.e. the rate of change) for spline based
            # surface definitions
            instructions = instructions + '\tswin {0}\n'.format(self.swin)
            # Temperature
            instructions = instructions + '\ttemp {0}\n'.format(self.temp)
            instructions = instructions + '\tcalcenergy {0}\n'.format(self.calcenergy)
            if self.calcforce == 'total':
                instructions = instructions + '\tcalcforce {0}\n'.format(self.calcforce)
            else:
                instructions = instructions + '\tcalcforce no\n'
            instructions = instructions + 'end\n'
        instructions = instructions + 'print elecEnergy complex - ligand - fixed end\n'
        instructions = instructions + 'quit\n'

        with open(file_path, 'w') as f:
            f.write(instructions)

    def write_binding_energy_long(self, file_path):
        '''
        This is the function that writes an infile for binding energy
        calculations.

        Args:
            file_path -> Path to the new infile.

        Returns:
            None.
        '''
        # check if everything is okay:
        self.fixGridSize()
        if self.checkInFile() == False:
            raise NameError('Error in consistency check')
        instructions = ''

        #### read statement
        instructions = instructions + 'read\n'
        # complex
        instructions = instructions + '\tmol pqr {0}\n'.format(self.pqr_path)
        # ligand
        instructions = instructions + '\tmol pqr {0}\n'.format(self.ligand_pqr_path)
        # fixed
        instructions = instructions + '\tmol pqr {0}\n'.format(self.fixed_pqr_path)
        instructions = instructions + 'end\n'

        ## solvated
        solv_names = ["complex_solv", "ligand_solv", "fixed_solv"]
        for i in range(1, 4):
            elec_name = solv_names[i-1]

            #### elec complex
            instructions = instructions + 'elec name {0}\n'.format(elec_name)
            instructions = instructions + '\t{0}\n'.format(self.elec_type)
            # grid size, needs to be: c*2^(l+1) + 1
            instructions = instructions + '\tdime {0} {1} {2}\n'.format(self.gridSizeX,
                                                                        self.gridSizeY,
                                                                        self.gridSizeZ)
            # starting grid length
            instructions = instructions + "\tcglen {0:.2f} {1:.2f} {2:.2f}\n".format(self.gridSizeX * self.meshSizeX,
                                                            self.gridSizeY * self.meshSizeY,
                                                            self.gridSizeZ * self.meshSizeZ)
            # fine mesh domain length in a multigrid focusing calculation
            shrink_fac = 1.
            instructions = instructions + "\tfglen {0:.2f} {1:.2f} {2:.2f}\n".format(
                    self.gridSizeX * self.meshSizeX * shrink_fac,
                    self.gridSizeY * self.meshSizeY * shrink_fac,
                    self.gridSizeZ * self.meshSizeZ * shrink_fac)

            # center grid at position given by the center of the pdb
            instructions = instructions + '\tcgcent {x} {y} {z}\n'.format(
                    x=self.centerX,
                    y=self.centerY,
                    z=self.centerZ)
            # center of the fine grid (in a focusing calculation)
            instructions = instructions + '\tfgcent {x} {y} {z}\n'.format(
                    x=self.centerX,
                    y=self.centerY,
                    z=self.centerZ)

            # molecule for which the PBE is to be solved, definition according
            # to the list of molecules given after 'read' starting with 1
            instructions = instructions + '\tmol {0}\n'.format(i)
            # type of equation to be solved, lpbe, npbe, smpbe
            instructions = instructions + '\t{0}\n'.format(self.pbeType)
            instructions = instructions + '\tbcfl {0}\n'.format(self.bcfl)
            # dielectric constant of the protein
            instructions = instructions + '\tpdie {0}\n'.format(self.pdie)
            # dielectric constant of the solvent
            instructions = instructions + '\tsdie {0}\n'.format(self.sdie)
            # bulk concentration of mobile ions, total bulk system of ions has
            # to be electroneutral!
            instructions = instructions + '\tion charge {0} conc {1:.3f} radius {2:.1f}\n'.format(self.ion_pos_charge,
                                                                                                  self.ion_pos_conc,
                                                                                                  self.ion_pos_radius)
            instructions = instructions + '\tion charge {0} conc {1:.3f} radius {2:.1f}\n'.format(self.ion_neg_charge,
                                                                                                  self.ion_neg_conc,
                                                                                                  self.ion_neg_radius)
            # model to construct dielectric constant and ion accessibility
            instructions = instructions + '\tsrfm {0}\n'.format(self.srfm)
            # method by which point charges are mapped to the grid
            instructions = instructions + '\tchgm {0}\n'.format(self.chgm)
            # number of grid points per square angstrom, ignored when srad is
            # 0.0 or srfm is sp12
            instructions = instructions + '\tsdens {0}\n'.format(self.sdens)
            # radius of solvent molceules, 1.4 for water, ignored for sp12
            instructions = instructions + '\tsrad {0}\n'.format(self.srad)
            # size of support (i.e. the rate of change) for spline based
            # surface definitions
            instructions = instructions + '\tswin {0}\n'.format(self.swin)
            # Temperature
            instructions = instructions + '\ttemp {0}\n'.format(self.temp)
            instructions = instructions + '\tcalcenergy {0}\n'.format(self.calcenergy)
            if self.calcforce == 'total':
                instructions = instructions + '\tcalcforce {0}\n'.format(self.calcforce)
            else:
                instructions = instructions + '\tcalcforce no\n'
            instructions = instructions + 'end\n'

        ## reference
        ref_names = ["complex_ref", "ligand_ref", "fixed_ref"]
        for i in range(1, 4):
            elec_name = ref_names[i-1]

            #### elec complex
            instructions = instructions + 'elec name {0}\n'.format(elec_name)
            instructions = instructions + '\t{0}\n'.format(self.elec_type)
            # grid size, needs to be: c*2^(l+1) + 1
            instructions = instructions + '\tdime {0} {1} {2}\n'.format(self.gridSizeX,
                                                                        self.gridSizeY,
                                                                        self.gridSizeZ)
            # starting grid length
            instructions = instructions + "\tcglen {0:.2f} {1:.2f} {2:.2f}\n".format(self.gridSizeX * self.meshSizeX,
                                                            self.gridSizeY * self.meshSizeY,
                                                            self.gridSizeZ * self.meshSizeZ)
            # fine mesh domain length in a multigrid focusing calculation
            shrink_fac = 1.
            instructions = instructions + "\tfglen {0:.2f} {1:.2f} {2:.2f}\n".format(
                    self.gridSizeX * self.meshSizeX * shrink_fac,
                    self.gridSizeY * self.meshSizeY * shrink_fac,
                    self.gridSizeZ * self.meshSizeZ * shrink_fac)

            # center grid at position given by the center of the pdb
            instructions = instructions + '\tcgcent {x} {y} {z}\n'.format(
                    x=self.centerX,
                    y=self.centerY,
                    z=self.centerZ)
            # center of the fine grid (in a focusing calculation)
            instructions = instructions + '\tfgcent {x} {y} {z}\n'.format(
                    x=self.centerX,
                    y=self.centerY,
                    z=self.centerZ)

            # molecule for which the PBE is to be solved, definition according
            # to the list of molecules given after 'read' starting with 1
            instructions = instructions + '\tmol {0}\n'.format(i)
            # type of equation to be solved, lpbe, npbe, smpbe
            instructions = instructions + '\t{0}\n'.format(self.pbeType)
            instructions = instructions + '\tbcfl {0}\n'.format(self.bcfl)
            # dielectric constant of the protein
            instructions = instructions + '\tpdie {0}\n'.format(self.pdie)
            # dielectric constant of the solvent
            instructions = instructions + '\tsdie {0}\n'.format(self.pdie)
            # bulk concentration of mobile ions, total bulk system of ions has
            # to be electroneutral!
            instructions = instructions + '\tion charge {0} conc {1:.3f} radius {2:.1f}\n'.format(self.ion_pos_charge,
                                                                                                  self.ion_pos_conc,
                                                                                                  self.ion_pos_radius)
            instructions = instructions + '\tion charge {0} conc {1:.3f} radius {2:.1f}\n'.format(self.ion_neg_charge,
                                                                                                  self.ion_neg_conc,
                                                                                                  self.ion_neg_radius)
            # model to construct dielectric constant and ion accessibility
            instructions = instructions + '\tsrfm {0}\n'.format(self.srfm)
            # method by which point charges are mapped to the grid
            instructions = instructions + '\tchgm {0}\n'.format(self.chgm)
            # number of grid points per square angstrom, ignored when srad is
            # 0.0 or srfm is sp12
            instructions = instructions + '\tsdens {0}\n'.format(self.sdens)
            # radius of solvent molceules, 1.4 for water, ignored for sp12
            instructions = instructions + '\tsrad {0}\n'.format(self.srad)
            # size of support (i.e. the rate of change) for spline based
            # surface definitions
            instructions = instructions + '\tswin {0}\n'.format(self.swin)
            # Temperature
            instructions = instructions + '\ttemp {0}\n'.format(self.temp)
            instructions = instructions + '\tcalcenergy {0}\n'.format(self.calcenergy)
            if self.calcforce == 'total':
                instructions = instructions + '\tcalcforce {0}\n'.format(self.calcforce)
            else:
                instructions = instructions + '\tcalcforce no\n'
            instructions = instructions + 'end\n'

#        ref_names = ["complex_ref", "ligand_ref", "fixed_ref"]
        instructions = instructions + 'print elecEnergy {0} - {1} - {2} + {3} - {4} + {5} end\n'.format(
            solv_names[0],ref_names[0],solv_names[1],ref_names[1],solv_names[2],ref_names[2])
        instructions = instructions + 'quit\n'

        with open(file_path, 'w') as f:
            f.write(instructions)

    def checkInFile(self):
        '''
        Checks the information in the infile, java leftover ... pointless!
        '''
        checkStatus = True
        if self.pqr_path == '':
            checkStatus = False
            print("ERROR: No input pqr file specified!")

        elif (self.meshSizeX == 0.0 or self.meshSizeY == 0.0 or
              self.meshSizeZ == 0.0):
            checkStatus = False
            print("ERROR: Illegal mesh sizes!")
        return checkStatus

    def set_options(self, apbs_input_dict):
        '''
        This method accepts a dictionary and sets the elements from the
        dictionary (keys). If there is an element which could not be set, it
        will raise an error.
        The given values of the dictionary are not checked for validation.
        If some options are not given the default values will be used.
        '''
        data_dict = apbs_input_dict.copy()
        if 'elec_type' in data_dict:
            #     mg-auto -> Automatically-configured sequential focusing
            #                multigrid Poisson-Boltzmann calculations.
            #     mg-para -> Automatically-configured parallel focusing
            #                multigrid Poisson-Boltzmann calculations.
            #     mg-manual -> Manually-configured multigrid Poisson-Boltzmann
            #                  calculations.
            #     fe-manual -> Manually-configured adaptive finite element
            #                  Poisson-Boltzmann calculations.
            #     mg-dummy -> Calculations of surface and charge distribution
            #                 properties which do not require solution of the PBE.
            self.elec_type = data_dict['elec_type']
            data_dict.pop('elec_type')

        if 'pbe_type' in data_dict:
            self.pbeType = data_dict['pbe_type']
            data_dict.pop('pbe_type')

        if 'bcfl' in data_dict:
            # type of boundary conditions
            #     sdh -> fast, works best, if boundaries are far away
            #     mdh -> slow, works better then sdh, if boundaries are close
            self.bcfl = data_dict['bcfl']
            data_dict.pop('bcfl')

        if 'pdie' in data_dict:
            # protein interior dielectric constant
            self.pdie = data_dict['pdie']
            data_dict.pop('pdie')

        if 'sdie' in data_dict:
            # solvent dielectric constant
            self.sdie = data_dict['sdie']
            data_dict.pop('sdie')

        if 'nlev' in data_dict:
            # depth of multilevel hierarchy? marked as deprecated
            self.nlev = data_dict['nlev']
            data_dict.pop('nlev')

        if 'srfm' in data_dict:
            # model to construct the dielectric and ion-accessibiity
            #     mol -> dielectric coefficient is defined based on molecular
            #            surface definition
            #     smol -> same as mol, but smoothed with a 9-point harmonic
            #            averaging to reduce sensitivity to grid setup
            #     spl2 -> dielectric and ion-accessibility are defined by a cubic
            #            spline surface, can generate unphysical results
            #     spl4 -> dielectric and ion-accessibility are defined a 7th order
            #            polynomial
            self.srfm = data_dict['srfm']
            data_dict.pop('srfm')

        if 'chgm' in data_dict:
            # method by which point charges (i.e. Dirac delta functions) are mapped
            # to the grid
            #    spl0 -> trilinear interpolation (linera splines), very sensitive
            #            to grid setup
            #    spl2 -> cubic B-spline interpolation, not so sensitive
            #    spl4 -> quintic B-spline interpolation (125 grid points receive
            #            charge density, next-next-nearest neighbors)
            self.chgm = data_dict['chgm']
            data_dict.pop('chgm')

        if 'sdens' in data_dict:
            # number of grid points per square angstrom to use in discontinuos
            # surface constructions
            self.sdens = data_dict['sdens']
            data_dict.pop('sdens')

        if 'srad' in data_dict:
            # radius of the solvent molecules, water = 1.4
            self.srad = data_dict['srad']
            data_dict.pop('srad')

        if 'swin' in data_dict:
            # size of the support (i.e., rate of change) for spline based surface
            # definitions
            self.swin = data_dict['swin']
            data_dict.pop('swin')

        if 'temp' in data_dict:
            # Temperature
            self.temp = data_dict['temp']
            data_dict.pop('temp')

        if('ion_pos_charge' in data_dict and 'ion_pos_conc' in data_dict and
           'ion_pos_radius' in data_dict and 'ion_neg_charge' in data_dict
           and 'ion_neg_conc' in data_dict and 'ion_neg_radius' in data_dict):
            # ion charge -> bulk system must be electroneutral!
            #     charge in e_c
            #     conc in M
            #     radius in A
            self.ion_pos_charge = data_dict['ion_pos_charge']
            self.ion_pos_conc = data_dict['ion_pos_conc']
            self.ion_pos_radius = data_dict['ion_pos_radius']
            self.ion_neg_charge = data_dict['ion_neg_charge']
            self.ion_neg_conc = data_dict['ion_neg_conc']
            self.ion_neg_radius = data_dict['ion_neg_radius']
            data_dict.pop('ion_pos_charge')
            data_dict.pop('ion_pos_conc')
            data_dict.pop('ion_pos_radius')
            data_dict.pop('ion_neg_charge')
            data_dict.pop('ion_neg_conc')
            data_dict.pop('ion_neg_radius')
        elif('ion_pos_charge' in data_dict and 'ion_pos_conc' in data_dict and
             'ion_pos_radius' in data_dict):
            # set positive and negative ions to the same values
            self.ion_pos_charge = data_dict['ion_pos_charge']
            self.ion_pos_conc = data_dict['ion_pos_conc']
            self.ion_pos_radius = data_dict['ion_pos_radius']
            self.ion_neg_charge = -data_dict['ion_pos_charge']
            self.ion_neg_conc = data_dict['ion_pos_conc']
            self.ion_neg_radius = data_dict['ion_pos_radius']
            data_dict.pop('ion_pos_charge')
            data_dict.pop('ion_pos_conc')
            data_dict.pop('ion_pos_radius')

        if 'etol' in data_dict:
            # tolerance for iterations of the pmg partial differential equation
            # solver
            self.etol = data_dict['etol']
            data_dict.pop('etol')

        if 'calcenergy' in data_dict:
            self.calcenergy = data_dict['calcenergy']
            data_dict.pop('calcenergy')

        if 'calcforce' in data_dict:
            self.calcforce = data_dict['calcforce']
            data_dict.pop('calcforce')

        if len(data_dict) != 0:
            for k, v in data_dict.iteritems():
                print(k, v)
            raise AttributeError('Invalid arguments!')


def fix_grid_size(proposed_box_dim, nlev=4):
    '''
    Due to a multilevel approach APBS requires the grid to be of certain sizes.
    More specifically, given a non-zero integer *c* and a multilevel hierarchy
    depth *l*, the APBS grid dimension is calculated as :math:`n = c2^{l+1}+1`.
    The proposed box dimensions will be inflated as necessary.
    
    :param proposed_box_dim: proposed APBS grid dimensions
    :type  proposed_box_dim: tuple(int,int,int)
    :param nlev: depth of the multilevel hierarchy
    :type  nlev: int
    :returns: Valid APBS grid dimensions
    :returntype: :class:`numpy.ndarray[3]`
    
    Examples::
    
        >>> from epitopsy.APBS import fix_grid_size
        >>> # valid APBS dimensions are not affected
        >>> fix_grid_size([33, 65, 129])
        array([ 33,  65, 129])
        >>> # invalid APBS dimensions are inflated
        >>> fix_grid_size([32, 66, 120])
        array([ 33,  97, 129])
    '''
    calc_dime = lambda c, l=4: c * 2 ** (l + 1) + 1
    
    dime = []
    for i in range(len(proposed_box_dim)):
        c = 1
        while proposed_box_dim[i] > calc_dime(c, nlev):
            c += 1
        dime.append(calc_dime(c, nlev))
    
    return np.array(dime)

