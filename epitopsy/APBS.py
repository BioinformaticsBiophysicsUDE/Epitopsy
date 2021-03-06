__author__     = "Christoph Wilms, Jean-Noel Grad"
__copyright__  = "Copyright 2012, Epitopsy"
__date__       = "2012-01-10"
__credits__    = ["Christoph Wilms", "Jean-Noel Grad"]

import os
import re
import stat
import subprocess
import numpy as np
import psutil
import epitopsy
from epitopsy.Structure import PDBFile, PQRFile
from epitopsy.DXFile import DXReader


class APBSWrapper(object):
    '''
    Wrapper for the APBS software.

    :param APBSPath: path to the APBS executable
    :type  APBSPath: str
    :param PDB2PQRPath: path to the PDB2PQR executable
    :type  PDB2PQRPath: str

    .. attribute:: APBSPath

        (**str**) Path to the APBS executable.

    .. attribute:: PDB2PQRPath

        (**str**) Path to the PDB2PQR executable.

    .. attribute:: binding_energy

        (**float**) Binding energy.

    .. attribute:: debye_length

        (**list(float)**) Iterable of the Debye length.

    '''

    def __init__(self, APBSPath='apbs', PDB2PQRPath='pdb2pqr'):
        self.APBSPath = APBSPath
        self.PDB2PQRPath = PDB2PQRPath
        self.binding_energy = None
        self.debye_length = []

    def runPDB2PQR(self, pdb_path, pqr_path, force_field='amber',
                   pdb2pqr_argv=None):
        '''
        Call pdb2pqr to replace the b-factor and the occupancy information
        in the pdb file with the charges and the vdw radii.

        :param pdb_path: path to the pdb file
        :type  pdb_path: str
        :param pqr_path: path for the new pqr file
        :type  pqr_path: str
        :param force_field: forcefield from which charges and radii should be
            taken, default is 'amber'
        :type  force_field: str
        :param pdb2pqr_argv: additional arguments for pdb2pqr, must be stored
            in a list, even if there is only one argument (e.g.
            ``['--assign-only']`` or ``['--assign-only', '--noopt']``).
        :type pdb2pqr_argv: list(str)
        '''
        pdb_object = PDBFile(pdb_path)
        pqr_object = pdb_object.get_pqr_structure(pqr_path,
                                                  force_field=force_field,
                                                  pdb2pqr_argv=pdb2pqr_argv)

    def runAPBS(self, apbsInParameters, inFilename, APBS_RAM_threshold=1.):
        '''
        Run APBS.
        
        :param apbsInParameters: all neccessary parameters
        :type  apbsInParameters: :class:`InFile`
        :param inFilename: name of the inputfile
        :type  inFilename: str
        :param APBS_RAM_threshold: how many gigabytes of RAM should remain
           available to the operating system during the APBS calculation
        :type  APBS_RAM_threshold: float
        :raises MemoryError: if not enough free RAM for the calculation to run
        '''
        if not(os.path.exists(apbsInParameters.pqr_path)):
            raise NameError('Error running APBS: File not found: {0}'.format(
                                           apbsInParameters.getPQRFilePath()))
        apbsInParameters.write(inFilename)
        
        # RAM manager
        # coarse grid: at most 256 bytes per grid cell
        # fine grid: at most 224 bytes per grid cell
        mem_required = 256 / 1024.**3 * float(apbsInParameters.gridSizeX) \
                                      * float(apbsInParameters.gridSizeY) \
                                      * float(apbsInParameters.gridSizeZ)
        mem_available = psutil.virtual_memory().available / 1024.**3
        if mem_available < mem_required:
            raise MemoryError('APBS job requires {:.1f} GB of RAM ({:.1f} GB '
                              'available), cannot run the calculation.'
                              .format(mem_required, mem_available))
        elif mem_available < mem_required + APBS_RAM_threshold:
            raise MemoryError('APBS job requires {:.1f} GB of RAM ({:.1f} GB '
                              'available), the calculation might freeze the OS'
                              .format(mem_required, mem_available))
        
        # fix HETATM records in large proteins
        pqr_path = apbsInParameters.pqr_path
        pqr_original = open(pqr_path).read()
        pqr_fixed = re.sub('HETATM(?=\d)', 'ATOM  ', pqr_original)
        if pqr_fixed != pqr_original:
            note = '''\
REMARK   5 NOTE: several "HETATM" records were changed to type "ATOM".
REMARK   5       This change is necessary for APBS version <= 1.4.1 where
REMARK   5       a whitespace character must be present between "HETATM"
REMARK   5       and the atom_number, otherwise the line is ignored.
REMARK   5       This is usually the case for atom_number >= 10000.
REMARK   5'''
            i = 0
            pattern = '(?:\n|^)(REMARK +6 |ATOM  |HETATM)'
            for m in re.finditer(pattern, pqr_fixed):
                i = m.start()
                if pqr_fixed[i] == '\n':
                    i += 1
                break
            pqr_fixed = pqr_fixed[:i] + note + '\n' + pqr_fixed[i:]
            open(pqr_path, 'w').write(pqr_fixed)
        
        # run APBS
        p = subprocess.Popen([self.APBSPath, inFilename],
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        apbs_output, apbs_error = p.communicate()

        apbs_error = apbs_error.split('\n')
        apbs_output = apbs_output.split('\n')
        # there are some error messages (something with fscan) but they
        # are not relevant as they just say that the programm encountered
        # an end of file during the calculation
        for line in apbs_error:
            if line.startswith('Error while parsing input file.'):
                raise AttributeError('Error while parsing input file to APBS.')

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

    def get_dxbox(self, protein_path, mesh_size, no_return=False, **kwargs):
        '''
        This method starts pdb2pqr and apbs and returns the specified boxes.

        :param protein_path: path to the structure
        :type  protein_path: str
        :param mesh_size: grid mesh size in all three dimensions
        :type  mesh_size: tuple(int)
        :param temperature: temperature for APBS (optional), default is 310 K
        :type  temperature: float
        :param ph: pH for protonation state determination in PDB2PQR
            (optional), default is ``None``
        :type  ph: float
        :param pdb2pqr_argv: additional arguments to PDB2PQR (optional), for
            example ``["--assign-only"]`` or ``["--assign-only", "--noopt"]``
        :type  pdb2pqr_argv: list(str)
        :param box_center: APBS box center (optional), default is [0,0,0]
        :type  box_center: list(float)
        :param box_dim: override box dimensions with custom values (optional),
            but values must comply to the APBS *dime* format
        :type  box_dim: list(int)
        :param extend: custom box size extension (optional)
        :type  extend: float
        :param cubic_box: use a cubic box if ``True`` (optional)
        :type  cubic_box: bool
        :param close_boundaries: use multiple Debye-Huckel model to initialize
            potential at the boundaries if ``True`` (optional)
        :type  close_boundaries: bool
        :param box_type: types of boxes APBS should write to disk (optional),
            by default ['esp','vdw']
        :type  box_type: list(str)

        :returns: DXBox objects stored by box types.
        :rtype: dict(str, :class:`epitopsy.DXFile.DXBox`)
        '''
        if not os.path.exists(protein_path):
            raise ValueError('Could not find file "{0}"'.format(protein_path))
        
        box_dim = kwargs.get('box_dim')
        box_center = kwargs.get('box_center')
        extend = kwargs.get('extend')
        box_type = kwargs.get('box_type', ['esp', 'vdw'])
        cubic_box = kwargs.get('cubic_box', True)
        close_boundaries = kwargs.get('close_boundaries', True)
        temperature = kwargs.get('temperature', 310)
        ph = kwargs.get('ph')
        pdb2pqr_argv = kwargs.get('pdb2pqr_argv', [])
        if isinstance(pdb2pqr_argv, basestring):
            pdb2pqr_argv = [pdb2pqr_argv]
        
        # process protein structure
        ext = protein_path.split('.')[-1].lower().rstrip('~')
        if ext == 'pqr':
            pqr_path = protein_path
        elif ext == 'pdb':
            pdb = PDBFile(protein_path)
            pqr_path = protein_path.replace('.pdb', '.pqr')
            if ph is not None:
                pdb2pqr_argv.append('--with-ph={0}'.format(ph))
            pqr = pdb.get_pqr_structure(pqr_path, pdb2pqr_argv = pdb2pqr_argv)
            pqr.save_to_file(pqr_path)
        else:
            raise ValueError('Unknown file extension "{}"'.format(ext))
        
        # generate APBS input file
        template_in = InFile(pqr_path = pqr_path,
                             calculation_type = 'potential',
                             box_mesh_size = mesh_size,
                             box_dim = box_dim,
                             box_center = box_center,
                             extend = extend, cubic_box = cubic_box,
                             box_type = box_type)

        template_in.temp = temperature

        if close_boundaries:
            template_in.bcfl = 'mdh'
        else:
            template_in.bcfl = 'sdh'

        self.runAPBS(template_in, '{0}.in'.format(pqr_path.replace('.pqr', '')))
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
                           fixed_pqr_path, box_mesh_size, extend=None,
                           **kwargs):
        '''
        Run APBS on the complex, protein and ligand structures and extract the
        binding energy (in kJ/mol) from the APBS output, according to the
        following equation:

            :math:`\\Delta_\\text{bind,solv}G = \
                        \\Delta_\\text{solv}G_\\text{complex} \
                      - \\Delta_\\text{solv}G_\\text{protein} \
                      - \\Delta_\\text{solv}G_\\text{ligand}`

        The total polar solvation energy includes both the reaction field and
        the coulombic contributions. Positive values are favorable.

        :param complex_pqr_path: path of the complex
        :type  complex_pqr_path: str
        :param ligand_pqr_path: path of the isolated ligand
        :type  ligand_pqr_path: str
        :param fixed_pqr_path: path of the isolated protein
        :type  fixed_pqr_path: str
        :param box_mesh_size: specify the grid mesh size (Angstroms)
        :type  box_mesh_size: array
        :param extend: increase the box dimensions along each axis (Angstroms)
        :type  extend: float
        :param \*\*kwargs: optional keywords for :meth:`InFile.set_options`
        :type  \*\*kwargs: any

        :returns: Binding energy in kJ/mol.
        :rtype: float
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

        apbsInParameters.set_options(kwargs)

        self.runAPBS(apbsInParameters, inFilename)
        return self.binding_energy

    def get_binding_energy_long(self, complex_pqr_path, ligand_pqr_path,
                           fixed_pqr_path, box_mesh_size, extend = None,
                           **kwargs):
        '''
        Run APBS on the complex, protein and ligand structures and extract the
        total binding energy (in kJ/mol) from the APBS output, according to
        the following equations:

            :math:`\\Delta_\\text{bind,solv}G = \
                \\left( \\Delta_\\text{solv}G_\\text{complex} \
                      - \\Delta G^\\text{ref}_\\text{complex} \\right) \
              - \\left( \\Delta_\\text{solv}G_\\text{protein} \
                      - \\Delta G^\\text{ref}_\\text{protein} \\right) \
              - \\left( \\Delta_\\text{solv}G_\\text{ligand} \
                      - \\Delta G^\\text{ref}_\\text{ligand} \\right)`

            :math:`\\Delta_\\text{coul}G = \
               \\dfrac{ \\Delta_\\text{coul}G_\\text{complex} \
                      - \\Delta_\\text{coul}G_\\text{protein} \
                      - \\Delta_\\text{coul}G_\\text{ligand} } \
                      { \\epsilon_{p} }`

            :math:`\\Delta_\\text{bind}G = \\Delta_\\text{bind,solv}G \
                                         + \\Delta_\\text{coul}G`

        Where
        :math:`\\Delta_\\text{bind,solv}G` is the binding energy in water,
        :math:`\\Delta_\\text{coul}G` is the electrostatic contribution,
        :math:`\\Delta_\\text{bind}G` is the total binding energy.
        :math:`\\Delta_\\text{solv}G_i` is the solvation energy of species *i*
        in water (:math:`\\epsilon_s = 79`),
        :math:`\\Delta G^\\text{ref}_i` is the free energy of species *i* in a
        reference medium (:math:`\epsilon_p = 2`),
        :math:`\\Delta_\\text{coul}G_i` is the coulombic free energy of
        species *i* in a reference medium (:math:`\epsilon_p = 2`).

        Positive values of :math:`\\Delta_\\text{bind}G` are favorable.

        :param complex_pqr_path: path of the complex
        :type  complex_pqr_path: str
        :param ligand_pqr_path: path of the isolated ligand
        :type  ligand_pqr_path: str
        :param fixed_pqr_path: path of the isolated protein
        :type  fixed_pqr_path: str
        :param box_mesh_size: specify the grid mesh size (Angstroms)
        :type  box_mesh_size: array
        :param extend: increase the box dimensions along each axis (Angstroms)
        :type  extend: float
        :param \*\*kwargs: optional keywords for :meth:`InFile.set_options`
        :type  \*\*kwargs: any

        :returns: Binding energy in kJ/mol.
        :rtype: float
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

        apbsInParameters.set_options(kwargs)

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
        # print("solv: {0}, coul: {1} -> bind: {2}".format(solvation_energy, coulomb_energy, total_energy))

        return total_energy

    def get_dissociation_energy(self, complex_pqr_path, ligand_pqr_path,
                                fixed_pqr_path, box_mesh_size, extend = None,
                                **kwargs):
        '''
        Run APBS on the complex, protein and ligand structures and extract the
        dissociation energy (in kJ/mol) from the APBS output, according to the
        following equation:

            :math:`\\Delta_\\text{diss,solv}G = - \\left( \
               \\Delta_\\text{solv}G_\\text{complex} \
             - \\Delta_\\text{solv}G_\\text{protein} \
             - \\Delta_\\text{solv}G_\\text{ligand} \\right)`

        The total polar solvation energy includes both the reaction field and
        the coulombic contributions. Positive values are favorable.

        :param complex_pqr_path: path of the complex
        :type  complex_pqr_path: str
        :param ligand_pqr_path: path of the isolated ligand
        :type  ligand_pqr_path: str
        :param fixed_pqr_path: path of the isolated protein
        :type  fixed_pqr_path: str
        :param box_mesh_size: specify the grid mesh size (Angstroms)
        :type  box_mesh_size: array
        :param extend: increase the box dimensions along each axis (Angstroms)
        :type  extend: float
        :param \*\*kwargs: optional keywords for :meth:`InFile.set_options`
        :type  \*\*kwargs: any

        :returns: Dissociation energy in kJ/mol.
        :rtype: float
        '''
        inFilename = complex_pqr_path.replace('.pqr', '.in')
        apbsInParameters = InFile(pqr_path = complex_pqr_path,
                                  calculation_type = 'binding_energy',
                                  box_mesh_size = box_mesh_size,
                                  extend = extend,
                                  cubic_box = True,
                                  ligand_pqr_path = ligand_pqr_path,
                                  fixed_pqr_path = fixed_pqr_path)

        apbsInParameters.set_options(kwargs)

        self.runAPBS(apbsInParameters, inFilename)
        return -self.binding_energy

def get_coulomb_energy(protein_pqr_path):
    '''
    :param protein_pqr_path: path to the pqr of the protein
    :type protein_pqr_path: str

    :returns: Coulomb energy in kJ/mol.
    :rtype: float
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


class InFile(object):
    '''
    APBS input file.
    Only the most basic parameters are covered here, intended for the
    generation of basic grids holding information of electrostatic
    potential and/or the van der Waals surface.

    :param pqr_path: pqr path, from which the potential energy should be
                calculated or the path of the complex for which the binding
                energy should be calculated
    :type  pqr_path: str
    :param calculation_type: potential or binding_energy
    :type  calculation_type: str

    :param box_dim: override box dimensions with custom values (optional), but
       values must comply to the APBS dime format
    :type  box_dim: tuple(float)
    :param box_center: center of the APBS box
    :type  box_center: tuple(float)
    :param mesh_size: grid mesh size in all three dimensions (Angstroms)
    :type  mesh_size: tuple(float)
    :param extend: extend box size along each axis (Angstroms),
       overriding **box_dim**
    :type  extend: float
    :param cubic_box: use a cubic box if ``True``
    :type  cubic_box: bool
    :param box_type: types of boxes APBS should write to disk:
       ``'esp'``: electrostatic potential,
       ``'vdw'``: van der waals based solvent accessibility,
       ``'smol'``: solvent accessibility,
       ``'ndens'``: total mobile ion number density in units of M,
       ``'qdens'``: total mobile charge density in units of e_c M
    :type  box_type: list(str)
    :param ligand_pqr_path: path to ligand PQR file if the calculation is
       a binding energy calculation
    :type  ligand_pqr_path: str
    :param fixed_pqr_path: path to the protein PQR file if the calculation is
       a binding energy calculation
    :type  fixed_pqr_path: str
    :raises ValueError: if **ligand_pqr_path** or **fixed_pqr_path** is
       ``None``
    :raises ValueError: if unknown value in argument **calculation_type**
    :raises AttributeError: if both **box_dim** and **extend** parameters are
       provided

    .. attribute:: pqr_path

        pqr path, from which the potential energy should be calculated or the
        path of the complex for which the binding energy should be calculated

    .. attribute:: calculation_type

        potential or binding_energy

    .. attribute:: box_mesh_size

        mesh size

    .. attribute:: box_dim

        special dimension can be supplied (although they might be
        fixed to work properly with apbs), if ``None`` is given it
        will be calculated from the size of the protein.

    .. attribute:: box_center

        center of the box, if ``None`` is given if will be set to the
        geometric center of the protein

    .. attribute:: extend

        can only be used when no box_dim is supplied, extends the box
        size by the given amount (in Angstroem)

    .. attribute:: cubic_box

        determine wheter it is a cubic box or not, in the case of
        a cubic box the largest dimension is used for all dimensions

    .. attribute:: box_type

        * esp: electrostatic potential
        * vdw: van der waals based solvent accessibility
        * smol: solvent accessibility
        * ndens: total mobile ion number density in units of M
        * qdens: total mobile charge density in units of e_c M

    .. attribute:: ligand_pqr_path

        path to the ligand PQR file if the calculation is a binding energy
        calculation

    .. attribute:: fixed_pqr_path

        path to the protein PQR file if the calculation is a binding energy
        calculation

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
        # type of calculation
        if calculation_type in (self.binding_energy, self.binding_energy_long):
            self.calculation_type = calculation_type
            if ligand_pqr_path is None or fixed_pqr_path is None:
                raise ValueError('For a binding energy calculation a ligand '
                                 'and a fixed protein are necessary!')
        elif calculation_type == self.potential:
            self.calculation_type = self.potential
        else:
            raise ValueError("Unknown Argument for 'calculation_type': '{0}'!"
                             .format(calculation_type))

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
            box_center = PQRFile(self.pqr_path).determine_geometric_center()
        self.setGridCenter(box_center)

        self.setMeshSize(box_mesh_size)

        if box_dim is not None and extend is not None:
            raise AttributeError('Input error! It does not make sense to '
                         'provide the box dimension AND an extend parameter!')

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
                box_dim = np.ceil(box_dim / np.array(box_mesh_size)) \
                          + 2 * extend / np.array(box_mesh_size)
            else:
                box_dim = np.ceil(box_dim / np.array(box_mesh_size))

        self.setGridSize(box_dim)
        self.fixGridSize()

    def setGridCenter(self, center):
        '''
        :setter: Sets grid geometrical center (Angstroms)
        :type: tuple(float)
        '''
        self.centerX = center[0]
        self.centerY = center[1]
        self.centerZ = center[2]


    def setGridSize(self, size):
        '''
        :setter: Sets grid size (Angstroms)
        :type: tuple(float)
        '''
        self.gridSizeX = size[0]
        self.gridSizeY = size[1]
        self.gridSizeZ = size[2]


    def setMeshSize(self, meshSize):
        '''
        :setter: Sets mesh size (Angstroms)
        :type: tuple(float)
        '''
        self.meshSizeX = meshSize[0]
        self.meshSizeY = meshSize[1]
        self.meshSizeZ = meshSize[2]

    def generateFromPDB(self, pdb, padding, cubicBox):
        self.generateFromPDB2(pdb, padding, 0.0, cubicBox)

    def generateFromPDB2(self, pdb, padding, minDiameter, cubicBox):
        if(self.meshSizeX == 0.0 or self.meshSizeY == 0.0 or self.meshSizeX == 0.0):
            raise NameError("ERROR: Cannot calculate grid size mesh sizes < =  0.0!")
        # * TODO: use method to determine geometric center of protein-associated atoms only!
        # * Same goes for determination of coordinate extremes. Only use protein-associated atoms!
        # ???
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

    def fixGridSize(self):
        '''
        See :func:`fix_grid_size`.
        '''
        fixed_size = fix_grid_size([self.gridSizeX,
                self.gridSizeY, self.gridSizeZ])
        self.gridSizeX = fixed_size[0]
        self.gridSizeY = fixed_size[1]
        self.gridSizeZ = fixed_size[2]

    def write(self, file_path):
        '''
        Write an APBS input file.

        :param file_path: path to the input file
        :type  file_path: str
        '''
        if self.calculation_type == self.potential:
            self._write_potential(file_path)
        elif self.calculation_type == self.binding_energy:
            self._write_binding_energy(file_path)
        elif self.calculation_type == self.binding_energy_long:
            self._write_binding_energy_long(file_path)


    def _write_potential(self, file_path):
        '''
        Write an APBS input file for a potential calculation.

        :param file_path: path to the input file
        :type  file_path: str
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
        instructions = instructions + "\tcglen {0:.2f} {1:.2f} {2:.2f}\n".format(
                (self.gridSizeX - 1) * self.meshSizeX,
                (self.gridSizeY - 1) * self.meshSizeY,
                (self.gridSizeZ - 1) * self.meshSizeZ)
        # fine mesh domain length in a multigrid focusing calculation
        shrink_fac = 1.
        instructions = instructions + "\tfglen {0:.2f} {1:.2f} {2:.2f}\n".format(
                (self.gridSizeX - 1) * self.meshSizeX * shrink_fac,
                (self.gridSizeY - 1) * self.meshSizeY * shrink_fac,
                (self.gridSizeZ - 1) * self.meshSizeZ * shrink_fac)
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

    def _write_binding_energy(self, file_path):
        '''
        Write an APBS input file for an energy calculation.

        :param file_path: path to the input file
        :type  file_path: str
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
            instructions = instructions + "\tcglen {0:.2f} {1:.2f} {2:.2f}\n".format(
                    (self.gridSizeX - 1) * self.meshSizeX,
                    (self.gridSizeY - 1) * self.meshSizeY,
                    (self.gridSizeZ - 1) * self.meshSizeZ)
            # fine mesh domain length in a multigrid focusing calculation
            shrink_fac = 1.
            instructions = instructions + "\tfglen {0:.2f} {1:.2f} {2:.2f}\n".format(
                    (self.gridSizeX - 1) * self.meshSizeX * shrink_fac,
                    (self.gridSizeY - 1) * self.meshSizeY * shrink_fac,
                    (self.gridSizeZ - 1) * self.meshSizeZ * shrink_fac)

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

    def _write_binding_energy_long(self, file_path):
        '''
        Write an APBS input file for an energy calculation.

        :param file_path: path to the input file
        :type  file_path: str
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
            instructions = instructions + "\tcglen {0:.2f} {1:.2f} {2:.2f}\n".format(
                    (self.gridSizeX - 1) * self.meshSizeX,
                    (self.gridSizeY - 1) * self.meshSizeY,
                    (self.gridSizeZ - 1) * self.meshSizeZ)
            # fine mesh domain length in a multigrid focusing calculation
            shrink_fac = 1.
            instructions = instructions + "\tfglen {0:.2f} {1:.2f} {2:.2f}\n".format(
                    (self.gridSizeX - 1) * self.meshSizeX * shrink_fac,
                    (self.gridSizeY - 1) * self.meshSizeY * shrink_fac,
                    (self.gridSizeZ - 1) * self.meshSizeZ * shrink_fac)

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
            instructions = instructions + '\tdime {0} {1} {2}\n'.format(
                    (self.gridSizeX - 1) * self.meshSizeX,
                    (self.gridSizeY - 1) * self.meshSizeY,
                    (self.gridSizeZ - 1) * self.meshSizeZ)
            # starting grid length
            instructions = instructions + "\tcglen {0:.2f} {1:.2f} {2:.2f}\n".format(self.gridSizeX * self.meshSizeX,
                                                            self.gridSizeY * self.meshSizeY,
                                                            self.gridSizeZ * self.meshSizeZ)
            # fine mesh domain length in a multigrid focusing calculation
            shrink_fac = 1.
            instructions = instructions + "\tfglen {0:.2f} {1:.2f} {2:.2f}\n".format(
                    (self.gridSizeX - 1) * self.meshSizeX * shrink_fac,
                    (self.gridSizeY - 1) * self.meshSizeY * shrink_fac,
                    (self.gridSizeZ - 1) * self.meshSizeZ * shrink_fac)

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

        :param apbs_input_dict: APBS instructions as key/value pairs
        :type  apbs_input_dict: dict(str,str)

        :raises AttributeError: if unknown key
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
    :type  proposed_box_dim: tuple(int)
    :param nlev: depth of the multilevel hierarchy
    :type  nlev: int
    :returns: Valid APBS grid dimensions
    :rtype: :class:`np.array[3]`

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

