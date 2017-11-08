__author__     = "Christoph Wilms, Jean-Noel Grad"
__copyright__  = "Copyright 2013, Epitopsy"
__date__       = "2013-03-07"
__credits__    = ["Christoph Wilms", "Jean-Noel Grad"]
__doc__        = '''Surface and charge complementarity screening'''

import os
import time
import shutil
import numpy as np

from epitopsy.Structure import PDBFile, PQRFile
from epitopsy import APBS
from epitopsy.DXFile import DXBox, DXReader, VDWBox
from epitopsy.cython import interaction_explicit_sampling
from epitopsy.tools import MathTools
from epitopsy.tools.UtilityClasses import Progress_bar_countdown
from epitopsy.EnergyGrid.FFT import FFT_correlation_scoring
from epitopsy.EnergyGrid.analysis import calc_vol_and_surface, FFT_Result


def electrostatics(protein,
             ligands,
             mesh_size    = (1.,1.,1.),
             temp         = 310.,
             ph           = None,
             cubic_box    = False,
             extend       = None,
             box_dim      = None,
             transform    = 'center+rotate',
             align        = None,
             box_center   = (0, 0, 0),
             box_type     = ('esp', 'vdw'),
             pdb2pqr_argv = None,
             verbose      = True):
    '''
    Compute the electrostatic potential of a protein in preparation for a
    correlation with :func:`scan`.
    
    PDB2PQR will convert the PDB file to a PQR file while preserving common
    ions (calcium, zinc, sulfate, see
    :meth:`Structure.PDBFile.get_pqr_structure`). For uncommon ions
    and non-standard residues, please supply a custom PQR file in **protein**.
    APBS will carry out the electrostatic calculation 
    (see :class:`APBS.APBSWrapper`).
    
    By default, the box dimensions (*dime*) are calculated from the structure
    in **protein** and from the diameter of the largest ligand in **ligands**,
    plus 2 Angstroms in the 6 directions. Custom box dimensions can be supplied
    to **box_dim**, but the custom box must be large enough to contain the
    protein + ligand. Compliance to the `APBS dime format
    <http://www.poissonboltzmann.org/docs/apbs-overview/#elec-keyword-dime>`_
    is assessed by :func:`APBS.fix_grid_size` and the number of grid
    points will be adjusted upwards if necessary.
    
    :param protein: path to the protein PDB or PQR file
    :type  protein: str
    :param mesh_size: grid mesh size in all three dimensions (optional),
        default is 1.0 Angstrom
    :type  mesh_size: list
    :param temp: temperature for APBS (optional), default is 310 K
    :type  temp: float
    :param ph: pH for protonation state determination in PDB2PQR (optional),
        default is ``None``
    :type  ph: float
    :param pdb2pqr_argv: additional arguments to PDB2PQR (optional), namely
        ["--assign-only"], or ["--assign-only", "--noopt"]
    :type  pdb2pqr_argv: list
    :param transform: structure transformation: `'center+rotate'` to center and
        rotate, `'center'` to center, ``None`` to leave the structure intact,
        `'superimpose'` to superimpose onto another protein provided in
        **align**
    :type  transform: str
    :param align: structure onto which to superimpose **protein** (optional),
       can be a tuple giving the structure filename and the atom names to use
       for alignment, only the structure filename or another protein object
    :type  align: str or :class:`Structure.PDBFile` or
       tuple(str, tuple(str))
    :param box_center: center of the APBS box (optional), default is [0,0,0]
    :type  box_center: array
    :param box_dim: override box dimensions with custom values (optional), but
        values must comply to the APBS dime format
    :type  box_dim: list
    :param extend: override box size extension with custom values (optional),
        by default the 6 sides are extended by ligand max diameter + 2 Angstroms
    :type  extend: float
    :param cubic_box: use a cubic box if ``True`` (optional)
    :type  cubic_box: bool
    :param box_type: types of boxes APBS should write to disk (optional), by
        default ['esp', 'vdw'] as they are mandatory for :func:`scan`
    :type  box_type: list(str)
    :param verbose: print calculation details to screen if ``True`` (optional)
    :type  verbose: bool
    
    Example::
    
        >>> from epitopsy import EnergyGrid
        >>> m = 0.8
        >>> protein_path = '5b1l-protein.pdb'
        >>> ligand_path  = '5b1l-DNA-fragment.pqr'
        >>> EnergyGrid.electrostatics(protein_path, [ligand_path],
        ...                           mesh_size=3*[m,], verbose=True)
        protein:	5b1l-protein.pdb
        centered:	False
        mesh size:	(0.8,0.8,0.8) Angstroms
        extension:	32 Angstroms
        box dim:	(225,225,225)
        Running APBS...
        Done.
    
    will produce:
    
    * :file:`5b1l-protein.pqr`: PQR file of the protein
    * :file:`5b1l-protein.in`: APBS input file
    * :file:`5b1l-protein_esp.dx`: electrostatic potential DXBox
    * :file:`5b1l-protein_vdw.dx`: van der Waals DXBox
    * :file:`io.mc`: APBS log file
    '''
    
    if isinstance(ligands, basestring):
        ligands = [ligands]
    for filename in [protein] + ligands:
        if not os.path.isfile(filename):
            raise IOError("File {0} cannot be read.".format(filename))
    
    # if PDB file is located outside current working directory, copy it here
    protein_local = os.path.basename(protein)
    if protein_local != os.path.relpath(protein):
        shutil.copy(protein, protein_local)
    
    ligands = [os.path.abspath(x) for x in ligands]
    protein = os.path.abspath(protein_local)
    
    # load PDB/PQR structure
    ext = protein.split('.')[-1].lower().rstrip('~')
    if ext == 'pdb':
        pro_struct = PDBFile(protein)
    elif ext == 'pqr':
        pro_struct = PQRFile(protein)
    else:
        raise ValueError('Unknown file extension "{}"'.format(ext))
    
    # center, rotate or align protein structure
    if transform in ('super', 'superimpose', 'align'):
        if not (isinstance(align, list) or isinstance(align, tuple)):
            align_template = align
            align_atoms = None
        elif len(align) == 1:
            align_template = align[0]
            align_atoms = None
        else:
            align_template = align[0]
            align_atoms = align[1]
        if isinstance(align_template, basestring):
            ext = align_template.split('.')[-1].lower().rstrip('~')
            if ext == 'pdb':
                ref = PDBFile(align_template)
            elif ext == 'pqr':
                ref = PQRFile(align_template)
            else:
                raise ValueError('Unknown file extension "{}"'.format(ext))
        else:
            ref = align_template
        if not (isinstance(ref, PDBFile) or isinstance(ref, PQRFile)):
            raise TypeError('Argument align should be str or PDBFile')
        pro_struct.superimpose_self_onto_given_pdb(ref, align_atoms)
        pro_struct.save_to_file(pro_struct.structure_path)
    elif transform in ('center', 'center+rotate'):
        if transform == 'center+rotate':
            pro_struct.apply_PCA_projection()
        pro_coord = pro_struct.determine_geometric_center()
        new_coord = np.array(box_center)
        pro_struct.translate(new_coord - pro_coord)
        pro_struct.save_to_file(pro_struct.structure_path)
    
    # calculate DXBox dime from the PDB structure, or use supplied dime
    if not box_dim:
        # extension of the 6 faces with the largest ligand + 2 Angstroms
        if extend is None:
            extend = max([int(PQRFile(ligand).determine_max_diameter()) 
                          for ligand in ligands]) + 2
        box_dim = pro_struct.get_dxbox_dim(mesh_size, extend, cubic_box)
    else:
        new_box_dim = APBS.fix_grid_size(box_dim) # check supplied dime
        if np.any(new_box_dim != box_dim):
            if verbose:
                print('fixed box_dim: {0} -> {1}'.format(box_dim, new_box_dim))
            box_dim = new_box_dim
    
    # print details
    if verbose:
        print("protein:\t{0}".format(os.path.split(protein)[-1]))
        print("transform:\t{0}".format(transform))
        print("mesh size:\t({0},{1},{2}) Angstroms".format(*mesh_size))
        print("extension:\t{0} Angstroms".format(extend))
        print("box dim:\t({0:.0f},{1:.0f},{2:.0f})".format(*box_dim))
        print("Running APBS...")
    
    # run APBS
    apbs = APBS.APBSWrapper()
    apbs.get_dxbox(protein, mesh_size,
                   box_dim=box_dim,
                   box_center=box_center,
                   box_type=box_type,
                   cubic_box=cubic_box,
                   temperature=temp, ph=ph,
                   pdb2pqr_argv=pdb2pqr_argv,
                   no_return=True)
    
    if verbose:
        print("Done.")


def scan(pdb_path,
             ligand_path,
             APBS_dx_path        = ".",
             number_of_rotations = 150,
             write_top_hits      = False,
             explicit_sampling   = False,
             zipped              = ("mic",),
             verbose             = True,
             interior            = -15.,
             protein_conc        = 0.001,
             no_return           = True,
             normalization_mode  = None,
             flood               = 'xy'):
    '''
    Surface and Charge Complementarity Screening: screen a protein surface with
    a molecular probe using the Fast Fourier Transform.
    Generate several files, including:
    
    * :file:`protein_epi.dx`: the difference in Gibbs free energy in
      units of |kbT| for each grid point,
    * :file:`protein_mic.dx.gz`: the number of free rotations on each grid
      point (available microstates)
    
    The ESP calculated in APBS (:math:`\\Phi`) is dimensionless (units of
    :math:`k_{\\text{B}}T/|e|`). We therefore don't have to multiply
    :math:`\\Phi` by :math:`\\beta = 1/k_{\\text{B}}T` nor have to multiply 
    the PQR charges (:math:`q`) by :math:`|e|` when computing the
    Maxwell-Boltzmann probability of presence, i.e. :math:`\\Delta G =
    e^{-\\Phi \\cdot q}`.
    
    :param pdb_path: path to the protein PDB file
    :type  pdb_path: str
    :param ligand_path: path to the ligand PQR file, if PDB is supplied,
        **ligand_path** extension is changed to ".pqr"
    :type  ligand_path: str
    :param APBS_dx_path: path to the directory where the .dx files are stored
        (optional), by default search in the current directory
    :type  APBS_dx_path: str
    :param number_of_rotations: how many ligand rotations to sample
    :type  number_of_rotations: int
    :param write_top_hits: write the top scoring conformations if ``True``
    :type  write_top_hits: bool
    :param explicit_sampling: if ``True`` it does not use FFT correlation to
        detect overlapping orientations (optional), default is ``False``; this
        option does not work with **write_top_hits**
    :type  explicit_sampling: bool
    :param zipped: gzip all DXBoxes if ``True``, none if ``False``, use "epi"
        or "mic" to select individual DXBoxes to gzip
    :type  zipped: bool or tuple
    :param verbose: print calculation details and progress bar to screen if
        ``True``
    :type  verbose: bool
    :type  protein_conc: protein concentration in mol/L (for energy summary),
       default is 1 mM
    :type  protein_conc: float
    :param interior: VdW penalty score delta (optional), must be negative
    :type  interior: float
    :param flood: default is 'xy' to flood by plane, use 'xyz' to start
       from a corner
    :type  flood: str

    :raises AttributeError: if **ligand_path** has no vdw information
    :raises IOError: if **APBS_dx_path**, **pdb_path** or **ligand_path**
        cannot be read
    :raises ValueError: if **explicit_sampling** is not a boolean
    
    Example::
    
        >>> from epitopsy import EnergyGrid
        >>> m = 0.8
        >>> protein_path = '5b1l-protein.pdb'
        >>> ligand_path  = '5b1l-DNA-fragment.pqr'
        >>> EnergyGrid.electrostatics(protein_path, [ligand_path],
        ...                           mesh_size=3*[m,],
        ...                           center_pdb=False, verbose=True)
        protein:	5b1l-protein.pdb
        centered:	False
        mesh size:	(0.8,0.8,0.8) Angstroms
        extension:	32 Angstroms
        box dim:	(225,225,225)
        Running APBS...
        Done.
        >>> EnergyGrid.scan(protein_path, ligand_path, APBS_dx_path='.')
        This is the setup on Tue May 16 21:27:15 2017...
        non-zero vdw atoms of ligand:	441 / 445
        box dim:	(225,225,225)
        rotations:	150
        FFT sampling:	True
        Start screening...
        0%                                            100%   time left
        #+++++++++++++++++++   DONE   +++++++++++++++++++#   0:00:00        
        Writing energy
        Writing counter matrix
        Analyzing data
    
    '''
    
    # settings
    vdw_correlation_eps = 0.0001
    num_of_best_scores = 3
    pymol_threshold = 0.9e37
    basename = os.path.basename(pdb_path).replace(".pdb", "_{0}.dx")
    path_epi = basename.format("epi")
    path_mic = basename.format("mic")
    path_esp = basename.format("esp")
    path_vdw = basename.format("vdw")
    if APBS_dx_path != '.':
        path_esp = os.path.realpath(os.path.join(APBS_dx_path,
                                                 os.path.basename(path_esp)))
        path_vdw = os.path.realpath(os.path.join(APBS_dx_path,
                                                 os.path.basename(path_vdw)))
    
    # rangecheck
    if not APBS_dx_path or not os.path.isdir(APBS_dx_path):
        raise IOError("APBS_dx_path cannot be read.")
    if not os.path.isfile(pdb_path):
        raise IOError("pdb_path cannot be read.")
    if not os.path.isfile(ligand_path):
        raise IOError("ligand_path cannot be read.")
    
    # select files to compress
    if isinstance(zipped, basestring):
        zipped = (zipped,)
    if isinstance(zipped, list) or isinstance(zipped, tuple):
        if "epi" in zipped:
            path_epi += ".gz"
        if "mic" in zipped:
            path_mic += ".gz"
    elif zipped:
        path_epi = path_epi + ".gz"
        path_mic = path_mic + ".gz"
    
    # load the ligand and count how many atoms have non-zero radius
    lig = PQRFile(ligand_path.replace(".pdb",".pqr"))
    vdw_radii_list = [atom.get_info_2() != 0 for atom in lig.get_all_atoms()]
    if sum(vdw_radii_list) == 0:
        raise AttributeError("The ligand has no van der Waals radius "
                             "information! Check PDB2PQR log file.")
    
    # remove old files if any
    if os.path.exists(path_epi):
        os.remove(path_epi)
    if os.path.exists(path_mic):
        os.remove(path_mic)
    
    # read the APBS ESP (dimensionless, already divided by kT/|e|)
    espbox = DXReader().parse(path_esp, "esp")
    # read the APBS VDW volume, 1's in the solvent and 0's inside protein
    vdwbox = DXReader().parse(path_vdw, "vdw")
    
    # create a list of independent rotations using a golden section algorithm
    # (warning: with the Gonzalez algorithm, N + 1 coordinates may be returned)
    angle_stack = MathTools.get_euler_angles_for_equal_distributed_rotation(number_of_rotations)
    number_of_rotations = len(angle_stack)
    
    # print details
    if verbose:
        print("This is the setup on {0}...".format(time.ctime()))
        print("non-zero vdw atoms of ligand:\t{0} / {1}".format(
                                  sum(vdw_radii_list), len(vdw_radii_list)))
        print("box dim:\t({0:.0f},{1:.0f},{2:.0f})".format(*vdwbox.box_dim))
        print("rotations:\t{0}".format(number_of_rotations))
        print("FFT sampling:\t{0}".format(not explicit_sampling))
    
    # flood closed cavities in the VDW box to prevent the ligand from fitting
    # inside them (physically inaccessible)
    vdwbox.flood(method=flood) # 1's outside protein, 0's inside
    
    # set protein core to a negative value: overlap between the ligand and the
    # protein is forbidden; vdw is now 0's outside, 1 on surface, -15 inside
    vdwbox.prepare_for_geometric_matching(interior = interior)
    
    # zero out ESP inside the protein (no correlation, avoid FFT artifacts)
    espbox.box[np.nonzero(vdwbox.box == interior)] = 0.
    
    # additional matrices and DXBoxes
    partition_function = np.zeros(espbox.box.shape) # Maxwell-Boltzmann stat.
    counter_matrix = np.zeros(espbox.box.shape)     # number of microstates
    vdw_fft = FFT_correlation_scoring(vdwbox.box)   # FFT of the protein VDW
    esp_fft = FFT_correlation_scoring(espbox.box)   # FFT of the protein ESP
    
    # store position the highest scoring ligand poses
    if write_top_hits:
        shape_results = FFT_Result()
        esp_results = FFT_Result()
        combined_result = FFT_Result()
    
    if verbose:
        print("Start screening...")
        progress_bar = Progress_bar_countdown(len(angle_stack), countdown=True, refresh=5)
        start_time = time.time()
    
    # start screening
    for angle_element in angle_stack:
        # get angles
        phi, theta, psi = angle_element
        # clone and rotate the ligand
        current_ligand_clone = PQRFile(ligand_path)
        current_ligand_clone.translate_origin_and_rotate(phi, theta, psi)
        
        # do the energy calculation in a cython script if requested
        if explicit_sampling == True:
            microstates_energy, counter_matrix = (interaction_explicit_sampling
                .explicit_sampling(current_ligand_clone, vdwbox.box, espbox,
                                   counter_matrix, 0., interior))
            
            partition_function += np.exp(-microstates_energy)
        
        else:
            # move ligand center to the box center
            current_ligand_clone.translate(-current_ligand_clone.determine_geometric_center())
            #current_ligand_clone.translate(vdwbox.box_center)
            
            # snap ligand shape to grid and perform correlation
            pqr_vdw = current_ligand_clone.snap_vdw_to_box(vdwbox.box_mesh_size,
                                              vdwbox.box_dim, vdwbox.box_offset)
            vdw_correlation = vdw_fft.get_correlation(pqr_vdw)
            vdw_correlation = vdw_fft.shift_fft(vdw_correlation)
            
            # snap ligand charges to grid and perform correlation
            pqr_esp = current_ligand_clone.snap_esp_to_dxbox(espbox)
            esp_correlation = esp_fft.get_correlation(pqr_esp)
            esp_correlation = esp_fft.shift_fft(esp_correlation)
            
            # increment partition function with microstates energy
            heaviside_mask = (vdw_correlation >= 0 - vdw_correlation_eps)
            counter_matrix += heaviside_mask
            partition_function += np.exp(-esp_correlation * heaviside_mask)
            
            # find best position of the ligand if requested
            if write_top_hits:
                # store the best scoring positions
                shape_results.find_scores(vdw_correlation, phi, theta, psi, num_of_best_scores, vdwbox)
                # omit inaccessible positions
                esp_correlation[np.nonzero(vdw_correlation < vdw_correlation_eps)] = 0
                # find the lowest energies, negative because 'find_scores' looks for the highest energy
                esp_results.find_scores(-esp_correlation, phi, theta, psi, num_of_best_scores, espbox)
                vdw_correlation[np.nonzero(vdw_correlation < 0)] = 0
                # set positive, i.e. repelling, esp to 0
                esp_correlation[np.nonzero(esp_correlation > 0)] = 0
                # normalize them and add them together
                vdw_correlation = vdw_correlation / np.amax(vdw_correlation)
                esp_correlation = esp_correlation / np.amin(esp_correlation)
                combined_distance = np.sqrt(vdw_correlation ** 2 + esp_correlation ** 2)
                combined_result.find_scores(combined_distance, phi, theta, psi, num_of_best_scores, espbox)
        
        if verbose:
            progress_bar.add()
    
    # end of calculation, update progress bar
    if verbose:
        progress_bar.terminate()
    
    # compute binding affinity in units of kbT
    probability = partition_function / float(number_of_rotations) # normalize
    #probability = partition_function / np.maximum(counter_matrix,1) # bad idea
    zero_indices = np.nonzero(probability == 0)   # flag zeros (log undefined)
    probability[zero_indices] = np.e              # zeros -> dummy values
    affinity = -np.log(probability)               # compute binding affinity
    affinity[zero_indices] = 0.                   # remove dummy values
    affinity[np.nonzero(np.isnan(affinity))] = 0  # remove NaN
    affinity[np.nonzero(counter_matrix == 0)] = 0 # no rotation -> no affinity
    affinity[np.nonzero(vdwbox.box)] = 0.         # no affinity inside protein
    
    if write_top_hits:
        # write docking results to disk
        shape_results.write_to_disk("shape_docking.txt")
        esp_results.write_to_disk("esp_docking.txt")
        combined_result.write_to_disk("combined_docking.txt")
        if verbose:
            print("Writing docked poses")
        # create docked structures
        ligand_pdb_path = ligand_path.replace(".pqr", ".pdb")
        shape_results.make_pdb_results(ligand_pdb_path, "pdb_pool", "shape", 10)
        esp_results.make_pdb_results(ligand_pdb_path, "pdb_pool", "esp", 10)
        combined_result.make_pdb_results(ligand_pdb_path, "pdb_pool", "combined", 10)
    
    # PyMOL crashes for values exceeding x < -1e37 or x > 1e37, so these
    # values have to be set manually to the maximal allowed value
    affinity[np.nonzero(affinity >  pymol_threshold)] =  pymol_threshold
    affinity[np.nonzero(affinity < -pymol_threshold)] = -pymol_threshold
    energy_box = DXBox(affinity, espbox.box_mesh_size, espbox.box_offset)
    
    # write DXBoxes to disk
    if verbose:
        print("Writing energy")
    energy_box.setCommentHeader([" OpenDX file created by {} on {}".format(
                                          os.getenv("USER"), time.ctime()),
        " using function Epitopsy.EnergyGrid.scan()",
        "   Interaction energy in kT (negative energies are attractive)",
        "     protein:     {}".format(os.path.relpath(pdb_path)),
        "     ligand:      {}".format(os.path.relpath(ligand_path)),
        "     APBS esp:    {}".format(os.path.relpath(path_esp)),
        "     APBS vdw:    {}".format(os.path.relpath(path_vdw)),
        "     microstates: {}".format(os.path.relpath(path_mic))])
    energy_box.write(path_epi)
    
    # save counter_matrix
    if verbose:
        print("Writing counter matrix")
    counter_box = DXBox(counter_matrix, espbox.box_mesh_size, espbox.box_offset)
    counter_box.setCommentHeader([" OpenDX file created by {} on {}".format(
                                          os.getenv("USER"), time.ctime()),
        " using function Epitopsy.EnergyGrid.scan()",
        "   Number of available microstates (allowed rotations), integer value"
                                                                    " between",
        "   0 (no rotation) and {} (free rotation)".format(number_of_rotations),
        "     protein:  {}".format(os.path.relpath(pdb_path)),
        "     ligand:   {}".format(os.path.relpath(ligand_path)),
        "     APBS esp: {}".format(os.path.relpath(path_esp)),
        "     APBS vdw: {}".format(os.path.relpath(path_vdw)),
        "     energies: {}".format(os.path.relpath(path_epi))])
    counter_box.write(path_mic)
    
    # compute binding energies and write to file
    if verbose:
        print("Analyzing data")
    calc_vol_and_surface(1, energy_box, counter_box, raise_error = False,
                         conc = protein_conc)
    
    if verbose:
        total_time = time.time() - start_time
        print("Total time elapsed: {0} min".format(round(total_time / 60., 2)))
    
    if not no_return:
        return affinity


def scan_multiconformational(protein_paths, ligand_paths,
                             elec_kwargs=None, scan_kwargs=None,
                             transform_on_first='center+rotate'):
    '''
    Compute energy grids for multiple protein and ligand conformations.
    
    Merged energies are only meaningful when the protein structures have been
    superimposed. By default the first structure is centered and rotated, while
    subsequent structures are superimposed on the first. Make sure the atom
    order is identical in all proteins, otherwise superimposition will fail.
    If you have already carried out the superimposition outside Epitopsy,
    please set **transform_on_first** to ``None``.
    
    :param protein_paths: paths to the protein PDB files
    :type  protein_paths: list(str) or str
    :param ligand_paths: paths to the ligand PQR files
    :type  ligand_paths: list(str) or str
    :param elec_kwargs: optional arguments for :func:`electrostatics`
    :type  elec_kwargs: dict
    :param scan_kwargs: optional arguments for :func:`scan`
    :type  scan_kwargs: dict
    :param transform_on_first: structure transformation: `'center+rotate'` to
        center and rotate, `'center'` to center, ``None`` to leave the
        structure intact
    :type  transform_on_first: str
    
    Example::
    
        >>> from epitopsy import EnergyGrid
        >>> protein_paths = ['protein1.pdb', 'protein2.pdb']
        >>> ligand_paths = ['ligand1.pqr', 'ligand2.pqr']
        >>> elec_kwargs = {'mesh_size': 3 * [0.8, ], 'verbose': False}
        >>> scan_kwargs = {'verbose': False}
        >>> EnergyGrid.scan_multiconformational(protein_paths, ligand_paths,
        ...                                     elec_kwargs=elec_kwargs,
        ...                                     scan_kwargs=scan_kwargs)
    
    The output files will be stored on disk following this structure:
    
    .. code-block:: none
    
        ./protein1/
          protein1.pdb
          protein1.pqr
          protein1_esp.dx
          protein1_vdw.dx
          ligand1/
            ligand1.pqr
            protein1_epi.dx
            protein1_mic.dx.gz
          ligand2/
            ligand2.pqr
            protein1_epi.dx
            protein1_mic.dx.gz
        ./protein2/
          protein2.pdb
          protein2.pqr
          protein2_esp.dx
          protein2_vdw.dx
          ligand1/
            ligand1.pqr
            protein2_epi.dx
            protein2_mic.dx.gz
          ligand2/
            ligand2.pqr
            protein2_epi.dx
            protein2_mic.dx.gz
    
    '''
    if transform_on_first not in ('center+rotate', 'center', None):
        raise ValueError('Unknown transform "{}"'.format(transform_on_first))
    # keyword arguments
    if elec_kwargs is None:
        elec_kwargs = {}
    if scan_kwargs is None:
        scan_kwargs = {}
    if len(protein_paths) >= 2:
        # with 2 proteins or more, the 'transform' keyword is important
        if transform_on_first:
            for key in ('transform', 'align'):
                if key in elec_kwargs:
                    raise KeyError('Cannot forward argument "{}" to function '
                                   'electrostatics(): use argument "transform_'
                                   'on_first" of function scan_multiconformati'
                                   'onal()'.format(key))
        else:
            if 'transform' not in elec_kwargs:
                elec_kwargs['transform'] = None
                print('WARNING: "transform" in elec_kwargs set to None')
    # file paths
    if isinstance(protein_paths, str):
        protein_paths = [protein_paths]
    if isinstance(ligand_paths, str):
        ligand_paths = [ligand_paths]
    protein_paths = [os.path.realpath(path) for path in protein_paths]
    ligand_paths  = [os.path.realpath(path) for path in  ligand_paths]
    reference_protein = None
    
    # compute APBS and energy grids
    for i, protein_path in enumerate(protein_paths):
        # create a subfolder for the protein
        dirname = os.path.basename(protein_path).rstrip('~')[:-4]
        if os.path.isdir(dirname): # erase folder if it already exists
            shutil.rmtree(dirname)
        os.makedirs(dirname)
        os.chdir(dirname)
        # copy the PDB/PQR file
        shutil.copyfile(protein_path, os.path.basename(protein_path))
        protein_path = os.path.realpath(os.path.basename(protein_path))
        # handle structure transformation
        kwargs = elec_kwargs.copy()
        if transform_on_first:
            # first structure is transformed, others are superimposed
            if i == 0:
                kwargs['transform'] = transform_on_first
            else:
                kwargs['transform'] = 'superimpose'
                kwargs['align'] = reference_protein
        # compute the electrostatic grid
        electrostatics(protein_path, ligand_paths, **kwargs)
        if transform_on_first and i == 0:
            reference_protein = protein_path
        # carry out energy calculations
        for ligand_path in ligand_paths:
            # create a subfolder for the ligand
            dirname = os.path.basename(ligand_path).split('.')[0]
            if not os.path.isdir(dirname):
                os.makedirs(dirname)
            os.chdir(dirname)
            # scan protein
            scan(protein_path, ligand_path, APBS_dx_path='..', **scan_kwargs)
            os.chdir('..')
        os.chdir('..')


def merge_multiconformational(protein_paths, ligand_paths,
                              protein_weights=None, ligand_weights=None):
    '''
    Compute energy grids for multiple protein and ligand conformations.
    
    :param protein_paths: paths to the protein PDB files
    :type  protein_paths: list(str) or str
    :param ligand_paths: paths to the ligand PQR files
    :type  ligand_paths: list(str) or str
    :param protein_weights: custom weights for the proteins, in the same order
       as in **protein_paths**, by default each protein has the same weight
    :type  protein_weights: list(float)
    :param ligand_weights: custom weights for the ligands, in the same order
       as in **ligand_paths**, by default each ligand has the same weight
    :type  ligand_weights: list(float)

    Example::
    
        >>> from epitopsy import EnergyGrid
        >>> protein_paths = ['protein1.pdb', 'protein2.pdb']
        >>> ligand_paths = ['ligand1.pqr', 'ligand2.pqr']
        >>> electrostatics_kwargs = {'mesh_size': 3 * [0.8, ], 'verbose': False}
        >>> scan_kwargs = {'verbose': False}
        >>> EnergyGrid.scan_multiconformational(protein_paths, ligand_paths,
        ...                                     electrostatics_kwargs,\
 scan_kwargs)
        >>> EnergyGrid.merge_multiconformational(protein_paths, ligand_paths)
    
    Two new files will be created in the current directory, `merge_epi.dx` and
    `merge_mic.dx.gz`.
    The multiconformation scan assumes by default the same weight for every
    protein/ligand combination. You may change this behavior by defining custom
    weights in the function arguments. If you need to recompute the merged grid
    with different weights in the future, it is not necessary to call this
    function again, as all grids would be recalculated. Simply call
    :func:`merge` with the new weights. The exact syntax is provided in the
    merged file ``merge_epi.dx``:
    
    .. code-block:: none
    
        $> head merge_epi.dx
        # OpenDX file created by grad on Tue May 16 20:21:45 2017
        # using the following function call:
        #   EnergyGrid.merge(["1krn/Ahx-open/1krn_epi.dx",
        #                     "1krn/Ahx-folded/1krn_epi.dx",
        #                     "2pk4/Ahx-open/2pk4_epi.dx",
        #                     "2pk4/Ahx-folded/2pk4_epi.dx",
        #                     "4duu/Ahx-open/4duu_epi.dx",
        #                     "4duu/Ahx-folded/4duu_epi.dx"],
        #                     weights=[1, 1, 1, 1, 1, 1])

    Simply copy-paste the code a Python terminal, without the ``#`` symbols and
    with the new weights to overwrite ``merge_epi.dx``::
    
        >>> from epitopsy import EnergyGrid
        >>> EnergyGrid.merge(["1krn/Ahx-open/1krn_epi.dx",
        ...                   "1krn/Ahx-folded/1krn_epi.dx",
        ...                   "2pk4/Ahx-open/2pk4_epi.dx",
        ...                   "2pk4/Ahx-folded/2pk4_epi.dx",
        ...                   "4duu/Ahx-open/4duu_epi.dx",
        ...                   "4duu/Ahx-folded/4duu_epi.dx"],
        ...                   weights=[2, 2, 5, 5, 1, 1])
    
    '''
    dxbox_paths = []
    dxbox_weights = []
    
    # file paths
    if isinstance(protein_paths, str):
        protein_paths = [protein_paths]
    if isinstance(ligand_paths, str):
        ligand_paths = [ligand_paths]
    protein_paths = [os.path.realpath(path) for path in protein_paths]
    ligand_paths = [os.path.realpath(path) for path in ligand_paths]
    
    # protein weights
    if isinstance(protein_weights, int) or isinstance(protein_weights, float):
        if len(protein_paths) > 1:
            protein_weights = len(protein_paths) * [protein_weights]
        else:
            protein_weights = [protein_weights]
    elif protein_weights is None:
        protein_weights = len(protein_paths) * [1,]
    
    # ligand weights
    if isinstance(ligand_weights, int) or isinstance(ligand_weights, float):
        if len(ligand_paths) > 1:
            ligand_weights = len(ligand_paths) * [ligand_weights]
        else:
            ligand_weights = [ligand_weights]
    elif ligand_weights is None:
        ligand_weights = len(ligand_paths) * [1,]
    
    # compute grid weights
    for protein_path, protein_weight in zip(protein_paths, protein_weights):
        protein_dirname = os.path.basename(protein_path).split('.')[0]
        for ligand_path, ligand_weight in zip(ligand_paths, ligand_weights):
            ligand_dirname = os.path.basename(ligand_path).split('.')[0]
            dxbox_weights.append(protein_weight * ligand_weight)
            dxbox_paths.append(os.path.join(protein_dirname, ligand_dirname,
                                            os.path.basename(protein_path)
                                            .replace('.pdb', '_epi.dx')))
    
    # merge grids
    merge(dxbox_paths, dxbox_weights)
    

def merge(energy_paths, weights=None, output_fmt='merge_{}.dx'):
    '''
    Merge multiple energetic maps using the Maxwell-Boltzmann formula:
    :math:`E^{\\text{avg}} = -\\log\\left(\\sum_i w_i e^{E_i}\\right)` with 
    :math:`w_i = e^{-U_i}` the weights. For example, in a MD simulation where
    the ligand has conformation A in 200 frames and conformation B in 300
    frames, the weights :math:`w_i` are 200 and 300.
    
    :param energy_paths: paths to the energy DXboxes
    :type  energy_paths: list(str)
    :param weights: Maxwell-Boltzmann weights, default is 1 for each box
    :type  weights: list(float)
    :param output_fmt: path to the output DXBoxes, with format
    :type  output_fmt: str
    
    Example::
    
        >>> from epitopsy import EnergyGrid
        >>> EnergyGrid.merge(["1krn/Ahx-open/1krn_epi.dx",
        ...                   "1krn/Ahx-folded/1krn_epi.dx",
        ...                   "2pk4/Ahx-open/2pk4_epi.dx",
        ...                   "2pk4/Ahx-folded/2pk4_epi.dx",
        ...                   "4duu/Ahx-open/4duu_epi.dx",
        ...                   "4duu/Ahx-folded/4duu_epi.dx"],
        ...                  weights=[1, 1, 1, 1, 1, 1])

    '''
    if not weights:
        weights = len(energy_paths) * [1,]
    merge_epi = output_fmt.format('epi')
    merge_mic = output_fmt.format('mic')
    if merge_epi == merge_mic:
        raise ValueError('arg. output_fmt is incorrect, please read the doc')
    
    # prepare DXBox comment
    comments = ['Output of function epitopsy.EnergyGrid.calculation.merge(',
                '  ["{}"'.format(os.path.relpath(energy_paths[0]))]
    for path in energy_paths[1:]:
        comments[-1] = comments[-1] + ','
        comments.append('   "{}"'.format(os.path.relpath(path)))
    if not weights:
        comments[-1] = comments[-1] + '])'
    else:
        comments[-1] = comments[-1] + '],'
        comments.append('  weights=[{}])'.format(
                        ', '.join('{:.3f}'.format(w) for w in weights)))
    
    weights = np.array(weights, dtype=float) / np.sum(weights)
    
    # merge energies
    summary = []
    dxb = DXReader().parse(energy_paths[0], 'esp')
    box = weights[0] * np.exp(-dxb.box)
    for filename, weight in zip(energy_paths, weights)[1:]:
        dxb = DXReader().parse(filename, 'esp')
        box += weight * np.exp(-dxb.box)
    dxb.box = -np.log(box)
    dxb.setCommentHeader(comments)
    dxb.write(merge_epi)
    
    # find microstate DXBoxes
    microstate_paths = []
    for i in range(len(energy_paths)):
        path = energy_paths[i].replace('_epi.dx', '_mic.dx')
        if not os.path.isfile(path):
            if os.path.isfile(path + '.gz'):
                microstate_paths.append(path + '.gz')
            else:
                raise IOError('Cannot read file "{}"'.format(path))
        else:
            microstate_paths.append(path)
    
    # merge microstates
    dxb = DXReader().parse(microstate_paths[0], 'vdw')
    box = dxb.box
    #box = weights[0] * dxb.box
    for filename, weight in zip(microstate_paths, weights)[1:]:
        dxb = DXReader().parse(filename, 'vdw')
        box += dxb.box
        #box += weight * dxb.box
    #box /= sum(weights)
    dxb.box = box
    dxb.setCommentHeader(comments)
    dxb.write(merge_mic)


def Calculate_Ligand_Interaction_Energy(pdb_path, ligand_path, mesh_size):
    '''
    This function calculates the Energy for the position of the ligand.
    The Complex has to be Centered first! Then take the ligand create a pqr
    file and use this function.
    '''
    kb = 1.3806504e-23
    Temperature = 310
    # load pdb
    pdb_struct = PDBFile(pdb_path)
    '''
    # center it for fft-correlation
    pdb_struct.center()
    pdb_struct.writeToFile(pdb_struct.PDBFilename)
    '''
    # to display the results we need a centered pdb, because the calculations
    # are performed with a centered ligand
    ligand_pdb_path = ligand_path.replace('.pqr', '.pdb')
    ligand_pdb_path_center = ligand_pdb_path.replace('.pdb', '_centered.pdb')
    ligand_pdb_struct = PDBFile(ligand_pdb_path)
    ligand_pdb_struct.center()
    ligand_pdb_struct.writeToFile(ligand_pdb_path_center)
    # load pqr
    ligand_struct = PQRFile(ligand_path)
    ligand_struct.read_pqr_structure()
    '''
    Run electrostatic calculations:
    '''
    apbs = APBS.APBSWrapper("apbs", "pdb2pqr")
    template_in = APBS.InFile("", "", "", mesh_size)
    '''
    Padding is some kind of extension, so i need to pad it by
    the length of the pqr-structure and add 2, because snapping the 
    structure to the box, may enlarge the ligand by 1 unit.
    '''
    padding = int(ligand_struct.determineMaxDiameter() / max(mesh_size) + 2)

    template_in.generateFromPDB(pdb_struct, padding, True)
    template_in.setTemp(Temperature)

    apbs.runPDB2PQR(pdb_path, pdb_path[:-4] + '.pqr')
    
    template_in.setPQRFilePath(pdb_path[:-4] + '.pqr')
    template_in.setOutSurfaceDXPath(pdb_path[:-4] + "_vdw")
    template_in.setOutPotentialDXPath(pdb_path[:-4] + "_esp")
    apbs.runAPBS(template_in, pdb_path[:-4] + '.in')
    '''
    Read the grids as boxes, no need to use the grid class, because i 
    need arrays:
    '''
    dx_esp = DXReader().parse(pdb_path[:-4] + "_esp-PE0.dx", DXReader().ESP,
                              mesh_size)
    espbox = dx_esp.getBox()
    dx_vdw = DXReader().parse(pdb_path[:-4] + "_vdw-PE0.dx", DXReader().VDW,
                              mesh_size)
    vdwbox = dx_vdw.getBox()
    print("Read apbs calculations!")
    
    vdwbox.flood()
    print("Flooded the vdw structure!")

    vdwbox.prepare_for_geometric_matching(interior = -15)
    
    
    # get a clone of the ligand, center it, store the old coordinates
    # and find the energy of the old position
    ligand_clone = ligand_struct.clone()
    shift_vector = ligand_struct.determineCenterOfMass()
    shift_vector_box = vdwbox.transform_real_to_box_space(shift_vector)
    ligand_clone.translate(-shift_vector)
    pqr_vdw, pqr_esp = ligand_clone.snap_to_box(vdwbox.getDimensions(),
                                                vdwbox.getMeshSize(),
                                                vdwbox.getOffset(), True)
    
    shape_scoring = FFT_correlation_scoring(vdwbox.box.astype(float))
    esp_scoring = FFT_correlation_scoring(espbox.box)
    
    shape_results = FFT_Result()
    esp_results = FFT_Result()
    num_of_best_scores = 5
    
    phi = 0
    theta = 0
    psi = 0
    
    # calculate the fft correlation
    shape_correlation = shape_scoring.get_correlation(pqr_vdw.astype(float))
    esp_correlation = esp_scoring.get_correlation(pqr_esp)
    
    # convert to kbT
    esp_correlation = esp_correlation / (kb * Temperature)
    
    # store the best scoring positions
    shape_results.find_scores(shape_correlation, phi, theta, psi,
                              num_of_best_scores, vdwbox)
    # omit positions, that are not accessible
    esp_correlation[np.nonzero(shape_correlation < -0.0001)] = 0
    # find the lowest energies, '-' because 'find_scores' looks
    # for the highest energy
    esp_results.find_scores(-esp_correlation, phi, theta, psi,
                            num_of_best_scores, espbox)
    
    # shift fft
    esp_correlation = esp_scoring.shift_fft(esp_correlation)
    
    # write docking results to disk
    shape_results.write_to_disk('shape_docking.txt')
    esp_results.write_to_disk('esp_docking.txt')
    
    # create docked structures
    # this method needs a pdb at the same positon as the ligand!
    shape_results.make_pdb_results(ligand_pdb_path_center, 'pdb_pool',
                                   'shape', num_of_best_scores)
    esp_results.make_pdb_results(ligand_pdb_path_center, 'pdb_pool',
                                'esp', num_of_best_scores)
    
    ligand_score = esp_correlation[shift_vector_box[0],
                                   shift_vector_box[1],
                                   shift_vector_box[2]]
    print('Score of the ligand at {0} is {1}!'.format(shift_vector,
                                                      ligand_score))


