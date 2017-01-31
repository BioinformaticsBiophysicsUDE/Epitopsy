__author__     = "Christoph Wilms, Jean-Noel Grad"
__copyright__  = "Copyright 2013, Epitopsy"
__date__       = "2013-03-07"
__credits__    = ["Christoph Wilms", "Jean-Noel Grad"]
__doc__        = """Surface complementarity screening"""



import os
import time
import shutil
import numpy as np

from epitopsy.Structure import PDBFile, PQRFile
from epitopsy.APBS import APBSWrapper
from epitopsy.DXFile import DXBox, DXReader, VDWBox
from epitopsy.scoring.FFTCorrelationScoring import FFT_correlation_scoring
from epitopsy.result.FFT_Result import FFT_Result
from epitopsy.cython import CalculateInteractionEnergyCython
from epitopsy.tools import MathTools
from epitopsy.tools.UtilityClasses import Progress_bar_countdown
from epitopsy.tools.style import style
from epitopsy.tools.calc_vol_and_surface import calc_vol_and_surface


def run_APBS(pdb_path,
             ligands,
             mesh_size    = [1.,1.,1.],
             temp         = 310.,
             ph           = None,
             cubic_box    = True,
             extend       = None,
             box_dim      = None,
             center_pdb   = False,
             box_center   = [0, 0, 0],
             box_type     = ["esp", "vdw"],
             use_pdb2pqr  = True,
             pdb2pqr_argv = None,
             verbose      = True):
    '''
    Compute the electrostatic potential of a protein using
    :class:`APBS.APBSWrapper`, in preparation to a surface complementarity
    screening with :func:`run_SCS`.
    Following files will be produced, with **pdb_path** = ``protein.pdb`` and
    all other arguments to default:
    
    * :file:`protein.in`: the APBS script
    * :file:`protein.pqr`: PQR file of the protein, unless already provided
    * :file:`protein_esp.dx`: electrostatic potential DXBox
    * :file:`protein_vdw.dx`: van der Waals DXBox
    * :file:`io.mc`: APBS log file
    
    By default, the box dimensions (*dime*) are calculated from the protein in
    **pdb_path** and from the diameter of the largest ligand in **ligands**,
    plus 2 Angstroms in the 6 directions. Custom box dimensions can be supplied
    to **box_dim**, but the custom box must be large enough to contain the
    protein + ligand. Compliance to the `APBS dime format
    <http://www.poissonboltzmann.org/docs/apbs-overview/#elec-keyword-dime>`_
    is assessed by :func:`MathTools.fix_grid_size` and the number of grid
    points will be adjusted upwards if necessary.
    
    With **center_pdb**, the protein is ensured to be at the center of the
    DXBox (this might result in a smaller DXBox if **cubic_box** is ``True``). 
    If one wants to compare the electrostatic potential of different proteins
    which have been already aligned, it is better to leave **center_pdb** to
    ``False``.
    
    If **use_pdb2pqr** is ``True``, PDB2PQR is called to generate a PQR file
    from the protein PDB file (optionally using **ph** to change the
    protonation state), otherwise the user should supply a PQR file in
    the working directory, with identical basename but with the extension
    ".pdb" changed to ".pqr".
    
    :param pdb_path: path to the protein PDB file
    :type  pdb_path: str
    :param mesh_size: grid mesh size in all three dimensions (optional),
        default is 1.0 Angstrom
    :type  mesh_size: list
    :param temp: temperature for APBS (optional), default is 310 K
    :type  temp: float
    :param ph: pH for protonation state determination in PDB2PQR (optional),
        default is ``None``
    :type  ph: float
    :param use_pdb2pqr: use pdb2pqr to generate the PQR from **pdb_path**;
        if ``False``, a PQR with the same name, but ending with '.pqr' has
        to be in the same folder
    :type  use_pdb2pqr: bool
    :param pdb2pqr_argv: additional arguments to PDB2PQR (optional), namely
        ["--assign-only"], or ["--assign-only", "--noopt"]
    :type  pdb2pqr_argv: list
    :param center_pdb: if ``True`` the PDB structure will be centered at
        **box_center** (optional), default is ``False``
    :type  center_pdb: bool
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
        default ["esp","vdw"] as they are mandatory for :func:`run_SCS`,
        "smol" can be useful for further processing outside Epitopsy
        (["esp","vdw","smol"])
    :type  box_type: list
    :param verbose: print calculation details on screen if ``True`` (optional)
    :type  verbose: bool

    :returns: ``None``, but files are created in the directory of **pdb_path**
    :raises AttributeError: if the protein in **pdb_path** has a chain break
    
    Example::
    
        >>> from epitopsy import calculate_partition_function3 as cpf
        >>> cpf.run_APBS(pdb_path = "WT.B99990003.pdb", mesh_size = 3*[1.0,],
        ...              ligands = ["8025.pqr"], verbose = True)
        centered pdb:	False
        fixed structure:	WT.B99990003.pdb
        mesh size:	(1.0,1.0,1.0)
        value of extend:	45
        box dim:	(161,161,161)
        Starting APBS ...
        Finished.
    
    '''
    
    if isinstance(ligands, basestring):
        ligands = [ligands]
    pdb_path = os.path.abspath(pdb_path)
    pqr_path = pdb_path.replace('.pdb', '.pqr')
    pdb_name = os.path.basename(pdb_path)
    pqr_name = os.path.basename(pqr_path)
    ligands = [os.path.abspath(x) for x in ligands]
    for filename in [pdb_path] + ligands:
        if not os.path.isfile(filename):
            raise IOError("File {0} cannot be read.".format(filename))
    
    # if PDB file is located outside current working directory, copy it here
    if os.path.abspath(pdb_name) != pdb_path:
        shutil.copy(pdb_path, pdb_name)
    
    # load PDB file
    pdb_struct = PDBFile(pdb_path)
    
    # center PDB at **box_center** if requested
    if center_pdb:
        pdb_coord = pdb_struct.determine_geometric_center()
        new_coord = np.array(box_center)
        pdb_struct.translate(new_coord - pdb_coord)
        pdb_struct.save_to_file(pdb_struct.structure_path)
    
    # calculate DXBox dime from the PDB structure, or use supplied dime
    if not box_dim:
        # extension of the 6 faces with the largest ligand + 2 Angstroms
        if extend is None:
            extend = max([int(PQRFile(ligand).determine_max_diameter()) 
                          for ligand in ligands]) + 2
        box_dim = pdb_struct.get_dxbox_dim(mesh_size, extend, cubic_box)
    else:
        new_box_dim = MathTools.fix_grid_size(box_dim) # check supplied dime
        if np.any(new_box_dim != box_dim):
            if verbose:
                print("fixed grid size {0} -> {1}".format(box_dim, new_box_dim))
            box_dim = new_box_dim
    
    # use PDB2PQR if requested
    if use_pdb2pqr is True or center_pdb is True:
        # APBS wrapper will automatically call PDB2PQR
        pqr_path = None
        pqr_name = None
    elif os.path.abspath(pqr_name) != pqr_path:
        shutil.copy(pqr_path, pqr_name)
    
    # print details
    if verbose:
        print("protein:\t{0}".format(os.path.split(pdb_path)[-1]))
        print("centered:\t{0}".format(center_pdb))
        print("mesh size:\t({0},{1},{2}) Angstroms".format(*mesh_size))
        print("extension:\t{0} Angstroms".format(extend))
        print("box dim:\t({0:.0f},{1:.0f},{2:.0f})".format(*box_dim))
        print("Starting APBS...")
    
    # run APBS
    apbs = APBSWrapper()
    apbs.get_dxbox(pdb_path=pdb_name, mesh_size=mesh_size,
                   pqr_path=pqr_name, box_dim=box_dim,
                   box_center=box_center,
                   box_type=box_type,
                   cubic_box=cubic_box,
                   temperature=temp, ph=ph,
                   pdb2pqr_argv=pdb2pqr_argv,
                   no_return = True)
    
    if verbose:
        print("Done.")


def run_SCS(pdb_path,
            ligand_path,
            APBS_dx_path        = '.',
            number_of_rotations = 150,
            write_top_hits      = False,
            explicit_sampling   = False,
            zipped              = False,
            verbose             = True,
            interior            = -15.,
            no_return           = True):
    '''
    Screen a protein surface with a molecular probe.
    
    :class:`cython.FFTCorrelationScoringCython` is called to compute the FFT
    of the DXfiles and :func:`CalculateInteractionEnergyCython.count_rotations`
    to calculate the free binding energy. If **explicit_sampling** is ``True``,
    :func:`CalculateInteractionEnergyCython.speed_up_brute_force` is used
    instead. If **write_top_hits** is ``True``, a list of the highest scoring
    ligand orientations are written in the working directory.
    Generate several files, including:
    
    * :file:`protein_epi.dx`: the difference in Gibbs free energy in
      units of |kbT| for each grid point (),
    * :file:`protein_mic.dx.gz`: the number of free rotations on each grid
      point
    
    The ESP calculated in APBS is dimensionless (it's divided by kT/|e|).
    We therefore don't have to multiply Phi by <math>\\beta = 1/kT</math>
    nor have to multiply the PQR charges by <math>|e|</math> when computing
    the Maxwell-Boltzmann probability of presence, i.e. <math>e^{-\\beta
    \\cdot \\Phi \\cdot q} = e^{-\\Phi^{APBS} \\cdot q^{PQR}}</math>.
    
    :param pdb_path: path to the protein PDB file
    :type  pdb_path: str
    :param ligand_path: path to the ligand PQR file, if PDB is supplied,
        **ligand_path** extension is changed to ".pqr"
    :type  ligand_path: str
    :param APBS_dx_path: path to the directory where the .dx files are stored
        (optional), by default search in the current directory
    :type  APBS_dx_path: str
    :param number_of_rotations: how many orientations allowed for the ligand
        (optional), default is 150
    :type  number_of_rotations: int
    :param write_top_hits: write the top scoring conformations
    :type  write_top_hits: bool
    :param explicit_sampling: if ``True`` it does not use FFT correlation to
        detect overlapping orientations (optional), default is ``False``; this
        option does not work with **write_top_hits**
    :type  explicit_sampling: bool
    :param zipped: if ``True`` all DXBoxes will be gzipped, if ``False``, the
        energy matrix and the LAS surface are not zipped (optional), default
        is ``False``
    :type  zipped: bool
    :param verbose: print calculation details and progress bar on screen
        (optional), default is ``True``
    :type  verbose: bool
    :param interior: VdW penalty score delta (optional), negative, default -15
    :type  interior: float

    :returns: ``None``, but create 2 files in the current working directory.
    :raises AttributeError: if **ligand_path** has no vdw information
    :raises IOError: if **APBS_dx_path**, **pdb_path** or **ligand_path**
        cannot be read
    :raises ValueError: if **explicit_sampling** is not a boolean
    
    Example::
    
        >>> from epitopsy import calculate_partition_function3 as cpf
        >>> cpf.run_APBS(pdb_path = "WT.B99990003.pdb", mesh_size = 3*[1.0,],
        ...              ligands = ["8025.pqr"], verbose = False)
        >>> cpf.run_SCS(pdb_path = "WT.B99990003.pdb", ligand_path = "8025.pqr")
        This is the setup on Thu Sep 11 14:49:35 2014 ...
        working dir:	/home/user/screening
        non-zero vdw atoms of ligand:	209 / 209
        box dim:	(97,97,97)
        rotations:	150
        FFT sampling:	True
        Starting screening...
        0%                                            100%   time left
        #++++++++++++++++++++++++++++++++++++++++++++++++#   0:00:00        
        Writing grids.
        Writing energy.
        Writing counter matrix.
        Analyzing data.
        total time elapsed: 2.09 min

    '''
    
    # settings
    vdw_correlation_eps = 0.0001
    num_of_best_scores = 3
    pymol_threshold = 0.9e37
    basename = os.path.basename(pdb_path).replace(".pdb", "_{0}.dx")
    path_epi = basename.format("epi")
    path_mic = basename.format("mic") + ".gz"
    path_esp = os.path.join(APBS_dx_path, basename.format("esp"))
    path_vdw = os.path.join(APBS_dx_path, basename.format("vdw"))
    LAS_energy_path = style["LAS_box_file_path"]
    
    # rangecheck
    if not APBS_dx_path or not os.path.isdir(APBS_dx_path):
        raise IOError("APBS_dx_path cannot be read.")
    if not os.path.isfile(pdb_path):
        raise IOError("pdb_path cannot be read.")
    if not os.path.isfile(ligand_path):
        raise IOError("ligand_path cannot be read.")
    
    if zipped:
        if not LAS_energy_path.endswith(".gz"):
            LAS_energy_path += ".gz"
        if not path_epi.endswith(".gz"):
            path_epi += ".gz"
    else:
        LAS_energy_path.rstrip(".gz")
        path_epi.rstrip(".gz")
    
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
    vdwbox.flood() # 1's outside protein, 0's inside
    
    # create a list of independent rotations using a golden section algorithm
    angle_stack = MathTools.get_euler_angles_for_equal_distributed_rotation(number_of_rotations)
    
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
            interaction_energy, counter_matrix = CalculateInteractionEnergyCython.speed_up_brute_force(
                    current_ligand_clone, vdwbox.box, espbox, counter_matrix, 0., interior)
            
            partition_function += np.exp(-interaction_energy)
        
        else:
            # move the ligand to the center (0,0,0)
            current_ligand_clone.translate(-current_ligand_clone.determine_geometric_center())
            
            # snap ligand shape to grid, correlate and shift FFT
            pqr_vdw = current_ligand_clone.snap_vdw_to_box(vdwbox.box_mesh_size,
                                              vdwbox.box_dim, vdwbox.box_offset)
            vdw_correlation = vdw_fft.get_correlation(pqr_vdw)
            vdw_correlation = vdw_fft.shift_fft(vdw_correlation)
            
            # snap ligand charges to grid, correlate and shift FFT
            pqr_esp = current_ligand_clone.snap_esp_to_dxbox(espbox)
            esp_correlation = esp_fft.get_correlation(pqr_esp)
            esp_correlation = esp_fft.shift_fft(esp_correlation)
            
            # energy calculation in a cython script
            interaction_energy, counter_matrix = CalculateInteractionEnergyCython.count_rotations(
                vdw_correlation, esp_correlation, counter_matrix)
            
            # convert binding affinity to equilibrium constant
            partition_function += np.exp(-interaction_energy)
            
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
        " using Epitopsy function run_SCS()",
        "   Interaction energy in kT (negative energies are attractive)",
        "     protein:     {}".format(pdb_path),
        "     ligand:      {}".format(ligand_path),
        "     APBS esp:    {}".format(path_esp),
        "     APBS vdw:    {}".format(path_vdw),
        "     microstates: {}".format(path_mic)])
    energy_box.write(path_epi)
    
    # save counter_matrix
    if verbose:
        print("Writing counter matrix")
    counter_box = DXBox(counter_matrix, espbox.box_mesh_size, espbox.box_offset)
    energy_box.setCommentHeader(["OpenDX file created by {} on {}".format(
                                          os.getenv("USER"), time.ctime()),
        " using Epitopsy function run_SCS()",
        "   Number of available microstates (allowed rotations), integer value"
                                                                    " between",
        "   0 (no rotation) and {} (free rotation)".format(number_of_rotations),
        "     protein:  {}".format(pdb_path),
        "     ligand:   {}".format(ligand_path),
        "     APBS esp: {}".format(path_esp),
        "     APBS vdw: {}".format(path_vdw),
        "     energies: {}".format(path_epi)])
    counter_box.write(path_mic)
    
    # compute binding energies and write to file
    if verbose:
        print("Analyzing data")
    calc_vol_and_surface(pdb_path, 1., zipped = zipped, energy_box = energy_box,
                         counter_box = counter_box, raise_error = False)
    
    if verbose:
        total_time = time.time() - start_time
        print("Total time elapsed: {0} min".format(round(total_time / 60., 2)))
    
    if not no_return:
        return affinity



def calculate_interaction_energy(pdb_path,
                                 ligand_path,
                                 mesh_size,
                                 number_of_rotations = 150,
                                 Temperature         = 310.,
                                 extend              = None,
                                 ph                  = None,
                                 use_pdb2pqr         = True,
                                 cubic_box           = True,
                                 center_pdb          = False,
                                 box_center          = [0, 0, 0],
                                 box_dim             = None,
                                 box_type            = ["esp", "vdw"],
                                 write_top_hits      = False,
                                 explicit_sampling   = False,
                                 zipped              = False,
                                 pdb2pqr_argv        = None):
    '''
    Wrapper for :func:`run_APBS` and :func:`run_SCS`.
    Used for retrocompatibility.
    '''
    
    run_APBS(pdb_path,
             ligands      = [ligand_path],
             mesh_size    = mesh_size,
             temp         = Temperature,
             ph           = ph,
             cubic_box    = cubic_box,
             extend       = extend,
             box_dim      = box_dim,
             center_pdb   = center_pdb,
             box_center   = box_center,
             box_type     = box_type,
             use_pdb2pqr  = use_pdb2pqr,
             pdb2pqr_argv = pdb2pqr_argv)
    
    affinity = run_SCS(pdb_path,
            ligand_path,
            number_of_rotations = number_of_rotations,
            APBS_dx_path        = '.',
            write_top_hits      = write_top_hits,
            explicit_sampling   = explicit_sampling,
            zipped              = zipped,
            no_return           = False)
    
    return affinity


