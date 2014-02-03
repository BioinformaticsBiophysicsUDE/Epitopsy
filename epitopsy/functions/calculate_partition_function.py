"""
Created on March 7, 2013

@author: Christoph Wilms
"""
import os
import time
import numpy as np

from epitopsy.Structure import PDBFile, PQRFile
from epitopsy.APBS import APBSWrapper
from epitopsy.DXFile import DXBox
from epitopsy.scoring.FFTCorrelationScoring import FFT_correlation_scoring
from epitopsy.result.FFT_Result import FFT_Result
from epitopsy.cython import CalculateInteractionEnergyCython
from epitopsy.tools import MathTools
from epitopsy.tools.UtilityClasses import Progress_bar
from epitopsy.tools.style import style
from epitopsy.tools.calc_vol_and_surface import calc_vol_and_surface

def calculate_interaction_energy(pdb_path, ligand_path, mesh_size,
        number_of_rotations=150,
        Temperature=310.,
        extend=None,
        ph=None,
        use_pdb2pqr=True,
        cubic_box=True,
        center_pdb=False,
        box_center=[0, 0, 0],
        box_dim=None,
        box_type=["esp", "vdw"],
        write_top_hits=False,
        explicit_sampling=False,
        zipped=False,
        pdb2pqr_argv=None):
    '''
    This function calculates the difference in gibbs free energy in units
    of kbT.

    probability of all states at position \vec{r}:
    p(r) = 1/Z Sum_i_N * exp(-E_i/kb/T)

    probability of all states at position \vec{r} in solution:
    p(r_solv) = 1/Z Sum_i_N * exp(-E_i/kb/T) = N / Z,
    with E_i = 0 for all rotations

    K = p(r) / p(r_solv) = Sum_i_N * exp(-E_i/kb/T) / N

    gibbs free energy:
    Delta G = - kb T ln(K)

    Args:
        pdb_path -> path to the fixed pdb file
        ligand_path -> path to the .pqr file of the ligand. the .pdb file has
                        to be in the same folder
        mesh_size -> grid mesh size in all three dimensions [mx,my,mz]
        number_of_rotations -> how many different orientations for the ligand?
        Temperature -> temperature for the electrostatic potential calculations
        extend -> extend the box dimensions by this factor
        ph -> Can be used to calculate the pqr for a specific ph, default is
                None
        use_pdb2pqr -> use pdb2pqr to calculate the pqr of pdb_path, if False
                        a pqr with the same name, but ending with '.pqr' has
                        to be in the same folder
        cubic_box -> cubic box or not, hence True or False
        center_pdb -> if True the pdb will be centered at box_center,
                    default is False
        box_center -> center of the apbs calculations, default is [0,0,0]
        box_dim -> dimension of the box, default is None and it is calculated
                    automatically
        box_type -> types of boxes APBS should write to disk. We need at least
                    ["esp","vdw"], "smol" is optional -> ["esp","vdw","smol"]
        write_top_hits -> write the top scoring conformations
        explicit_sampling -> if True it does not use FFT correlation for
                            the detection of overlapping orienations.
                            This option does not work with write_top_hits.
                            default is False.
        zipped -> if True all dxboxes will be gzipped, else the energy
                    matrix and the energy matrix with the LAS surface
                    are not zipped. default is False.

    Returns:
        The energy.
    '''
    ## settings
    shape_correlation_cutoff = 0.0001
    num_of_best_scores = 3
#    np.seterr('raise')
#    kB = 1.3806504e-23
    pymol_threshold = 0.9e37
    energy_path = style["energy_box_file_path"]
    LAS_energy_path = style["LAS_box_file_path"]
    counter_matrix_path = style["counter_box_file_path"]

    wd = os.getcwd()

    if zipped:
        if not LAS_energy_path.endswith(".gz"):
            LAS_energy_path += ".gz"

        if not energy_path.endswith(".gz"):
            energy_path += ".gz"
    else:
        LAS_energy_path.rstrip(".gz")
        energy_path.rstrip(".gz")
    
    if not explicit_sampling in [True,False]:
        raise ValueError("Unkown input for 'explicit_sampling' = {0}".format(explicit_sampling))
    if not use_pdb2pqr in [True,False]:
        raise ValueError("Unkown input for 'use_pdb2pqr' = {0}".format(use_pdb2pqr))

    # remove old energy files, otherwise it may cause errors when the
    # calculation terminates and an old energy file is still there
    # (been there ...)
    if os.path.exists(energy_path):
        os.remove(energy_path)
    # same for the counter_matrix
    if os.path.exists(counter_matrix_path):
        os.remove(counter_matrix_path)

    start_time = time.time()
    # read the ligand
    ligand_struct = PQRFile(ligand_path)

    vdw_radii_list = []
    # check vdw
    for atom in ligand_struct.get_all_atoms():
        if atom.get_info_2() != 0: # atom.get_info_2() returns the vdw radius
            vdw_radii_list.append(1)
        else:
            vdw_radii_list.append(0)

    if sum(vdw_radii_list) == 0:
        raise AttributeError("The ligand has no vdw information!")

    '''
    The parameter <extend> can be used to extend the box in each dimension.
    If extend is none, we will use the maximal diameter of the ligand to
    extend the box, because otherwise not all rotation will fit in the box.
    '''
    if extend is None:
        extend = int(ligand_struct.determine_max_diameter() + 2)

    '''
    Read the pdb file. If box dimensions have been supplied they will be used,
    but one has to make sure, that the box is big enough to contain the
    protein. If none box dimensions have been supplied, they will be
    calculated from the protein.
    '''
    pdb_struct = PDBFile(pdb_path)
    if pdb_struct.contains_chain_break():
        raise AttributeError("The supplied protein structure has a chain break!")

    '''
    It could be useful to center the pdb to (0,0,0). This way the protein
    is definetly at the center of the box. If one wants to compare the
    potential of different structurally aligned proteins, it is better to have
    a fixed grid and use the position of the structural alignment. Otherwise
    the potentials will be shifted by the translation to (0,0,0).
    '''
    if center_pdb:
        # center the structure at the given box_center
        pdb_coord = pdb_struct.determine_geometric_center()
        new_coord = np.array(box_center)
        pdb_struct.translate(new_coord - pdb_coord)
        pdb_struct.save_to_file(pdb_struct.structure_path)

    if box_dim is None:
        # calculate the box dimensions from the pdb structure
        box_dim = pdb_struct.get_dxbox_dim(mesh_size, extend, cubic_box)
    else:
        # check the supplied dimension
        new_box_dim = MathTools.fix_grid_size(box_dim)
        if np.any(new_box_dim != box_dim):
            print("fixed grid size {0} -> {1}".format(box_dim, new_box_dim))
            box_dim = new_box_dim

    # print everything
    print('This is the setup on {0} ...'.format(time.ctime()))
    print("working dir:\t{0}".format(wd))
    print("centered pdb:\t{0}".format(center_pdb))
    print('fixed structure:\t{0}'.format(os.path.split(pdb_path)[-1]))
    print('ligand structure:\t{0}'.format(os.path.split(ligand_path)[-1]))
    print('non zero vdw atoms of ligand:\t{0} / {1}'.format(sum(vdw_radii_list), len(vdw_radii_list)))
    print("value of extend:\t{0}".format(extend))
    print('mesh size:\t({0},{1},{2})'.format(mesh_size[0],
        mesh_size[1], mesh_size[2]))
    print('box dim:\t({0:.0f},{1:.0f},{2:.0f})'.format(box_dim[0],
        box_dim[1], box_dim[2]))
    print('rotations:\t{0}'.format(number_of_rotations))
    print("FFT sampling:\t{0}".format(not explicit_sampling))
    print('Starting APBS ...')

    '''
    Run electrostatic calculations with APBS.
    '''
    apbs = APBSWrapper("apbs", "pdb2pqr")

    '''
    If one wants to use his own pqr file, it can be supplied and the parameter
    <use_pdb2pqr> can be set to False.
    '''
    if use_pdb2pqr is True:
        pqr_path = None
    else:
        pqr_path = pdb_path.replace('.pdb', '.pqr')
#    else:
#        raise ValueError("Unkown input for 'use_pdb2pqr'!") # JN: moved to the top of the method

    dxbox_dict = apbs.get_dxbox(pdb_path=pdb_path, mesh_size=mesh_size,
                                pqr_path=pqr_path, box_dim=box_dim,
                                box_center=box_center,
                                box_type=box_type,
                                cubic_box=cubic_box,
                                temperature=Temperature, ph=ph,
                                pdb2pqr_argv=pdb2pqr_argv)

    espbox = dxbox_dict['esp']
    vdwbox = dxbox_dict['vdw']

    # dxbox_dict has done its job
    del dxbox_dict

    print("Read apbs calculations!")

    '''
    Flood the van der waals box, because it can sometimes happen that there
    are cavities inside the protein and we do not want to fit our ligand to
    these positions as they are not physical accessible.
    '''
    vdwbox.flood()

    '''
    Create a list of independent rotations using a golden section algorithm.
    '''
    angle_stack = MathTools.get_euler_angles_for_equal_distributed_rotation(
        number_of_rotations)

    '''
    Set the interior of the protein to a negative value, so that overlap
    between the ligand and the protein is forbidden.
    '''
    vdwbox.prepare_for_geometric_matching(interior=-15)

    '''
    <interaction_energy> contains the total energy
    <counter_matrix> counts how many rotations for each grid point where
        possible without producing any overlapp
    '''
    interaction_energy = np.zeros(espbox.box.shape)
    counter_matrix = np.zeros(espbox.box.shape)

    '''
    Store the FFT of the VDW surface and the FFT of the ESP (electrostatic
    potential) in two matrices shape_scoring and esp_scoring. They will be
    correlated to the FFT of the ligand VDW surface and ESP.
    '''
    shape_scoring = FFT_correlation_scoring(vdwbox.box)
    esp_scoring = FFT_correlation_scoring(espbox.box)

    '''
    If one is interested in the actual positions of the highest scoring
    ligand orientations they can be written to disk.
    '''
    if write_top_hits:
        # object which stores the best scoring positions
        shape_results = FFT_Result()
        esp_results = FFT_Result()
        combined_result = FFT_Result()

    print("Calculating energies...")
    progress_bar = Progress_bar(len(angle_stack))

    for angle_element in angle_stack:
        # get angles
        phi = angle_element[0]
        theta = angle_element[1]
        psi = angle_element[2]

        # clone the ligand
        current_ligand_clone = PQRFile(ligand_path)
        # rotate it
        current_ligand_clone.translate_origin_and_rotate(phi,
                                                         theta,
                                                         psi)

        if explicit_sampling == True:
            # do the energy calculation in a cython script, this is much faster!
            new_interaction_energy, counter_matrix = CalculateInteractionEnergyCython.speed_up_brute_force(
                    current_ligand_clone,
                    vdwbox.box,
                    espbox,
                    counter_matrix)

            interaction_energy += np.exp(-new_interaction_energy)

        else:
            # move the ligand to the center 0,0,0
            current_ligand_clone.translate(-
                current_ligand_clone.determine_geometric_center())
            # snap its shape to the grid
            pqr_vdw = current_ligand_clone.snap_vdw_to_box(vdwbox.box_mesh_size,
                                                           vdwbox.box_dim,
                                                           vdwbox.box_offset)
            # snap its charges to the grid
            pqr_esp = current_ligand_clone.snap_esp_to_dxbox(espbox)

            # calculate the fft correlation
            shape_correlation = shape_scoring.get_correlation(pqr_vdw)
            # apply fft shift
            shape_correlation = shape_scoring.shift_fft(shape_correlation)
            # also for the electrostatic potential
            # unit is kbT/e * e -> kbT
            esp_correlation = esp_scoring.get_correlation(pqr_esp)
            # apply fft shift
            esp_correlation = esp_scoring.shift_fft(esp_correlation)

            # do the energy calculation in a cython script, this is much faster!
            # the method yields exp(-energy)
            new_interaction_energy, counter_matrix = CalculateInteractionEnergyCython.count_rotations(
                shape_correlation, esp_correlation, counter_matrix,
                shape_correlation_cutoff)

            interaction_energy += np.exp(-new_interaction_energy)

            '''
            If the best scoring ligand orientations should be written to disk,
            the data has to be prepared in a special way.
            '''
            if write_top_hits:
                # store the best scoring positions
                shape_results.find_scores(shape_correlation, phi, theta, psi,
                                      num_of_best_scores, vdwbox)
                # omit positions, that are not accessible
                esp_correlation[np.nonzero(shape_correlation < shape_correlation_cutoff)] = 0
                # find the lowest energies, '-' because 'find_scores' looks
                # for the highest energy
                esp_results.find_scores(-esp_correlation, phi, theta, psi,
                                        num_of_best_scores, espbox)

                shape_correlation[np.nonzero(shape_correlation < 0)] = 0
                # set positive, i.e. reppeling, esp to 0
                esp_correlation[np.nonzero(esp_correlation > 0)] = 0
                # normalize them and add them together
                shape_correlation = shape_correlation / np.amax(shape_correlation)
                esp_correlation = esp_correlation / np.amin(esp_correlation)
                combined_distance = np.sqrt(shape_correlation ** 2
                    + esp_correlation ** 2)
                combined_result.find_scores(combined_distance, phi, theta, psi,
                                            num_of_best_scores, espbox)

#        else:
#            raise ValueError("Unkown input for 'explicit_sampling' = {0}".format(explicit_sampling)) # JN: moved to the top of the method

        progress_bar.add()

    '''
    Use an infinite reference in solution, for which the probability of
    all rotations is given as

        p(\vec{r}_{solv}) = 1/Z_K Sum_i_N ( exp(-E_i/kb/T) ) = N / Z_K

    The ratio of the probabilty around the protein and the one in solution
    yields a equilibrium constant K

        K = p(\vec{r}) / p(\vec{r}_{solv})

    This enable us to calculate the difference in Gibbs free energy

    \Delta G(solv -> Prot) = -k_B T ln(K)
    '''
    interaction_energy = interaction_energy / float(number_of_rotations)
    zero_indices = np.nonzero(interaction_energy == 0)
    interaction_energy[zero_indices] = np.e
    interaction_energy = - np.log(interaction_energy)
    interaction_energy[np.nonzero(counter_matrix == 0)] = 0.
    interaction_energy[zero_indices] = 0.

    if write_top_hits:
        # write docking results to disk
        shape_results.write_to_disk('shape_docking.txt')
        esp_results.write_to_disk('esp_docking.txt')
        combined_result.write_to_disk('combined_docking.txt')
        print('Writing docking results!')
        # create docked structures
        ligand_pdb_path = ligand_path.replace('.pqr', '.pdb')
        shape_results.make_pdb_results(ligand_pdb_path, 'pdb_pool',
                                       'shape', 10)
        esp_results.make_pdb_results(ligand_pdb_path, 'pdb_pool',
                                       'esp', 10)
        combined_result.make_pdb_results(ligand_pdb_path, 'pdb_pool',
                                         'combined', 10)

    print('Writing grids!')
    '''
    Pymol crashes for values exceeding x < -1e37 or x > 1e37, so these values
    have to be set manually to the maximal allowed value
    '''
    interaction_energy[np.nonzero(interaction_energy > pymol_threshold)] = pymol_threshold
    interaction_energy[np.nonzero(interaction_energy < - pymol_threshold)] = - pymol_threshold
    # check for nan
    interaction_energy[np.nonzero(np.isnan(interaction_energy))] = 0

    '''
    Set energies in the grid to 0 for which no ligand orientation was possible.
    '''
    interaction_energy[np.nonzero(counter_matrix == 0)] = 0

    '''
    Write results to disk:
    '''
    print("writing energy ...")
    energy_box = DXBox(interaction_energy,
                       espbox.box_mesh_size,
                       espbox.box_offset)
    energy_box.write(energy_path)

    print("writing counter matrix ...")
    # save counter_matrix
    counter_box = DXBox(counter_matrix,
                        espbox.box_mesh_size,
                        espbox.box_offset)
    counter_box.write(counter_matrix_path)

    '''
    Analyze the protein.
    '''
    print("analyzing the data ...")
    result_data = calc_vol_and_surface(pdb_path,1.,zipped=zipped,
                                       energy_box=energy_box,
                                       counter_box=counter_box,
                                       raise_error=False)

    end_time = time.time()
    total_time = end_time - start_time
    print("total time elapsed: {0} min".format(round( total_time / 60., 2)))
    return interaction_energy
