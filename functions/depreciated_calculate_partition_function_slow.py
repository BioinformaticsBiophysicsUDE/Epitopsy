"""
Created on March 7, 2013

@author: Christoph Wilms
"""
import os
import numpy as np

from epitopsy.Structure import PDBFile, PQRFile
from epitopsy.APBS import APBSWrapper
from epitopsy.DXFile import DXBox
from epitopsy.cython import CalculateInteractionEnergyCython
from epitopsy.tools import MathTools
from epitopsy.tools.UtilityClasses import Progress_bar
from epitopsy.tools.style import style


def calculate_interaction_energy(pdb_path, ligand_path, mesh_size,
        number_of_rotations=150,
        Temperature=310.,
        extend=None,
        ph=None,
        use_pdb2pqr=True,
        cubic_box=True,
        center_pdb=False,
        box_center=[0,0,0],
        box_dim=None,
        box_type=["esp","vdw"],
        write_top_hits=False):
    '''
    This function calculates the ensemble averaged energy in units kbT.

    partition function:
    Z_k = Sum_i_N exp(-E_i/kb/T)

    probability of state i:
    p_i = 1/Z_k * exp(-E_i/kb/T)

    averaged ensemble energy:
    <E> = Sum_i_N p_i * E_i
        = Sum_i_N 1./Z_k * E_i * exp(-E_i/kb/T)

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
        center_pdb -> if True the pdb will be centered at box_center, default is
                    False
        box_center -> center of the apbs calculations, default is [0,0,0]
        box_dim -> dimension of the box, default is None and it is calculated
                    automatically
        box_type -> types of boxes APBS should write to disk. We need at least
                    ["esp","vdw"], "smol" is optional -> ["esp","vdw","smol"]
        write_top_hits -> write the top scoring conformations

    Returns:
        None.
    '''
    ## settings
    energy_path = style["energy_box_file_path"]
    LAS_energy_path = style["LAS_box_file_path"]
    counter_matrix_path = style["counter_box_file_path"]
    np.seterr('raise')
    pymol_threshold = 0.9e37

    '''
    Read the pdb file. If box dimensions have been supplied they will be used,
    but one has to make sure, that the box is big enough to contain the
    protein. If none box dimensions have been supplied, they will be
    calculated from the protein.
    '''

    pdb_struct = PDBFile(pdb_path)
    if box_dim is None:
        # calculate the box dimensions from the pdb structure
        box_dim = pdb_struct.get_dxbox_dim(mesh_size, extend, cubic_box)
    else:
        # check the supplied dimension
        new_box_dim = MathTools.fix_grid_size(box_dim)
        if new_box_dim != box_dim:
            print("fixed grid size {0} -> {1}".format(box_dim, new_box_dim))
            box_dim = new_box_dim

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
        pdb_struct.translate(new_coord-pdb_coord)
        pdb_struct.save_to_file(pdb_struct.structure_path)

    ligand_struct = PQRFile(ligand_path)
    ligand_struct_pdb = PDBFile(ligand_path.replace('.pqr','.pdb'))

    '''
    The parameter <extend> can be used to extend the box in each dimension.
    If extend is none, we will use the maximal diameter of the ligand to
    extend the box, because otherwise not all rotation will fit in the box.
    '''
    if extend is None:
        extend = int(ligand_struct.determine_max_diameter() + 2)

    # print everything
    print('This is the setup ...')
    print("centered pdb:\t{0}".format(center_pdb))
    print('fixed structure:\t{0}'.format(os.path.split(pdb_path)[-1]))
    print('ligand structure:\t{0}'.format(os.path.split(ligand_path)[-1]))
    print("value of extend:\t{0}".format(extend))
    print('mesh size:\t({0},{1},{2})'.format(mesh_size[0], mesh_size[1], mesh_size[2]))
    print('box dim:\t({0:.0f},{1:.0f},{2:.0f})'.format(box_dim[0], box_dim[1], box_dim[2]))
    print('rotations:\t{0}'.format(number_of_rotations))
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
    elif use_pdb2pqr is False:
        pqr_path = pdb_path.replace('.pdb', '.pqr')
    else:
        raise ValueError("Unkown input for 'use_pdb2pqr'!")

    dxbox_dict = apbs.get_dxbox(pdb_path = pdb_path, mesh_size = mesh_size,
                                pqr_path = pqr_path, box_dim=box_dim,
                                box_center=box_center,
                                box_type = box_type,
                                cubic_box = cubic_box,
                                temperature = Temperature, ph = ph)


    espbox = dxbox_dict['esp']
    vdwbox = dxbox_dict['vdw']

    # read fixed pqr
    if pqr_path is None:
        pqr_path = pdb_path.replace('.pdb', '.pqr')

    pqr_struct = PQRFile(pqr_path)

    # dxbox_dict has done its job
    del dxbox_dict

    print("Read apbs calculations!")

    '''
    Flood the van der waals box, because it can sometimes happen that there
    are cavities inside the protein and we do not want to fit our ligand to
    these positions as they are not physical accssible.
    '''
    vdwbox.flood()

    '''
    Create a list of independet rotations using a golden section algorithm.
    '''
    angle_stack = MathTools.get_euler_angles_for_equal_distributed_rotation(number_of_rotations)

    '''
    <interaction_energy> contains the total energy
    <counter_matrix> counts how many rotations for each grid point where
        possible without producing any overlapp
    '''
    interaction_energy = np.zeros(espbox.box.shape)
    counter_matrix = np.zeros(espbox.box.shape)

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

        # do the energy calculation in a cython script, this is much faster!
        new_interaction_energy, counter_matrix = CalculateInteractionEnergyCython.speed_up_brute_force(
                current_ligand_clone,
                vdwbox.box,
                espbox,
                counter_matrix)

        interaction_energy += np.exp(-new_interaction_energy)

        '''
        If the best scoring ligand orientations should be written to disk,
        the data has to be prepared in a special way.
        '''

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
    zero_indices = np.nonzero(interaction_energy==0)
    interaction_energy[zero_indices] = np.e
    interaction_energy = - np.log(interaction_energy)
    interaction_energy[np.nonzero(counter_matrix==0)] = 0.
    interaction_energy[zero_indices] = 0.

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
    DXBox(interaction_energy,
          espbox.box_mesh_size,
          espbox.box_offset).write(energy_path)

    # save LAS surface
    interaction_energy[np.nonzero(counter_matrix < number_of_rotations)] = 0
    DXBox(interaction_energy,
          espbox.box_mesh_size,
          espbox.box_offset).write(LAS_energy_path)

    # save counter_matrix
    DXBox(counter_matrix,
          espbox.box_mesh_size,
          espbox.box_offset).write(counter_matrix_path)

    return interaction_energy


