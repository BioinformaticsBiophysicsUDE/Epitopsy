import numpy as np
from epitopsy.DXFile import VDWBox


def find_LAS(counter_box):
    """
    Detect the Ligand Accessible Surface (LAS), defined as the center of a
    molecular probe rolling on the protein van der Waals surface. The resulting
    DXBox encodes the solvent as 1, the LAS as 2 the protein interior as 0.
    
    :param counter_box: microstates (number of allowed ligand rotations)
    :type  counter_box: :class:`DXFile.VDWBox`
    :returns: (:class:`DXFile.VDWBox`) LAS protein volume
    """
    # create a VDW box with 1's in the solvent (= ligand was allowed to rotate)
    # and 0's inside the protein (= no rotation was allowed)
    vdw = VDWBox(np.zeros(counter_box.box.shape),
                 counter_box.box_mesh_size,
                 counter_box.box_offset)
    solv_pos = np.nonzero(counter_box.box > 0)
    vdw.box[solv_pos] = 1
    
    # encode the LAS with 2's (= ligand was allowed to rotate)
    vdw.flood()
    vdw.find_solvent_surface()
    
    return vdw


def calc_volume_solvent(protein_concentration, LAS_points_count,
                        protein_points_count, mesh_size):
    """
    Compute the volume of solution containing exactly one protein, subtract
    the volume of the protein to get the volume of solvent, divide by the mesh
    size to get the theoretical number of grid points containing the solvent.
    
    :param protein_concentration: protein concentration in mol/L
    :type  protein_concentration: float
    :param LAS_points_count: number of grid points on the LAS
    :type  LAS_points_count: int
    :param protein_points_count: number of grid points below the LAS
    :type  protein_points_count: int
    :param mesh_size: mesh size in angstroms
    :type  mesh_size: list
    """
    # given a concentration c, find the volume of solvent containing 1 molecule
    vol_solvation_sphere = 1. / (6.02214129e23 * protein_concentration)
    
    # find the LAS volume and protein volume (much smaller than vol_solvation)
    Angstrom_cube_to_liter = 1000 * (1e-10)**3  # volume of 1 A^3 in L
    vol_grid_point = np.prod(mesh_size) * Angstrom_cube_to_liter
    vol_LAS     = vol_grid_point * LAS_points_count
    vol_protein = vol_grid_point * protein_points_count
    
    # deduce the volume of solvent
    vol_outside_LAS = vol_solvation_sphere - vol_LAS - vol_protein
    vol_outside_LAS_count = vol_outside_LAS / vol_grid_point
    
    return vol_outside_LAS_count


def estimate_DG(energy_box, counter_box, protein_concentration=1., Temp=310.):
    '''
    Compute the approximate binding free energy in kJ/mol and dissociation
    constant, based on the energies found on the ligand accessible surface.

    :param energy_box: energies
    :type  energy_box: :class:`DXFile.DXBox`
    :param counter_box: microstates (number of allowed ligand rotations)
    :type  counter_box: :class:`DXFile.VDWBox`
    :param Temp: temperature
    :param Temp: float
    :param protein_concentration: protein concentration (mol/L)
    :param protein_concentration: float
    :returns: (*tuple*) Binding free energy in kJ/mol and dissociation constant
    '''
    # molar gas constant in J/mol/K
    R = 8.3144598
    
    # number of grid points contained on the LAS and below the LAS
    LAS = find_LAS(counter_box)
    protein_points_count = np.sum(LAS.box == 0)
    LAS_points_count     = np.sum(LAS.box == 2)
    LAS_points = np.nonzero(LAS.box == 2)
    
    # compute total number of solvent grid points if the DXBox dimensions were
    # extended to reach a concentration of *conc*
    solvent_points_count = calc_volume_solvent(protein_concentration,
              LAS_points_count, protein_points_count, energy_box.box_mesh_size)
    
    # DG = - R T ln( K )
    # K =  [E_LAS] / [nonLAS]
    # [LAS] = sum( exp( -E_{onLAS} / (k_B T) ) )
    # [nonLAS] = sum( exp( -E_{notonLAS} / (k_B T) ) )
    # E_{notonLAS} = 0 # in water
    LAS_dG = energy_box.box[LAS_points]
    LAS_K_sum = np.sum(np.exp(-LAS_dG))
    return (-R*Temp*np.log(LAS_K_sum / solvent_points_count) / 1000, LAS_K_sum)


