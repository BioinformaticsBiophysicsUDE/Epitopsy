import os
import sys
import numpy as np
from epitopsy.DXFile import DXReader, DXBox, VDWBox
from epitopsy.tools.estDG import estimate_DG


def check_indices(indices, energy_box, raise_error, energy_cutoff):
    """
    Make sure there is no grid point with energy above the energy cutoff used
    to compute statistics.
    
    :param indices: 
    """
    error_msg = ("#------------------------------------------------\n"
                 "# The box dimensions are too close to the\n"
                 "# protein and therefore the energy values near\n"
                 "# the borders are < -{0:.1f}k_BT or > {0:.1f}k_BT!\n"
                 "#------------------------------------------------")
    error_msg = error_msg.format(energy_cutoff)
    for i in range(3):
        if 0 in indices[i] or (energy_box.shape[i] - 1) in indices[i]:
            if raise_error:
                raise ValueError(error_msg)
            else:
                print error_msg


def process_vol(energy_box, vdw_box, energy_cutoff, result_dict,
                raise_error, vdw_type):
    """
    Compute statistics on grid points above the LAS.
    
    :param energy_box: energies
    :type  energy_box: :class:`DXFile.DXBox`
    :param counter_box: microstates (number of allowed ligand rotations)
    :type  counter_box: :class:`DXFile.VDWBox`
    :param energy_cutoff: energy cutoff in units of kbT
    :type  energy_cutoff: float
    :param result_dict: dictionary to update
    :type  result_dict: dict
    :param raise_error: when the box is to small, raise error if ``True``
    :type  raise_error: bool
    :param vdw_type: LAS definition
    :type  vdw_type: str
    :returns: (*dict*) Updated dictionary
    """
    # settings
    solvent_score = 1.
    protein_score = 0.

    # set interior of the protein
    vol_box = energy_box.box * vdw_box.box

    fav_energy_index = np.nonzero(vol_box < -energy_cutoff)
    check_indices(fav_energy_index, vol_box, raise_error, energy_cutoff)
    unfav_energy_index = np.nonzero(vol_box > energy_cutoff)
    check_indices(unfav_energy_index, vol_box, raise_error, energy_cutoff)

    total_volume = len(np.nonzero(vdw_box.box == solvent_score)[0])
    fav_volume = len(fav_energy_index[0])
    fav_volume_score = np.sum(vol_box[fav_energy_index])
    unfav_volume = len(unfav_energy_index[0])
    unfav_volume_score = np.sum(vol_box[unfav_energy_index])
    neutral_volume = total_volume - fav_volume - unfav_volume

    # update dictionary
    result_dict["{0}_fav_volume".format(vdw_type)] = fav_volume
    result_dict["{0}_unfav_volume".format(vdw_type)] = unfav_volume
    result_dict["{0}_fav_volume_score".format(vdw_type)] = fav_volume_score
    result_dict["{0}_unfav_volume_score".format(vdw_type)] = unfav_volume_score
    result_dict["{0}_total_volume".format(vdw_type)] = total_volume
    result_dict["{0}_neutral_volume".format(vdw_type)] = neutral_volume
    return result_dict


def process_surface(energy_box, surface_index, energy_cutoff, result_dict,
                    vdw_type):
    """
    Compute statistics on grid points on the LAS.
    
    :param energy_box: energies
    :type  energy_box: :class:`DXFile.DXBox`
    :param surface_index: indices of grid point on the LAS
    :type  surface_index: np.array
    :param energy_cutoff: energy cutoff in units of kbT
    :type  energy_cutoff: float
    :param result_dict: dictionary to update
    :type  result_dict: dict
    :param vdw_type: LAS definition
    :type  vdw_type: str
    :returns: (*dict*) Updated dictionary
    """
    surface_energies = energy_box.box[surface_index]
    
    fav_energy_index = np.nonzero(surface_energies < -energy_cutoff)
    unfav_energy_index = np.nonzero(surface_energies > energy_cutoff)
    
    total_surface = len(surface_index[0])
    fav_surface = len(fav_energy_index[0])
    unfav_surface = len(unfav_energy_index[0])
    fav_surface_score = np.sum(surface_energies[fav_energy_index])
    unfav_surface_score = np.sum(surface_energies[unfav_energy_index])
    neutral_surface = total_surface - fav_surface - unfav_surface
    
    # update dictionary
    result_dict["{0}_fav_surface".format(vdw_type)] = fav_surface
    result_dict["{0}_unfav_surface".format(vdw_type)] = unfav_surface
    result_dict["{0}_fav_surface_score".format(vdw_type)] = fav_surface_score
    result_dict["{0}_unfav_surface_score".format(vdw_type)] = unfav_surface_score
    result_dict["{0}_total_surface".format(vdw_type)] = total_surface
    result_dict["{0}_neutral_surface".format(vdw_type)] = neutral_surface
    return result_dict


def calc_vol_and_surface(energy_cutoff, energy_box, counter_box,
                         result_path="result_data.txt",
                         raise_error=True, conc=1.):
    """
    Compute statistics on the energy box.
    
    :param energy_cutoff: energy cutoff in units of kbT
    :type  energy_cutoff: float
    :param energy_box: energies
    :type  energy_box: :class:`DXFile.DXBox`
    :param counter_box: microstates (number of allowed ligand rotations)
    :type  counter_box: :class:`DXFile.VDWBox`
    :param result_path: path to output file
    :type  result_path: str
    :param raise_error: raise error if ``True`` when the box is too small
    :type  raise_error: bool
    :param conc: protein concentration
    :type  conc: float
    :returns: (*dict*) Values written to disk.
    """
    # rangecheck
    if energy_cutoff <= 0:
        raise ValueError("Argument energy_cutoff cannot be <= 0")
    if not isinstance(energy_box, DXBox):
        if isinstance(energy_box, basestring):
            energy_box = DXReader().parse(energy_box, "esp")
        else:
            raise ValueError("Argument energy_box cannot be read")
    if not isinstance(counter_box, DXBox):
        if isinstance(counter_box, basestring):
            counter_box = DXReader().parse(counter_box, "esp")
        else:
            raise ValueError("Argument counter_box cannot be read")
    
    solvent_score = 1.
    protein_score = 0.
    
    result_dict = {}
    
    # estimate Delta G
    DG = estimate_DG(energy_box, counter_box, protein_concentration = conc)[0]
    
    ## normal surface definition
    vdw_type = "normal"
    vdw_box = VDWBox(np.zeros(counter_box.box.shape),
                     counter_box.box_mesh_size,
                     counter_box.box_offset)
    
    # box setup: 1's on LAS and above LAS (= ligand was allowed to rotate
    # at least once), 0's below LAS (= no rotation allowed)
    vdw_box.box += solvent_score
    vdw_box.box[np.nonzero(counter_box.box == 0)] = protein_score
    vdw_box.flood()
    
    # calc volume
    result_dict = process_vol(energy_box, vdw_box, energy_cutoff,
                              result_dict, raise_error, vdw_type)
    
    # find surface
    vdw_box.find_solvent_surface()
    surface_index = np.nonzero(vdw_box.box == vdw_box.score_of_surface)
    
    # calc surface
    result_dict = process_surface(energy_box, surface_index, energy_cutoff,
                                  result_dict, vdw_type)
    
    ## LAS surface definition
    vdw_type = "LAS"
    vdw_box = VDWBox(np.zeros(counter_box.box.shape),
                     counter_box.box_mesh_size,
                     counter_box.box_offset)
    
    # box setup: 1's on extended LAS and above extended LAS (= no rejected
    # rotation), 0's below extended LAS (= at least 1 rejected rotation)
    vdw_box.box += solvent_score
    vdw_box.box[np.nonzero(counter_box.box < np.max(counter_box.box))] = protein_score
    vdw_box.flood()
    
    # calc volume
    result_dict = process_vol(energy_box, vdw_box, energy_cutoff,
                              result_dict, raise_error, vdw_type)
    
    # find surface
    vdw_box.find_solvent_surface()
    surface_index = np.nonzero(vdw_box.box == vdw_box.score_of_surface)
    
    # calc surface
    result_dict = process_surface(energy_box, surface_index, energy_cutoff,
                                  result_dict, vdw_type)
    
    with open(result_path, 'w') as f:
        f.write("energy cutoff:{0}\n".format(energy_cutoff))
        f.write("conc:{0}\n".format(conc))
        f.write("DG:{0}\n".format(DG))
        for key in sorted(result_dict.keys()):
            f.write("{0}:{1}\n".format(key,result_dict[key]))
    
    result_dict["DG"] = DG
    result_dict["conc"] = conc
    
    return result_dict


if __name__ == "__main__":
    usage = "python calc_vol_and_surface.py <pdb_path> <energy_cutoff> <conc>"
    try:
        pdb_path = sys.argv[1]
        energy_cutoff = float(sys.argv[2])
        conc = float(sys.argv[3])
    except:
        print(usage)
        sys.exit()
    
    if not pdb_path.endswith(".pdb"):
        print("Given pdb path argument does not end with '.pdb'!")
        sys.exit()
    
    energy_box_path = pdb_path.replace(".pdb", "_epi.dx")
    counter_box_path = pdb_path.replace(".pdb", "_mic.dx")
    if not os.path.isfile(energy_box_path):
        energy_box_path = pdb_path.replace(".pdb", "_epi.dx.gz")
    if not os.path.isfile(counter_box_path):
        counter_box_path = pdb_path.replace(".pdb", "_mic.dx.gz")
        
    result_dict = calc_vol_and_surface(energy_cutoff, energy_box_path,
                                       counter_box_path, conc = conc)


