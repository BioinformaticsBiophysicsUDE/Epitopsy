import os
import sys

import numpy as np

from epitopsy.DXFile import DXReader, VDWBox
from epitopsy.tools.style import style
from epitopsy.tools.estDG import calcDG

def check_indices(indices, box, raise_error):
    for i in range(3):
        if 0 in indices[i] or (box.shape[i]-1) in indices[i]:
            if raise_error:
                raise ValueError("The algorithm touched the borders of the box!")
            else:
                print("""#------------------------------------------------
# The box dimensions are too close to the
# protein and therefore the energy values near
# the borders are < -1k_BT or > 1k_BT!
#------------------------------------------------""")

def process_vol(energy_box, vdw_box, energy_cutoff, result_dict,
                raise_error, vdw_type):
    '''
    This function is not well designed.

    Args:
        probability_box -> DXBox contains the energies.
        vdw_box -> flooded VDWBox.
        energy_cutoff -> count grid points with an energy < limit
        result_dict -> dictionary

    Returns:
        Updated dictionary
    '''
    # settings
    solvent_score = 1.
    protein_score = 0.

    # set interior of the protein
    vol_box = energy_box.box * vdw_box.box

    fav_energy_index = np.nonzero(vol_box < -energy_cutoff)
    check_indices(fav_energy_index, vol_box, raise_error)
    unfav_energy_index = np.nonzero(vol_box > energy_cutoff)
    check_indices(unfav_energy_index, vol_box, raise_error)

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


def calc_vol_and_surface(pdb_path, energy_cutoff,
                         zipped=True,
                         result_path="result_data.txt",
                         energy_box=None,
                         counter_box=None,
                         raise_error=True,
                         conc=1.):
    '''
    Args:
        pdb_path -> path to the pdb file
        energy_cutoff -> positive float
        ----------------------------------- optional:
        zipped -> are the matrices gzipped or not? Default is True.
        result_path -> path to which the results are to be written.

    Returns:
        A dictionary with the keys:
            'total_volume'
            'fav_volume'
            'unfav_volume'
            'neutral_volume'
            'total_surface'
            'fav_surface'
            'unfav_surface'
            'neutral_surface'
    '''
    if energy_cutoff <= 0:
        raise ValueError("Given energy value is <= 0")
    energy_cutoff = np.abs(energy_cutoff)

    if energy_box is None and counter_box is None:
        energy_path = style["energy_box_file_path"]
        counter_path = style["counter_box_file_path"]
        if zipped:
            if not energy_path.endswith(".gz"):
                energy_path += ".gz"
        else:
            energy_path.rstrip(".gz")

        energy_box = DXReader().parse(energy_path, "esp")
        counter_box = DXReader().parse(counter_path, "esp")

    solvent_score = 1.
    protein_score = 0.

    result_dict = {}

    # estimate DG
    [DG,k_D] = calcDG(energy_box=energy_box,
        counter_box=counter_box,
        conc=conc)
    result_dict["DG"] = DG
#    result_dict["k_D"] = k_D
    result_dict["conc"] = conc

    ## normal surface definition
    vdw_type = "normal"
    vdw_box = VDWBox(np.zeros(energy_box.box.shape)+solvent_score,
                     energy_box.box_mesh_size,
                     energy_box.box_offset)

    # set normal
    vdw_box.box[np.nonzero(counter_box.box == 0)] = protein_score
    vdw_box.flood()

    # calc volume
    result_dict = process_vol(energy_box,
                              vdw_box,
                              energy_cutoff,
                              result_dict,
                              raise_error,
                              vdw_type)

    # find surface
    vdw_box.find_solvent_surface()
    surface_index = np.nonzero(vdw_box.box == vdw_box.score_of_surface)

    # calc surface
    result_dict = process_surface(energy_box,
                                  surface_index,
                                  energy_cutoff,
                                  result_dict,
                                  vdw_type)

    ## LAS surface definition
    vdw_type = "LAS"
    vdw_box = VDWBox(np.zeros(energy_box.box.shape)+solvent_score,
                     energy_box.box_mesh_size,
                     energy_box.box_offset)

    # set LAS
    vdw_box.box[np.nonzero(counter_box.box < np.max(counter_box.box))] = protein_score
    vdw_box.flood()

    # calc volume
    result_dict = process_vol(energy_box,
                              vdw_box,
                              energy_cutoff,
                              result_dict,
                              raise_error,
                              vdw_type)

    # find surface
    vdw_box.find_solvent_surface()
    surface_index = np.nonzero(vdw_box.box == vdw_box.score_of_surface)

    # calc surface
    result_dict = process_surface(energy_box,
                                  surface_index,
                                  energy_cutoff,
                                  result_dict,
                                  vdw_type)


    with open(result_path, 'w') as f:
        f.write("# energy cutoff:{0}\n".format(energy_cutoff))
        for key,value in result_dict.items():
            f.write("{0}:{1}\n".format(key,value))

    return result_dict

if __name__ == '__main__':
    usage = """python calc_vol_and_surface.py <pdb_path> <energy_cutoff> <conc>"""
    if len(sys.argv) != 4:
        print(usage)
        sys.exit()

    ## settings
    pdb_path = sys.argv[1]
    energy_cutoff = float(sys.argv[2])
    conc = float(sys.argv[3])
    if not pdb_path.endswith(".pdb"):
        print("Given pdb path argument does not end with '.pdb'!")
        sys.exit()

    result_dict = calc_vol_and_surface(pdb_path, energy_cutoff,conc=conc)
