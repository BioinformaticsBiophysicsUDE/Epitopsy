'''
Created on 14.11.2011

@author: chris
'''

import os
import sys
import operator
import numpy as np
from epitopsy.DXFile import DXReader, DXBox, VDWBox
from epitopsy.Structure import PDBFile


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
    :param raise_error: raise error when the box is too small if ``True``
    :type  raise_error: bool
    :param conc: protein concentration
    :type  conc: float
    :returns: Values written to disk
    :rettype: dict
    :raises ValueError: if the algorithm touched the borders of the box
    :raises ValueError: if **energy_cutoff** > 0
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
    max_concentration = 1. / vol_protein / 6.02214129e23
    if max_concentration < protein_concentration:
        print('Warning: given the protein size, its maximal concentration in '
              'solution is {:.4f} mol/L, which is below the concentration of '
              '{:.4f} mol/L chosen for the calculation.'.format(
              max_concentration, protein_concentration))
    
    # deduce the volume of solvent
    vol_outside_LAS = vol_solvation_sphere - vol_LAS - vol_protein
    vol_outside_LAS_count = vol_outside_LAS / vol_grid_point
    
    if vol_outside_LAS < 0:
        raise ValueError('Protein too large, volume outside LAS is negative.')
    
    return vol_outside_LAS_count


def estimate_DG(energy_box, counter_box, protein_concentration=0.001, Temp=310.):
    '''
    Compute the approximate binding free energy in kJ/mol and dissociation
    constant, based on the energies found on the ligand accessible surface.
    The value itself is meaningless, only a :math:`\\Delta\\Delta\\text{G}`
    between two proteins or two ligands is meaningful.

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


class FFT_Result(object):
    '''
    This class stores the results from FFT Correlation Docking.
    It stores them in a list, which contains a dictionary for every entry.
    
    The format looks like this:
    
    x_coord y_coord z_coord phi theta psi score.
    '''


    def __init__(self):
        '''
        some stuff
        '''
        self.results = []
    
    def _add_score(self, x, y, z, phi, theta, psi, score):
        '''
        This method adds a score to the list.
        '''
        new_element = {}
        new_element['x'] = x
        new_element['y'] = y
        new_element['z'] = z
        new_element['phi'] = phi
        new_element['theta'] = theta
        new_element['psi'] = psi
        new_element['score'] = score
        self.results.append(new_element)
        
    def _sort(self):
        '''
        Sort the results by their score, with the highest beeing the first.
        '''
        self.results.sort(key = operator.itemgetter('score'), reverse = True)
    
    def find_scores(self, box, phi, theta, psi, number_of_elements,
                    dxbox = None):
        '''
        This method finds the highest scoring elements from the box array. It
        needs the dxbox to transform the coordinates to real space, if none is
        given, it will store the box indices.
        '''
        # make a copy of the box
        mod_box = box.copy()
        
        for i in range(0, number_of_elements):
            # find score
            score, [i, j, k] = self._return_highest_score(mod_box)
            # transform indices [x, y, z] to real space
            if dxbox is not None:
                # transform the fft-shift
                [xdim, ydim, zdim] = dxbox.box_dim
                if i < xdim / 2:
                    x = i
                else:
                    x = i - xdim
                if j < ydim / 2:
                    y = j
                else:
                    y = j - ydim
                if k < zdim / 2:
                    z = k
                else:
                    z = k - zdim
                x = xdim / 2 - x - 1
                y = ydim / 2 - y - 1
                z = zdim / 2 - z - 1
                indices = [x, y, z]
                indices = dxbox.transform_box_to_real_space(indices) 
            
            # add score
            self._add_score(x = indices[0], y = indices[1], z = indices[2],
                            phi = phi, theta = theta, psi = psi, score = score)
            
        
        
    def _return_highest_score(self, mod_box):
        '''
        This method returns the highest score + the box coordinates of this 
        score. For further calls, it sets the highest score to the lowest, so 
        that the next highest score can be returned the next time.
        '''
        # find highest score
        highest_score = np.amax(mod_box)
        # find coordinates
        coord = np.nonzero(mod_box == highest_score)
        # save coordinates
        box_indices = [coord[0][0], coord[1][0], coord[2][0]]
        # set found score to the lowest score
        mod_box[box_indices[0], box_indices[1], box_indices[2]] = np.amin(mod_box)
        return highest_score, box_indices
    
    def write_to_disk(self, filename):
        '''
        Write results to disk. Before doing so, the list will be sorted by the
        score.
        '''
        # Make sure the results are sorted.
        self._sort()
        
        with open(filename, 'w') as f:
            for item in self.results:
                w_line = ('{0} {1} {2} {3} {4} {5} {6}\n'
                          .format(item['x'], item['y'], item['z'],
                                  item['phi'], item['theta'], item['psi'],
                                  item['score']))
                f.write(w_line)
        
    def read_from_filename(self, filename):
        '''
        Read the results from disk.
        '''
        with open(filename, 'r') as f:
            for line in f:
                if line.rstrip() != '':
                    content = line.rstrip().split(' ')
                    x = float(content[0])
                    y = float(content[1])
                    z = float(content[2])
                    phi = float(content[3])
                    theta = float(content[4])
                    psi = float(content[5])
                    score = float(content[6])
                    self._add_score(x = x, y = y, z = z,
                                    phi = phi, theta = theta, psi = psi,
                                    score = score)
    
    def make_pdb_results(self, pdb_path, storage_dir, new_name, best_x):
        '''
        This method takes the 'best_x' results and moves the pdb structure 
        to these positions. It stores the moved and rotated structures in the
        specified directory. The name of the pdb works like this:
        'new_name' + '_number.pdb', with number indexing the results from 
        0 as the best to the end.
        
        This method assuemes, that the given coordinates are in real space!!!
        '''
        # check if there are enough results
        if best_x > len(self.results):
            print('Not enough results found! Writing all results ...')
            best_x = len(self.results)
            
        # sort the list
        self._sort()
        
        # check if storage dir exists
        if not os.path.exists(storage_dir):
            os.mkdir(storage_dir)
        
        for i in range(0, best_x):
            x = self.results[i]['x']
            y = self.results[i]['y']
            z = self.results[i]['z']
            phi = self.results[i]['phi']
            theta = self.results[i]['theta']
            psi = self.results[i]['psi']
            # load pdb, its easier to load it every time again, than to create
            # copys of the pdb structure
            pdb = PDBFile(pdb_path)
            pdb.translate_origin_and_rotate(phi, theta, psi)
            pdb.translate([x, y, z])
            # make filename
            new_filename = '{0}_{1}.pdb'.format(new_name, i)
            # update filename
            new_filename = os.path.join(storage_dir, new_filename)
            pdb.save_to_file(new_filename)



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


