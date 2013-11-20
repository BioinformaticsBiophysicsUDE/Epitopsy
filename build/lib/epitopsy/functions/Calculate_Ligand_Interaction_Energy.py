'''
Created on 24.11.2011

@author: chris
'''

import sys
import numpy as np

from epitopsy.Structure import PQRFile, PDBFile
from epitopsy.APBS import APBSWrapper, InFile
from epitopsy.DXFile import DXBox, DXReader
from epitopsy.scoring.FFTCorrelationScoring import FFT_correlation_scoring
from epitopsy.result.FFT_Result import FFT_Result


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
    apbs = APBSWrapper("apbs", "pdb2pqr")
    template_in = InFile("", "", "", mesh_size)
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
    
if __name__ == '__main__':
    usage = ("Usage:\npython Calculate_Ligang_Interaction_Energy.py " + 
              "--protein=<pdb_path> --ligand=<pqr_path")
    arguments = {}
    for item in sys.argv[1:]:
        arguments[item.split('=')[0]] = item.split('=')[1]
    
    if '--protein' in arguments:
        pdb_path = arguments['--protein']
    else:
        print(usage)
        sys.exit(1)
    
    if '--ligand' in arguments:
        ligand_path = arguments['--ligand']
    else:
        print(usage)
        sys.exit(1)
    
    if '--mesh_size' in arguments:
        m = float(arguments['--mesh_size'])
        mesh_size = [m, m, m]
    else:
        mesh_size = [0.5, 0.5, 0.5]
    
    Calculate_Ligand_Interaction_Energy(pdb_path, ligand_path, mesh_size)
