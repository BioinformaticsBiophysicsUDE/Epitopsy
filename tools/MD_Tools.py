'''
Created on Feb 16, 2012

@author: chris
'''

import os
import sys
import shutil
import numpy as np
import subprocess

from epitopsy.Structure import PDBFile
from epitopsy.tools.PDB_Tools import extract_pdbs_from_pdb_models

def dump_gromacs_snapshots(trajectory_file, tpr_file, delta_t, new_pdb_path):
    '''
    This function dumps a pdb snapshot at the given time interval. The data
    will be stored in the new pdb path.
    '''
    process = subprocess.Popen(['trjconv_d', '-f', trajectory_file, '-s',
                                tpr_file, '-o', new_pdb_path, '-pbc',
                                'nojump', '-dt', str(delta_t)],
                               stdin = subprocess.PIPE)
    process.communicate('Protein')


def calculate_rmsd_from_pdbmodel(pdb_model_path, pdb_frame_dir, ref_pdb_path, atom_types, delta_t):
    extract_pdbs_from_pdb_models(pdb_model_path, pdb_frame_dir, 'dump')
    
    # calculate the rmsd
    time_list = []
    rmsd_list = []

    ref_pdb = PDBFile(ref_pdb_path)
    frame_ensemble = os.listdir(pdb_frame_dir)

    # sort the ensemble by time
    frame_ensemble.sort(key=lambda x:int(x.split('_')[-1].replace('.pdb','')))

    for pdb_item in frame_ensemble:
        pdb = PDBFile(os.path.join(pdb_frame_dir, pdb_item))
        rmsd = ref_pdb.get_rmsd_rotation_translation_from_superposition(pdb, atom_types)['rmsd']
        rmsd_list.append(rmsd)
        snap_time = (int(pdb_item.split('_')[-1].replace('.pdb','')) - 1) * delta_t
        time_list.append(snap_time)

    return time_list, rmsd_list


def calculate_rmsf_from_pdbmodel(pdb_model_path, pdb_frame_dir, fit_atom_types = None,
                                 rmsf_type = None):
    # extract each model from the MODEL for rmsf fitting
    extract_pdbs_from_pdb_models(pdb_model_path, pdb_frame_dir, 'dump')
    
    # calculate the rmsf
    coord_list = []

    frame_ensemble = os.listdir(pdb_frame_dir)
    # sort the ensemble by time
    frame_ensemble.sort(key=lambda x:int(x.split('_')[-1].replace('.pdb','')))
    
    # 
    ref_pdb = PDBFile(os.path.join(pdb_frame_dir, frame_ensemble[0]))
    
    # sanity check -> get all ids to compare them in the for loop
    if rmsf_type is None:
        # nothing to check, every protein has a protein, except None
        check_id_list = 1
    elif rmsf_type == 'CA' or rmsf_type == 'res':
        # record all CA residue ids
        check_id_list = ref_pdb.get_residue_id_list()
    else:
        raise AttributeError("Unkown input for the rmsf calculation '{0}'!".format(rmsf_type))
    
    coord_list = []
    for i,pdb_item in enumerate(frame_ensemble):
        pdb = PDBFile(os.path.join(pdb_frame_dir, pdb_item))
        pdb.superimpose_self_onto_given_pdb(ref_pdb, fit_atom_types)
        # sanity check part 2
        if i == 0:
            if rmsf_type is None:
                # again nothing to check
                check_id_list_2 = 1
            elif rmsf_type == 'CA' or rmsf_type == 'res':
                check_id_list_2 = pdb.get_residue_id_list()
            if check_id_list != check_id_list_2:
                raise ValueError("There is an error with the index lists!")

        if rmsf_type is None:
            # calculate geometric center
            pdb_center = pdb.determine_geometric_center()
            coord_list.append([pdb_center])
        elif rmsf_type == 'CA':
            # calculate the rmsf for each c alpha atom
            # -> get all c alpha atoms
            ca_atom_list = pdb.get_atoms_of_type(rmsf_type)
            ca_coord_list = []
            for atom in ca_atom_list:
                ca_coord_list.append(atom.get_coord())
            coord_list.append(ca_coord_list)
        elif rmsf_type == 'res':
            def calculate_center_of_mass(res):
                cms = []
                for atom in res:
                    cms.append(atom.get_coord())
                return np.mean(cms,0)
            # calculate the rmsf for each residue
            # -> get residues
            # -> calculate their geometrical center
            res_list = pdb.get_residues()
            res_coord_list = []
            for res in res_list:
                res_coord_list.append(calculate_center_of_mass(res))
            coord_list.append(res_coord_list)
            
    # array with dimensions
    # 0: number of snapshots
    # 1: number of particles (e.g. 'CA', 'res', etc.)
    # 2: 3d coordinates
    coord_list = np.array(coord_list)
    coord_mean = []
    for i in range(coord_list.shape[1]):
        coord_mean.append(np.mean(coord_list[:,i],0))
    coord_mean = np.array(coord_mean)
    
    rmsf_list = []
    # iterate over each particle
    for i in range(coord_list.shape[1]):
        dif_list = coord_list[:,i] - coord_mean[i]
        norm_dif_list = np.sqrt(dif_list[:,0]**2
                                + dif_list[:,1]**2
                                + dif_list[:,2]**2)
        rmsf = np.sqrt(np.mean(norm_dif_list**2))
        rmsf_list.append(rmsf)

    # this is returned as well, so that one can map the rmsf onto the
    # structure
    id_list = check_id_list
    return id_list, rmsf_list

def calculate_rmsd_from_trajectory(traj_file, tpr_file, delta_t, ref_pdb_path,
                                   atom_types = None):
    '''
    This function returns two lists, which contain the time and the
    rmsd to the given reference structure path. 
    The function extracts snapshots from the trajectory file and stores them 
    in /dev/shm.
    '''
    try:
        dir_path = os.path.join(os.sep, 'dev', 'shm', 'traj_dump')
        counter = 0
        # find a unique name in /dev/shm
        while os.path.exists(dir_path):
            if counter >= 10:
                raise ValueError("Something is wrong because the path '{0}' already exists!".format(dir_path))
            dir_path = '{0}{1}'.format(dir_path, np.random.randint(10))
            counter += 1
        
        os.mkdir(dir_path)
        
        # extract the frames from the trajectory
        pdb_model_path = os.path.join(dir_path, 'frames.pdb')
        dump_gromacs_snapshots(traj_file, tpr_file, delta_t, pdb_model_path)
        
        # extract each model from the MODEL for rmsd fitting
        pdb_frame_dir = os.path.join(dir_path, 'frame_dump')
        
        time_list, rmsd_list = calculate_rmsd_from_pdbmodel(pdb_model_path,
                                                            pdb_frame_dir,
                                                            ref_pdb_path,
                                                            atom_types,
                                                            delta_t)

        # remove the directory
        shutil.rmtree(dir_path)

    except (Exception, KeyboardInterrupt), e:
        import traceback
        print(traceback.format_exc())
        # remove the directory
        shutil.rmtree(dir_path)
        raise e

    return time_list, rmsd_list


def calculate_rmsf_from_trajectory(traj_file, tpr_file, delta_t,
                                   fit_atom_types = None,
                                   rmsf_type = None):
    '''
    This function returns two lists, one with a list of indieces (e.g. res 
    numbers) and one with the calculated rmsfs.
    The unit of the rmsf is Angstroem.

    Arguments:
        - fit_atom_types -> atom types that are used to fit the trajectory 
          onto the reference structure
        - rmsf_type -> 'res', 'CA', None; types that are used for the
          calculation of the rmsf. 'res' means rmsf per residue, 'CA' 
          calculates the rmsf for each c alpha atom and None for the
          geometric center of the protein.
    '''
    try:
        dir_path = os.path.join(os.sep, 'dev', 'shm', 'traj_dump')
        counter = 0
        # find a unique name in /dev/shm
        while os.path.exists(dir_path):
            if counter >= 10:
                raise ValueError("Something is wrong because the path '{0}' already exists!".format(dir_path))
            dir_path = '{0}{1}'.format(dir_path, np.random.randint(10))
            counter += 1
        
        os.mkdir(dir_path)

        
        # extract the frames from the trajectory
        pdb_model_path = os.path.join(dir_path, 'frames.pdb')
        dump_gromacs_snapshots(traj_file, tpr_file, delta_t, pdb_model_path)
        
        # extract each model from the MODEL for rmsf fitting
        pdb_frame_dir = os.path.join(dir_path, 'frame_dump')
        
        id_list, rmsf_list = calculate_rmsf_from_pdbmodel(pdb_model_path,
                                                          pdb_frame_dir,
                                                          fit_atom_types,
                                                          rmsf_type)
        # remove the directory
        shutil.rmtree(dir_path)
    
    except (Exception, KeyboardInterrupt), e:
        import traceback
        print(traceback.format_exc())
        # remove the directory
        shutil.rmtree(dir_path)
        raise e

    return id_list, rmsf_list

