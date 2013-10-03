'''
Created on Feb 16, 2012

@author: chris
'''

import os
import numpy as np

from epitopsy.Structure import PDBFile



def extract_pdbs_from_pdb_models(pdb_model_ensemble_path, new_dir,
                                 name_template, ref_pdb_path = None):
    '''
    This function reads all models from 'pdb_model_ensemble_path' and writes
    each model as a new pdb to the supplied directory. The new name of the 
    pdbs is: 'given_name' + '_x.pdb' where x is the model id.
    If a reference pdb is supplied, the function will fit the 'CA' atoms.
    '''
    if not os.path.exists(new_dir):
        os.mkdir(new_dir)

    with open(pdb_model_ensemble_path) as f:
        pdb_model_content = f.readlines()

    model_counter = -1
    for line in pdb_model_content:
        if line.startswith('MODEL'):
            line = line.strip()
            current_model = int(line[5:])
            model_counter = current_model
            f = open(os.path.join(new_dir, '{0}_{1}.pdb'.format(name_template,
                                                                current_model)), 'w')
            continue
        elif line.startswith('ENDMDL'):
            f.close()
            continue
        elif(line.startswith('ATOM') or line.startswith('HETATM')
                or line.startswith('TER')):
            f.write(line)

    # fit pdbs onto reference
    if ref_pdb_path is not None:
        ref_pdb = PDBFile(ref_pdb_path)
        pdb_ensemble = os.listdir(new_dir)
        for file_name in pdb_ensemble:
            new_pdb = PDBFile(os.path.join(new_dir, file_name))
            ref_pdb.superimpose_given_pdb_onto_self(new_pdb, 'CA')
            new_pdb.save_to_file(new_pdb.structure_path)


def get_one_atom_pdb(pdb_path, atom_coord):
    '''
    This function returns a pdb structure, which consists of only one atom.
    
    This is kind of a hack :-)
    '''
    atom_id = 1
    atom_type = 'CA'
    res_type = 'ALA'
    chain = 'A'
    res_id = 1
    with open(pdb_path, 'w') as f:
        f.write('ATOM      {0} {1}   {2} {3}   {4}       0.000   0.000   0.000  1.00  1.00             \n'.format(atom_id, atom_type, res_type, chain, res_id))
        f.write('TER')
        
    
    pdb = PDBFile(pdb_path)
    pdb.structure[chain][res_id][atom_type].set_coord(np.array(atom_coord))
    pdb.save_to_file(pdb.structure_path)
    return pdb
