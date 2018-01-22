# Mi 17. Mai 18:38:12 CEST 2017

import os
import shutil
import tarfile
import tempfile
import unittest
import numpy as np
from epitopsy import EnergyGrid as EG
from epitopsy.DXFile import read_dxfile
from epitopsy.tools.MathTools import RotationProtocolFibonacci

class TestCases(unittest.TestCase):
    def setUp(self):
        # remember current directory
        self.cwd = os.getcwd()
        # create an empty directory
        self.dirname = tempfile.mkdtemp(prefix='epitopsy-unittest-EG-',
                                        suffix='')
        self.ref = os.path.join(self.dirname, 'ref')
        # untar reference files
        reference_grids = 'sample_data/EnergyGrid/reference_grids.tar.bz2'
        with tarfile.open(reference_grids, 'r:bz2') as archive:
            archive.extract('ligand1.pqr', self.dirname)
            archive.extract('ligand2.pqr', self.dirname)
            archive.extract('protein1.pdb', self.dirname)
            archive.extract('protein2.pdb', self.dirname)
            archive.extract('protein1_esp.dx', self.ref)
            archive.extract('protein1_vdw.dx', self.ref)
            archive.extract('protein1_epi.dx', self.ref)
            archive.extract('protein1_mic.dx', self.ref)
            archive.extract('result_data.txt', self.ref)
            archive.extract('merge_epi.dx', self.ref)
            archive.extract('merge_mic.dx', self.ref)
    
    def tearDown(self):
        # delete empty directory
        shutil.rmtree(self.dirname)
        # move back to original directory
        os.chdir(self.cwd)
    
    def test_electrostatics_scan(self):
        '''
        Grid values should be identical up to five digits after the decimal
        dot. The sixth digit is prone to rounding errors.
        '''
        os.chdir(self.dirname)
        
        elec_kwargs = {'mesh_size': [1, 1, 1], 'verbose': False,
                       'box_center': [0, 0, 0], 'centering': None,
                       'box_type': ['esp', 'vdw']}
        scan_kwargs = {'verbose': False, 'interior': -15, 'protein_conc': 0.1,
                       'zipped': False,
                       'rotation_protocol': RotationProtocolFibonacci(90)}
        
        # run PDB2PQR and APBS
        EG.electrostatics('protein1.pdb', ['ligand1.pqr'], **elec_kwargs)
        
        # APBS Van der Waals volume (integer)
        vdw_ref = read_dxfile(os.path.join(self.ref, 'protein1_vdw.dx'), 'vdw')
        vdw = read_dxfile('protein1_vdw.dx', 'vdw')
        self.assertTrue(np.allclose(vdw.box, vdw_ref.box, rtol=0, atol=1e-5),
                        'Could not reproduce the APBS Van der Waals volume')
        
        # APBS potential (float)
        esp_ref = read_dxfile(os.path.join(self.ref, 'protein1_esp.dx'), 'esp')
        esp = read_dxfile('protein1_esp.dx', 'esp')
        self.assertTrue(np.allclose(esp.box, esp_ref.box, rtol=1e-5, atol=0),
                        'Could not reproduce the APBS potential')
        
        # scan the protein
        EG.scan('protein1.pdb', 'ligand1.pqr', APBS_dx_path='.', **scan_kwargs)
        
        # Epitopsy microstates (integer)
        mic_ref = read_dxfile(os.path.join(self.ref, 'protein1_mic.dx'), 'vdw')
        mic = read_dxfile('protein1_mic.dx', 'vdw')
        self.assertTrue(np.allclose(mic.box, mic_ref.box, rtol=0, atol=1e-5),
                        'Could not reproduce the Epitopsy microstates')
        
        # Epitopsy energy (float)
        epi_ref = read_dxfile(os.path.join(self.ref, 'protein1_epi.dx'), 'esp')
        epi = read_dxfile('protein1_epi.dx', 'esp')
        self.assertTrue(np.allclose(epi.box, epi_ref.box, rtol=1e-5, atol=0),
                        'Could not reproduce the Epitopsy energy')
        
        # Epitopsy energy statistics
        stats_ref = open(os.path.join(self.ref, 'result_data.txt')).read()
        stats_ref = dict([x.split(':') for x in stats_ref.strip().split('\n')])
        stats = open('result_data.txt').read()
        stats = dict([x.split(':') for x in stats.strip().split('\n')])
        for k in stats_ref.keys():
            if k.endswith('_score') or k in ('DG', 'conc', 'energy cutoff'):
                self.assertAlmostEqual(float(stats[k]), float(stats_ref[k]),
                    'Could not reproduce statistic "{}"'.format(k))
            else:
                self.assertEqual(int(stats[k]), int(stats_ref[k]),
                    'Could not reproduce statistic "{}"'.format(k))
        
        # multiconformational: APBS grids and scans take 1 s each -> total 6 s
        protein_paths = ['protein1.pdb', 'protein2.pdb']
        ligand_paths  = [ 'ligand1.pqr',  'ligand2.pqr']
        EG.scan_conformers(protein_paths, ligand_paths, merge=True,
                           elec_kwargs=elec_kwargs, scan_kwargs=scan_kwargs,
                           superpose_proteins=False,
                           protein_weights=[0.6, 0.4],
                           ligand_weights=[1, 1.5])
        
        # Epitopsy microstates (integer)
        mic_ref = read_dxfile(os.path.join(self.ref, 'merge_mic.dx'), 'vdw')
        mic = read_dxfile('merge_mic.dx.gz', 'vdw')
        self.assertTrue(np.allclose(mic.box, mic_ref.box, rtol=0, atol=1e-5),
           'Could not reproduce the multiconformational Epitopsy microstates')
        
        # Epitopsy energy (float)
        epi_ref = read_dxfile(os.path.join(self.ref, 'merge_epi.dx'), 'esp')
        epi = read_dxfile('merge_epi.dx', 'esp')
        self.assertTrue(np.allclose(epi.box, epi_ref.box, rtol=1e-5, atol=0),
           'Could not reproduce the multiconformational Epitopsy energy')
        
        os.chdir(self.cwd)


if __name__ == '__main__':
    unittest.main()

