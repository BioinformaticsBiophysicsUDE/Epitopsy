# Mi 17. Mai 18:38:12 CEST 2017

import os
import shutil
import tarfile
import tempfile
import unittest
import numpy as np
from epitopsy import EnergyGrid
from epitopsy.DXFile import read_dxfile

class TestCases(unittest.TestCase):
    def setUp(self):
        # remember current directory
        self.cwd = os.getcwd()
        # create an empty directory
        self.dirname = tempfile.mkdtemp(prefix='epitopsy-unittest-EnergyGrid-',
                                        suffix='')
        self.ref = os.path.join(self.dirname, 'ref')
        # untar reference files
        reference_grids = 'sample_data/EnergyGrid/reference_grids.tar.bz2'
        with tarfile.open(reference_grids, 'r:bz2') as archive:
            archive.extract('ligand.pqr', self.dirname)
            archive.extract('protein.pdb', self.dirname)
            archive.extract('protein_esp.dx', self.ref)
            archive.extract('protein_vdw.dx', self.ref)
            archive.extract('protein_epi.dx', self.ref)
            archive.extract('protein_mic.dx', self.ref)
            archive.extract('result_data.txt', self.ref)
    
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
        
        # run PDB2PQR and APBS
        EnergyGrid.electrostatics('protein.pdb', ['ligand.pqr'],
                                  mesh_size=[1,1,1], box_center=[0,0,0],
                                  center_pdb=False, use_pdb2pqr=True,
                                  box_type=['esp', 'vdw'], verbose=False)
        
        # APBS Van der Waals volume (integer)
        vdw_ref = read_dxfile(os.path.join(self.ref, 'protein_vdw.dx'), 'vdw')
        vdw = read_dxfile('protein_vdw.dx', 'vdw')
        self.assertTrue(np.allclose(vdw.box, vdw_ref.box, rtol=0, atol=1e-5),
                        'Could not reproduce the APBS Van der Waals volume')
        
        # APBS potential (float)
        esp_ref = read_dxfile(os.path.join(self.ref, 'protein_esp.dx'), 'esp')
        esp = read_dxfile('protein_esp.dx', 'esp')
        self.assertTrue(np.allclose(esp.box, esp_ref.box, rtol=1e-5, atol=0),
                        'Could not reproduce the APBS potential')
        
        # scan the protein
        EnergyGrid.scan('protein.pdb', 'ligand.pqr', number_of_rotations=90,
                        protein_conc=0.1, interior=-15., APBS_dx_path='.',
                        zipped=False, verbose=False)
        
        # Epitopsy microstates (integer)
        mic_ref = read_dxfile(os.path.join(self.ref, 'protein_mic.dx'), 'vdw')
        mic = read_dxfile('protein_mic.dx', 'vdw')
        self.assertTrue(np.allclose(mic.box, mic_ref.box, rtol=0, atol=1e-5),
                        'Could not reproduce the Epitopsy microstates')
        
        # Epitopsy energy (float)
        epi_ref = read_dxfile(os.path.join(self.ref, 'protein_epi.dx'), 'esp')
        epi = read_dxfile('protein_epi.dx', 'esp')
        self.assertTrue(np.allclose(epi.box, epi_ref.box, rtol=1e-5, atol=0),
                        'Could not reproduce the Epitopsy energy')
        
        # Epitopsy statistics
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
        
        os.chdir(self.cwd)


if __name__ == '__main__':
    unittest.main()

