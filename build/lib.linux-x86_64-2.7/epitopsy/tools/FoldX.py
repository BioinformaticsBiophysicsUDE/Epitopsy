'''
Created on Aug 31, 2011

@author: Christoph Wilms
'''

import csv
import os

from subprocess import call

class FoldX(object):
    '''
    This class implements a method to run a FoldX simulation.
    '''


    def __init__(self, pdb_name, sequence_parent = None, sequence_child = None,
                 repair_flag = True, new_name = None, clean_up = False,
                 Temp = 310, pH = 7):
        '''
        * pdb_name refers to a structure that should be mutated and analyzed#
            without '.pdb'!
        * sequence_parent is the sequence of the pdb (pdb_name)
        * sequence_child is a mutated version of sequence_parent
        * repair_flag indicates, if the pdb should be repaired at first
        * new_name is the name of the scored structure
        * Temp is the temperature of the experiment in K
        * pH is the pH of the experiment
        '''
        # repair or don't
        self.repair_flag = repair_flag
        
        # all names miss the *.pdb at the end
        if pdb_name.endswith('.pdb'):
            self.pdb_name = pdb_name.replace('.pdb', '')
        else:
            self.pdb_name = pdb_name
        
        # name of a repaired pdb
        self.repaired_pdb_name = "RepairPDB_{0}".format(self.pdb_name)
        
        # name for build model
        if self.repair_flag is True:
            self.build_model_pdb = self.repaired_pdb_name
        else:
            self.build_model_pdb = self.pdb_name
        
        # name after build model
        if new_name is not None:
            if new_name.endswith('.pdb'):
                self.new_pdb_name = new_name.replace('.pdb', '')
            else:
                self.new_pdb_name = new_name
        else:
            self.new_pdb_name = '{0}_1'.format(self.build_model_pdb)
        
        self.sequence_parent = sequence_parent
        self.sequence_child = sequence_child
        
        self.clean_up = clean_up
        
        self.Temp = Temp
        self.pH = pH
        
        # filenames
        # list.txt
        self.pdb_list = 'list_{0}.txt'.format(self.pdb_name)
        # mutant_file.txt
        self.mutant_file = 'mutant_file_{0}.txt'.format(self.pdb_name)
        # run-repair.txt
        self.run_repair_file = 'run_repair_{0}.txt'.format(self.pdb_name)
        # run-buildModel
        self.run_build_file = 'run_buildModel_{0}.txt'.format(self.pdb_name)
        
        
        # FoldX path
        self.foldx_path = './FoldX_3_5_1'
        
        #ddG
        # stable < less stable < totally unstable
        self.ddG = None
        
    def _make_FoldX_run_Repair(self):
        with open(self.run_repair_file, 'w') as outdata:
            # first part
            outdata.write('<TITLE>FOLDX_runscript;\n<JOBSTART>#;\n<PDBS>#;\n<BATCH>{0};\n<COMMANDS>FOLDX_commandfile;\n<RepairPDB>#;\n<END>#;\n<OPTIONS>FOLDX_optionfile;\n<Temperature>'.format(self.pdb_list))
            # write Temperature
            outdata.write(str(self.Temp) + ';\n')
            # write more stuff
            outdata.write('<R>#;\n<pH>')
            # write PH
            outdata.write(str(self.pH) + ';\n')
            # write the rest
            outdata.write('<IonStrength>0.050;\n<water>-CRYSTAL;\n<metal>-CRYSTAL;\n<VdWDesign>2;\n<OutPDB>true;\n<pdb_hydrogens>false;\n<END>#;\n<JOBEND>#;\n<ENDFILE>#;')
        
        
    def _make_FoldX_run_Build(self):
        with open(self.run_build_file, 'w') as outdata:
            # first part
            outdata.write('<TITLE>FOLDX_runscript;\n<JOBSTART>#;\n<PDBS>#;\n<BATCH>{0};\n<COMMANDS>FOLDX_commandfile;\n<BuildModel>#,{1};\n<END>#;\n<OPTIONS>FOLDX_optionfile;\n<Temperature>'.format(self.pdb_list, self.mutant_file))
            # write Temperature
            outdata.write(str(self.Temp) + ';\n')
            # write more stuff
            outdata.write('<R>#;\n<pH>')
            # write PH
            outdata.write(str(self.pH) + ';\n')
            # write the rest
            outdata.write('<IonStrength>0.050;\n<water>-CRYSTAL;\n<metal>-CRYSTAL;\n<VdWDesign>2;\n<OutPDB>true;\n<pdb_hydrogens>false;\n<complex_with_DNA>false;\n<END>#;\n<JOBEND>#;\n<ENDFILE>#;')
        
    
    # repair!
    def _make_FoldX_mut_List(self):
        with open(self.mutant_file, 'w') as individualList:
            individualList.write(self.sequence_parent + '\n')
            individualList.write(self.sequence_child)
        
    
    def _make_FoldX_pdb_List(self, pdb_name):
        with open(self.pdb_list, 'w') as pdblist:
            pdblist.write("{0}.pdb".format(pdb_name))
    
    def _extract_FoldX_ddG(self):
        ddGfilename = 'Dif_BuildModel_{0}.fxout'.format(self.build_model_pdb)
        ddGcontent = csv.reader(open(ddGfilename), delimiter = '\t')
        rowcount = 1
        for row in ddGcontent:
            if rowcount == 10:
                ddG = row[1]
            rowcount += 1
        self.ddG = float(ddG)
    
    def extract_FoldX_ddG(self):
        return self.ddG
    
    # repair 
    def run_FoldX(self):
        # check if rotabase exists
        if not os.path.exists('rotabase.txt'):
            raise AttributeError('Missing rotabase.txt in {0}!'.format(os.getcwd()))
        else:
            # trash goes nowhere
            foldX_trash_file = '/dev/null'
            
            # repair PDB
            if self.repair_flag is True:
                self._make_FoldX_pdb_List(self.pdb_name)
                self._make_FoldX_run_Repair()
                f = open(foldX_trash_file, 'w')
                call([self.foldx_path, '-runfile', self.run_repair_file], stdout = f)
                f.close()
                
            
            
            self._make_FoldX_pdb_List(self.build_model_pdb)
            self._make_FoldX_mut_List()
            self._make_FoldX_run_Build()
            f = open(foldX_trash_file, 'w')
            call([self.foldx_path, '-runfile', self.run_build_file], stdout = f)
            f.close()
            # remove FoldX header ... 
            f = open('{0}_1.pdb'.format(self.build_model_pdb), 'r')
            pdb_data = f.readlines()
            f.close()
            with open('{0}.pdb'.format(self.new_pdb_name), 'w') as f:
                for line in pdb_data:
                    if (line[0:4] == 'ATOM' or line[0:3] == 'TER'
                        or line[0:6] == 'HETATM' or line[0:3] == 'END'):
                        f.write(line)
            
            # get ddG
            self._extract_FoldX_ddG()
            
            # clean up
            if self.clean_up is True:
                os.remove(self.pdb_list)
                os.remove(self.mutant_file)
                os.remove(self.run_build_file)
                if self.repair_flag is True:
                    # names
                    repaired_name = self.repaired_pdb_name
                    # remove
                    os.remove(self.run_repair_file)
                    os.remove('RepairPDB_{0}.fxout'.format(self.pdb_name))
                    os.remove('{0}.pdb'.format(repaired_name))
                
                os.remove('BuildModel_{0}.fxout'.format(self.build_model_pdb))
                os.remove('PdbList_BuildModel_{0}.fxout'.format(self.build_model_pdb))
                os.remove('Dif_BuildModel_{0}.fxout'.format(self.build_model_pdb))
                os.remove('Average_BuildModel_{0}.fxout'.format(self.build_model_pdb))
                os.remove('Raw_BuildModel_{0}.fxout'.format(self.build_model_pdb))
                os.remove('WT_{0}_1.pdb'.format(self.build_model_pdb))
                os.remove('{0}_1.pdb'.format(self.build_model_pdb))
