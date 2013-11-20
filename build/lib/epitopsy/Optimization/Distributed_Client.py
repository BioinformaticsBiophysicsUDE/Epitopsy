'''
Created on Feb 13, 2012

@author: chris
'''

import os
import shutil
import time

from base_classes import *
from cmmanagers import CMTask
from test.test_commons import CMPTestWithManager

import run_evaluation

class CMOptimizer_Client(CMParticipant):
    '''
    This class will be run on the clients. 
    
    Usage:
        1) always sync a script 'run_evaluation.py' with the client manager
        2) this script should contain a function:
            'evaluate(path_to_dir, seq_data, pdb_path)'
           which returns a list with the calculated scores. If no pdb is used,
           pdb_path is set to 'None'. The function also has to make sure, that
           the modified pdb is also written to the pdd_path, because it will
           be read after the execution. Please switch the working directory
           always to the new directory, this way the source folder will not 
           be filled up with trash!
    '''
    
    def __init__(self, config, ip, basic = False, max_jobs = 50, *args, **kwargs):
        CMParticipant.__init__(self, config, ip, basic = False, max_jobs = 50, *args, **kwargs)
    
    def secure_run(self):
        '''
        Method that keeps the clients going.
        '''
        self.logger.info('CLIENT {0}: client is up'.format(self.ip))
        
        #--- Endless running method waiting to process jobs listed in job queue
        while not self.shutting_down.is_set():
            
            #--- get details from job queue
            time.sleep(2)
            '''
            0. individual.unique_id
            1. pdb_content
            2. individual.data
            '''
            job_details = self.get_task()
            if job_details is None:
                self.logger.info('CLIENT {0}: Job queue ist leer'.format(self.ip))
                continue                
            
            # assign parameters
            unique_id = job_details.data[0]
            pdb_content = job_details.data[1]
            gene_data = job_details.data[2]
            data_dict = job_details.data[3]
            
            self.logger.info('CLIENT {0}: Got job from job queue for individual {1}'.format(self.ip, unique_id))
            
            
            #--- Create jobs working directory
            try:
                work_dir = os.path.join(os.getcwd(), unique_id)
                os.mkdir(work_dir)
            except Exception as e:
                self.logger.error("CLIENT {0}: An error was raised while creating working directory for individual '{1}' with '{2}'".format(self.ip, unique_id, type(e)))
                self.task_failed(self.act_task, "CLIENT {0}: An error was raised while creating working directory for individual '{1}'".format(self.ip, unique_id))
                raise RuntimeError('Evil exception')
            
            #--- write pdb_to file, if it exists
            if pdb_content is not None:
                try:
                    pdb_path = os.path.join(work_dir, unique_id + '.pdb')
                    with open(pdb_path) as f:
                        f.writelines(pdb_content)
                
                except Exception as e:
                    self.logger.error("CLIENT {0}: An error was raised while writing the pdb '{1}.pdb' to disk with '{1}'!".format(self.ip, unique_id, type(e)))
                    self.task_failed(self.act_task, "CLIENT {0}: An error was raised while writing the pdb '{1}.pdb' to disk.".format(self.ip, unique_id))
            else:
                pdb_path = None
            
            #--- run evaluation and get a dictonary, that contains the 'scores'
            # and a new structure, if it is wished
            try:
                result_dict = run_evaluation.evaluate(unique_id, pdb_content,
                                                      gene_data, data_dict)
            except Exception as e:
                self.logger.error("CLIENT {0}: An error was raised while running the evaluation script!".format(self.ip, unique_id, type(e)))
                self.task_failed(self.act_task, "CLIENT {0}: An error was raised while running the evaluation script.".format(self.ip, unique_id))
            
            #--- read new pdb structure
            if pdb_content is not None:
                try:
                    new_pdb_content = result_dict['pdb_content']
                except:
                    new_pdb_content = None
            
            else:
                new_pdb_content = None
            
            #--- return results
            '''
            0. unique_id
            1. list of scores
            2. pdb_content
            '''
            self.data_queue.put(CMTask([unique_id, result_dict['scores'],
                                        new_pdb_content]))
            
            #--- remove job from job queue
            try:
                self.task_done(job_details)
                time.sleep(1)
            except KeyError:
                self.logger.error("CLIENT {0}: KeyError was raised while executing task_done method for individual '{1}'".format(self.ip, unique_id))
            
            
            #--- remove jobs working directory
            try:
                shutil.rmtree(work_dir)
            except Exception as e:
                self.logger.error("CLIENT {0}: An error was raised while removing jobs working directory '{1}' with '{2}'".format(self.ip, unique_id, type(e)))

            

class CMOptimizer_ClientProcess(Process, CMOptimizer_Client):
        def __init__(self, config, ip, basic = False, max_jobs = 500, **kwargs):
            super(CMOptimizer_ClientProcess, self).__init__(target = self.real_run)
            CMOptimizer_Client.__init__(self, config, ip, basic = False, max_jobs = 500, **kwargs)

