'''
Created on Feb 13, 2012

@author: chris
'''
import os
import sys
import time
import numpy as np

from epitopsy.Optimization.Optimizer import GA_Optimizer, Individual, Population
from epitopsy.Optimization.Selector import nsga_random_tournament, \
    nsga_tournament

from base_classes import CMObject
from cmmanagers import CMTask
from server_processes import CMServerProcess
from multiprocessing import Process
#from test.test_commons import CMPTestWithManager
from cmlogging import QueueHandler

from Queue import Empty

#===============================================================================
# HOW TO USE THIS OPTIMIZER:
#    1. supply a config file
#    2. the config is supposed to contain more parameters, please add them to 
#        'the data_dict'
#    3. modify the 'score_pop' class method of CMOptimizer
#    4. modify the 'mutate_pop' class method of CMOptimizer
#    5. supply a 'run_evaluator.py' script which will be copied to all clients
#        and calculates the fitness for a given individual. The script gets the
#        the following input:
#            unique_id, pdb_content, gene_data, data_dict
#===============================================================================

class CMOptimizer(GA_Optimizer, CMObject):
    def __init__(self, config, ip, logging_queue, job_queue, server, auth_key, **kwdargs):
        '''
        Properties, which have to be specified in the config-file:
        
        pop_size -> # individuals
        generations -> # generations
        mutate_pop -> supply a function, which copies the data of each
                        individual, mutates it and creates a new individual
                        with the new data and returns a new population
                        afterwards
        score_pop -> supply a function, which scores a given population
        select_pop -> supply a function, which selects a new pool from the
                        three given arguments: parents, children, pareto_front
                        a further argument is 'self.minimization' because the 
                        selection is sensitive to minimization or maximization
        analyze_pop -> supply a function, which analyzes a given population 
                        from the four given arguments: parents, children,
                        pareto_front and current_generation. The result should
                        be a dictionary e.g.
                        {'gen': current_gen, 'hyper_vol_parents' : x,
                         hyper_vol_pareto':y}
        minimization -> Is this a minimization or a maximization?
        init_data -> list of data that is used for the initialization
        crossover -> Can be True, False or a function. The function works on 
                        the parent population, therefore it is necessary to 
                        create new individuals, because the parents should not 
                        be modified. After the crossover this population will 
                        be mutated. If no function is supplied and crossover is
                        True, the implemented function will be used
        crossover_frac -> Fraction of the parent population, which will be
                            mutated with the built in crossover
        data_dict -> additional data, which can be used in any given function,
                    a dictionary seems to be a good idea.
        run_id -> can be used to create unique directory names for the same
                    generation and population sizes
        log_data -> log data or not
        print_generation -> either print the current generation or not
        '''
        CMObject.__init__(self, config, ip, queues = {})
        
        #self.logger.removeHandler(logging.StreamHandler)
        #self.logger.addHandler(QueueHandler(logging_queue))
        self.job_queue = job_queue
        self.server = server
        self.data_queue = self.server.queues['data']
        
        # dictionary that can contain all necesary data
        data_dict = {}
        #data_dict['job_queue'] = self.job_queue
        #data_dict['data_queue'] = self.data_queue
        
        #=======================================================================
        # PUT CODE HERE THAT READS IN THE INITIAL DATA AND STORES IT AS A LIST
        #=======================================================================
        init_data = []
        with open(config['init_data'], 'r') as f:
            raise NotImplementedError('The initial data needs to be read!')
            
        
        
        pop_size = int(config['pop_size'])
        generations = int(config['generations'])
        
        if config['minimize'] == 'True':
            minimize_problem = True
        else:
            minimize_problem = False
        
        if config['crossover'] == 'True':
            use_crossover = True
            crossover_frac = float(config['crossover_frac'])
        else:
            use_crossover = False
            crossover_frac = None
        
        if config['log_data'] == 'True':
            log_data = True
        else:
            log_data = False
            
        if config['print_generation'] == 'True':
            print_generation = True
        else:
            print_generation = False
        
        # optional, i.e. can be set to 'None'
        pdb_pool_dir = config['pdb_pool_dir']
        if pdb_pool_dir == 'None':
            data_dict['pdb_pool_dir'] = None
        else:
            data_dict['pdb_pool_dir'] = pdb_pool_dir

        GA_Optimizer.__init__(self, population_size = pop_size,
                              generations = generations,
                              mutate_pop = None,
                              score_pop = None,
                              select_pop = nsga_tournament,
                              analyze_pop = None,
                              minimization = minimize_problem,
                              init_data = init_data,
                              crossover = use_crossover,
                              crossover_frac = crossover_frac,
                              data_dict = data_dict,
                              run_id = 0,
                              log_data = log_data,
                              print_generation = print_generation)
        
    def score_pop(self, unscored_pop, data_dict):
        '''
        This method is used to score a population. Here the distributed computing
        is used! The actual scoring is performed by the supplied scoring script,
        which is run by the clients. For the proper function the client returns a
        list with the scores
        '''
        # this is a dictionary which maps the unique_id to the individuals instance
        individual_dict = {}
        
        scored_pop = Population()
        for indi in unscored_pop.individuals:
            # add individual to the dictionary
            individual_dict[indi.unique_id] = indi
            
            #=======================================================================
            # PUT CODE HERE TO SEND NECCESARY DATA TO CLIENTS
            #=======================================================================
            # is the pdb needed?
            pdb_content = None
            self.job_queue.put(CMTask([indi.unique_id, pdb_content, indi.data,
                                       data_dict]))
            
        
        while len(scored_pop.individuals) < len(unscored_pop.individuals):
            try:
                time.sleep(1)
                # the results should contain the unique_id, a list of scores
                # and a variable that contains a new pdb from the simulation,
                # otherwise (i.e. if no pdbs are used) None
                result_data = self.data_queue.get(True, 60)
            except Empty:
                continue
            
            unique_id = result_data.data[0]
            scores = result_data.data[1]
            pdb_content = result_data.data[2]
            
            scored_indi = individual_dict[unique_id]
            scored_indi.score = scores
            
            scored_pop.add(scored_indi)
            
            # update the data queue
            self.data_queue.task_done(result_data)
            
            if pdb_content is not None:
                #===================================================================
                # PUT CODE HERE TO WRITE THE NEW PDB CONTENT TO A FILE
                #===================================================================
                pass
        
        return scored_pop
    
    def mutate_pop(self, population, current_generation, data_dict):
        '''
        This function is used to mutate a given population.
        '''
        pass
    
    def secure_run(self):
        
        time.sleep(5)
        print('--------------------------------------------------')
        print('--------------------------------------------------')
        start = time.time()
        self.start_ga()
        end = time.time()
        print('--------------------------------------------------')
        print('Time needed for optimization: {0:.2f} min'.format((end - start) / 60))
        print('--------------------------------------------------')
        print('--------------------------------------------------')
        
        self.server.feeder_feeded()
        
        return
        
            
class CMOptimizerProcess(Process, CMOptimizer):
    def __init__(self, config, ip, logging_queue, job_queue, server, auth_key, **kwdargs):
        super(CMOptimizerProcess, self).__init__(target = self.real_run)
        CMOptimizer.__init__(self, config, ip, logging_queue, job_queue, server, auth_key, **kwdargs)


if __name__ == '__main__':
    config_file = sys.argv[1]
    myserver = CMServerProcess(config_file, terminate_childs = True)
    myserver.start()
    myserver.join()
