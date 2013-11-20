'''
Created on 17.01.2012

@author: chris
'''

import os
import sys
import time
import numpy as np
import random
import multiprocessing

from Queue import Empty

from epitopsy.cython import optimizer_cython
from epitopsy.Optimization.Analysis import hypervolume


def dominates(p, q, optimization_type):
    '''
    Does p dominate q?
    '''
    def _dominates_max(p, q):
        d = True
        for x, y in zip(p, q): 
            d = not (y > x)
        return d
    
    def _dominates_min(p, q):
        d = True
        for x, y in zip(p, q): 
            d = not (y < x)
        return d
    
    # this is only to get the same optimization type string as in population
    pop = Population()
    
    if optimization_type == pop.minimize:
        return _dominates_min(p, q)
    elif optimization_type == pop.maximize:
        return _dominates_max(p, q)
    else:
        raise ValueError("Unkown optimization type: '{0}'".format(optimization_type))
    
def is_absolutely_fitter(p, q, optimization_type):
    '''
    Is p fitter in ALL score values than q?
    '''
    def _is_bigger(p, q):
        for x, y in zip(p, q): 
            d = x >= y
            if d is  False:
                break

        return d
    
    def _is_smaller(p, q):
        for x, y in zip(p, q): 
            d = x <= y
            if d is False:
                break

        return d
    
    # this is only to get the same optimization type string as in population
    pop = Population()
    
    if optimization_type == pop.minimize:
        return _is_smaller(p, q)
    elif optimization_type == pop.maximize:
        return _is_bigger(p, q)
    else:
        raise ValueError("Unkown optimization type: '{0}'".format(optimization_type))
 

class SMS_EMOA(object):
    
    minimize = 'min'
    maximize = 'max'

    def __init__(self, pop_size, eval_steps, n_processes, init_data,
                 data_dict = None,
                 get_mutant = None,
                 score_indi = None,
                 optimization_type = 'min',
                 run_id = 0,
                 analyze_pop = None,
                 log_data = False,
                 print_eval_step = False):
        '''
        Args:
            pop_size -> size of the population
            eval_steps -> number of evaluations to perform
            n_processes -> number of process that should run in parallel
            init_data -> data for initialisation (list of genes -> list of
                lists) 
            data_dict -> dictionary that is given to each 'user-supplied'
                function like get_mutant. This allows sharing of special
                information.
            get_mutant -> function that selects an individual from a given
                population and mutates it.
                Args:
                    population
                    received_steps,
                    submitted_steps
                    data_dict
                The unique has to be set to:
                    '{0}_{0}'.format(received_individuals, submitted_individuals)
                -> otherwise it will lead to problems if one wants to restart
                a terminated optimization
            score_indi -> score a given individual, arguments are the
                individual and data_dict
            optimization_type -> min or max
            run_id -> use it to add your own flavour to the log dir path
            analyze_pop -> supply a function, which analyzes a given
                population from the five given arguments: population,
                pareto_front, eval_step, optimization_type and
                data_dict. The result should be a dictionary e.g.
                {'eval_step': eval_step, 'hyper_vol_pop' : x, hyper_vol_pareto':y}
            log_data -> True or False
            print_eval_step -> True or False

        Notice:
            Keep in mind, that this algorithm does not guarantee an ever
            increasing hypervolume, because there is at least one case, where
            the selection method leads to an decrease in hypervolume:
            # +                  x
            #  +                o
            #    +
            #      +
            #        +
            # # # # # # # # # # # # # #
            + -> pareto front
            x -> reference point
            o -> individual, which is removed because of the first criterion,
                 which uses pareto dominance
            The removal of the 'o'-individual leads to a decrease of the hypervolume.
            ---> WRONG! Hypervolume is defined by the pareto front! One
            source of a decrease in the hypervolume is a 'wrong' reference
            point. This leads to jumps in the hypervolume development.
        '''
        self.pop_size = pop_size
        self.eval_steps = eval_steps
        self.starting_iteration = 0
        self.n_processes = n_processes
        self.init_data = init_data
        self.data_dict = data_dict
        # If None is given, it is assumed, there is a derivative of this class
        # which implements this functions as class methods
        if get_mutant is not None:
            self.get_mutant = get_mutant
        # If None is given, it is assumed, there is a derivative of this class
        # which implements this functions as class methods
        if score_indi is not None:
            self.score_indi = score_indi
        self.optimization_type = optimization_type
        
        # is there a broken calculation?
        self.resume_work = False
        
        # logging stuff
        self.run_id = run_id
        self.log_data = log_data
        self.log_dir = 'log_pool_{0}_{1}_{2}_{3}'.format(self.pop_size,
                                                     self.eval_steps,
                                                     self.n_processes,
                                                     self.run_id)
        self.log_pool_path = os.path.join(self.log_dir, 'log_pool.LOG')
        self.log_pareto_path = os.path.join(self.log_dir, 'log_pareto.LOG')
        self.log_individual_path = os.path.join(self.log_dir, 'log_individual.LOG')
        self.log_individual_status = False
        
        self.analyze_pop = analyze_pop
        self.analysis_log_file_name = os.path.join(self.log_dir, 'log_analysis.LOG')

        self.print_eval_step = print_eval_step

        if self.log_data is True or self.analyze_pop is not None:
            if os.path.exists(self.log_dir):
                # maybe the calculation was interrupted?
                # -> check last pool
                if (os.path.exists(self.log_individual_path) and
                    os.path.exists(self.log_pool_path) and
                    os.path.exists(self.log_pareto_path)):
                    print('Found existing log files')
                    self._read_existing_log_files()
                    self.resume_work = True
                    
                    print('Resuming work from iteration {0}'.format(self.starting_iteration))
                else:    
                    print('Log dir already exists, but there is something wrong! Exiting ...!')
                    sys.exit(1)
            
            else:
                os.mkdir(self.log_dir)
        
        if self.log_data is True and self.resume_work is False:
            # touch the files
            with open(self.log_pool_path, 'w') as f:
                pass
                
            with open(self.log_individual_path, 'w') as f:
                pass
                
            with open(self.log_pareto_path, 'w') as f:
                pass
   
    def evolve(self):
        '''
        Start the optimization.
        '''
        if self.n_processes == 1:
            self._evolve_single()
        else:
            self._evolve_parallel()


    def _evolve_single(self):
        print('Setting up SMS-EMOA')
        start_time = time.time()
        
        # in case the calculations crash, it is nice if the opbject
        # still has accsess to the front
        self.backup_pop = Population()
        self.backup_pareto = Population()

        if self.resume_work is False:
            init_pop = Population()
            for i in range(self.pop_size):
                indi = Individual(self.init_data[i], 'x', i)
                init_pop.add(indi)
    
            # fill the queue with the initial population
            pop = Population()
            for indi in init_pop.individuals:
                indi_scored = self.score_indi(indi, self.data_dict)
                pop.add(indi_scored)
            
                # logging
                if self.log_data is True:
                    # log new individual
                    self._log_individual(indi_scored, 0)
    
    
            # determine the starting pareto_front
            pareto_front = pop.get_pareto_frontier(self.optimization_type)
            # logging
            if self.log_data is True:
                # log population
                self._log_population(pop, pareto_front, 0)
            
            # counts the number of submitted individuals
            eval_steps_submitted = 0
            # counte the number of received individuals
            eval_steps_received = 0
            
        else:
            pop = self.init_pop
            pareto_front = self.init_pareto_front
            # counts the number of submitted individuals
            eval_steps_submitted = self.starting_iteration
            # counte the number of received individuals
            eval_steps_received = self.starting_iteration
            
        
        for i in range(eval_steps_received, self.eval_steps):
            if self.print_eval_step is True:
                print('#{0}/{1}'.format(eval_steps_received, self.eval_steps))
            eval_steps_submitted += 1    
            eval_steps_received += 1
            
            new_indi = self.get_mutant(pop, eval_steps_received,
                    eval_steps_submitted,
                    self.data_dict)
            
            new_indi = self.score_indi(new_indi, self.data_dict)

            # this is a mu + 1 algorithm
            pop.add(new_indi)
            self.backup_pre_pop = pop
            self.backup_pre_pareto = pareto_front
            pop = self.select_pop(pop)

            # update pareto front
            mixed_pool = Population()
            mixed_pool.add_population(pop)
            mixed_pool.add_population(pareto_front)
            pareto_front = mixed_pool.get_pareto_frontier(self.optimization_type)
            
            # backup
            self.backup_pop = pop
            self.backup_pareto = pareto_front

            # logging
            if self.log_data is True:
                # log new individual
                self._log_individual(new_indi, eval_steps_received)
                # log populationdata
                self._log_population(pop, pareto_front, eval_steps_received)
                 
            # analyze if a function has been supplied
            if self.analyze_pop is not None:
                self._write_analysis(eval_steps_received,
                                     self.analyze_pop(pop,
                                                      pareto_front,
                                                      eval_steps_received,
                                                      self.optimization_type,
                                                      self.data_dict))

        # store last population and the pareto front
        self.final_pop = pop
        self.final_pareto_front = pareto_front
        
        print('SMS-EMOA finished after {0}s'.format(time.time() - start_time))
    
    def _evolve_parallel(self):
        print('Setting up SMS-EMOA')
        start_time = time.time()
        
        # in case the calculations crash, it is nice if the opbject
        # still has accsess to the front
        self.backup_pop = Population()
        self.backup_pareto = Population()

        # queues for communication
        tasks = multiprocessing.JoinableQueue()
        results = multiprocessing.Queue()
        
        # set up worker bees
        num_worker_bees = self.n_processes
        print('Creating {0} worker bees'.format(num_worker_bees))
        worker_bee_list = [Worker_Bee(tasks, results)
                         for i in range(num_worker_bees)]
        
        # let'em work
        for worker_bee in worker_bee_list :
            worker_bee.start()
        
        if self.resume_work is False:
            init_pop = Population()
            for i in range(self.pop_size):
                indi = Individual(self.init_data[i], 'x', i)
                init_pop.add(indi)
    
            # fill the queue with the initial population
            for indi in init_pop.individuals:
                tasks.put(Task(indi, self.score_indi, self.data_dict.copy()))
            
            tasks.join()
            
            # get the scored population
            pop = Population()
            for i in range(self.pop_size):
                indi = results.get()
                pop.add(indi)
    
                # logging
                if self.log_data is True:
                    # log new individual
                    self._log_individual(indi, 0)
    
    
            # determine the starting pareto_front
            pareto_front = pop.get_pareto_frontier(self.optimization_type)
            # logging
            if self.log_data is True:
                # log population
                self._log_population(pop, pareto_front, 0)
            
            # counts the number of submitted individuals
            eval_steps_submitted = 0
            # counte the number of received individuals
            eval_steps_received = 0
            
        else:
            pop = self.init_pop
            pareto_front = self.init_pareto_front
            # counts the number of submitted individuals
            eval_steps_submitted = self.starting_iteration
            # counte the number of received individuals
            eval_steps_received = self.starting_iteration
            
        finished = False
        while not finished:
            if eval_steps_submitted == self.starting_iteration:
                # this is the first generation and we need to fill the queue
                for i in range(num_worker_bees):
                    if eval_steps_submitted == self.eval_steps:
                        # this should not happen, but it could ... 
                        # num of workers > num of evaluation steps
                        break
                    eval_steps_submitted += 1
                    mutant_indi = self.get_mutant(pop,
                                                  eval_steps_received,
                                                  eval_steps_submitted,
                                                  self.data_dict)
                    # set mutant indi!
                    mutant_indi.update_unique_id(eval_steps_received,
                                                 eval_steps_submitted)
                    tasks.put(Task(mutant_indi, self.score_indi,
                                   self.data_dict.copy()))

            try:
                # test if the result queue is empty -> submit a new individual
                new_indi = results.get(timeout = 1)
                eval_steps_received += 1
                if self.print_eval_step is True:
                    print('#{0}/{1}'.format(eval_steps_received, self.eval_steps))
                
                # this is a mu + 1 algorithm
                pop.add(new_indi)
                pop = self.select_pop(pop)

                # submit a new individual
                if eval_steps_submitted < self.eval_steps:
                    eval_steps_submitted += 1
                    mutant_indi = self.get_mutant(pop,
                                                  eval_steps_submitted,
                                                  eval_steps_received,
                                                  self.data_dict)
                    # set mutant indi!
                    mutant_indi.update_unique_id(eval_steps_received,
                                                 eval_steps_submitted)
                    tasks.put(Task(mutant_indi, self.score_indi,
                                   self.data_dict.copy()))

                # update pareto front
                pareto_front.add_population(pop)
                pareto_front = pareto_front.get_pareto_frontier(self.optimization_type)
                
                # backup
                self.backup_pop = pop
                self.backup_pareto = pareto_front

                # logging
                if self.log_data is True:
                    # log new individual
                    self._log_individual(new_indi, eval_steps_received)
                    # log populationdata
                    self._log_population(pop, pareto_front, eval_steps_received)
                 
                # analyze if a function has been supplied
                if self.analyze_pop is not None:
                    self._write_analysis(eval_steps_received,
                                         self.analyze_pop(pop,
                                                          pareto_front,
                                                          eval_steps_received,
                                                          self.optimization_type,
                                                          self.data_dict))

            except Empty:
                # there is no individual -> wait a little bit
                time.sleep(0.1)

            if eval_steps_received == self.eval_steps:
                finished = True
        
        for i in range(num_worker_bees):
            tasks.put(None)
    
        # wait for the tasks to finish
        tasks.join()
        print('Killed {0} worker bees'.format(num_worker_bees))

        # store last population and the pareto front
        self.final_pop = pop
        self.final_pareto_front = pareto_front
        
        print('SMS-EMOA finished after {0}s'.format(time.time() - start_time))

    def select_pop(self, pop):
        '''
        This selection is based on the standard selection for a sms-emoa
        and reduces the given population by one individual.
        The sms-emoa uses two rules to select an individual:
            1. Calculate number of dominating individuals for each individual,
               select the one, which is dominated by the most
            2. If rule 1 is not clear use the hypervolume contribution of the 
               competing individuals and remove the individual, which
               contributes less to the hypervolume

        Notice: There is no guarantee for an ever increasing hypervolume!
        '''
        # set up kill id
        kill_id = None

        # check if there is a double unique_id
        unique_list = []
        for indi in pop.individuals:
            if indi.unique_id in unique_list:
                self.error_pop = pop
                print(unique_list)
                print(indi.unique_id)
                raise ValueError("Unique id '{0}' has been found twice!".format(indi.unique_id))
            unique_list.append(indi.unique_id)

        dummy_pop = Population()
        dummy_pop.add_population(pop)

        # calculate number of dominating individuals
        dominant_number_list = dummy_pop.get_dominant_number_list(self.optimization_type)
        max_num = np.amax(dominant_number_list)
        max_pos = np.nonzero(dominant_number_list == max_num)

        # There are now two cases:
        if len(max_pos[0]) == 1:
            # there is just one individual and this is the one, which will
            # be eliminated
            kill_indi = dummy_pop.individuals[max_pos[0]]
            kill_id = kill_indi.unique_id
        elif len(max_pos[0]) >= 2:
            # get reference point
            ref = dummy_pop.get_reference_point(self.optimization_type)

            # dict to store the hypervolume contribution
            kill_candidates = {}
            for pos in max_pos[0]:
                indi = np.array(dummy_pop.individuals)[pos]
                vol = hypervolume([indi.score], self.optimization_type, ref)
                # in case it sits at a corner
                if vol == 0:
                    vol = np.inf
                kill_candidates[indi.unique_id] = vol
            
            # find minimum contribution
            hyper_min = min(kill_candidates.values())
            hyper_min_pos = np.nonzero(np.array(kill_candidates.values()) == hyper_min)
            # two cases:
            if len(hyper_min_pos[0]) == 1:
                kill_pos = hyper_min_pos[0]
                kill_id = kill_candidates.keys()[kill_pos]
            elif len(hyper_min_pos[0]) > 1:
                kill_pos = random.choice(hyper_min_pos[0])
                kill_id = kill_candidates.keys()[kill_pos]
            elif np.inf in ref:
                # this algorithm crashes if one of the individuals has an infinity score
                infinity_list = []
                for indi in pop.individuals:
                    for score in indi.score:
                        if np.isinf(score):
                            infinity_list.append(indi.unique_id)
                kill_id = random.choice(infinity_list)
            else:
                for indi in pop.individuals:
                    print(indi.score)
                # I can not imagine what could go wrong ...
                # It should be related to the hypervolume
                self.error_pop = pop
                raise AttributeError('#sms_emoa_selection_error_5')
        
        else:
            # I can not imagine what could go wrong ...
            # It should be related to the dominated number list
            self.error_pop = pop
            raise AttributeError('#sms_emoa_selection_error_6')
    
        if kill_id is None:
            self.error_pop = pop
            raise AttributeError('Unkown error, error code: #sms_emoa_selection_error_3')
        else:
            dummy_pop.remove_individual(kill_id)
            return dummy_pop


    def _log_individual(self, indi, eval_step):
        '''
        This method logs a given individual.
        '''
        with open(self.log_individual_path, 'a') as f:
            if eval_step == 0 and self.log_individual_status is False:
                self.log_individual_status = True
                f.write('#************************************CONFIG************************************\n')
                f.write('#* population_size = {0}\n'.format(self.pop_size))
                f.write('#* iterations = {0}\n'.format(self.eval_steps))
                f.write('#* optimization_type = {0}\n'.format(self.optimization_type))
                f.write('#* time = {0}\n'.format(time.ctime()))
                f.write('#**************************************END*************************************\n')
                f.write('{0}\n'.format(indi.get_r_log_string()))
            f.write(indi.get_log_string(eval_step))

    def _log_population(self, population, pareto_pop, eval_step):
        self._write_population(population, eval_step, self.log_pool_path)
        self._write_population(pareto_pop, eval_step, self.log_pareto_path)

    def _write_population(self, population, eval_step, filename):
        with open(filename, 'a') as f:
            if eval_step == 0:
                number_of_genes = len(population.individuals[0].data)
                number_of_scores = len(population.individuals[0].score)
                f.write('#************************************CONFIG************************************\n')
                f.write('#* population_size = {0}\n'.format(self.pop_size))
                f.write('#* iterations = {0}\n'.format(self.eval_steps))
                f.write('#* optimization_type = {0}\n'.format(self.optimization_type))
                f.write('#* number_of_genes = {0}\n'.format(number_of_genes))
                f.write('#* number_of_scores = {0}\n'.format(number_of_scores))
                f.write('#* time = {0}\n'.format(time.ctime()))
                f.write('#**************************************END*************************************\n')
                f.write('{0}\n'.format(population.individuals[0].get_r_log_string()))
                
            f.write('#@ Iteration: {0}\n'.format(eval_step))
            for indi in population.individuals:
                f.write(indi.get_log_string(eval_step))
            f.write('#@ END\n')
    
    
    def _write_analysis(self, eval_step, analysis_dict):
        '''
        This function writes out the supplied dictionary. 
        '''
        # steady state cannot log first iteration
        if eval_step == 1:
            with open(self.analysis_log_file_name, 'w') as f:
                # keys should be strings!
                str2write = ' '.join(analysis_dict.keys())
                f.write(str2write + '\n')
        
        with open(self.analysis_log_file_name, 'a') as f:
            cache_list = []
            for value in analysis_dict.values():
                cache_list.append(str(value))
            cache_str = ' '.join(cache_list)
            f.write(cache_str + '\n')
    
    def _read_existing_log_files(self):
        '''
        This method reads the existing log files and sets 'self.init_pop'
        and 'self.init_pareto_front'.
        '''
        # read log file of all generations
        with open(self.log_pool_path) as f:
            log_content = f.readlines()
            
        if not log_content[0].startswith('#************************************CONFIG'):
            raise AttributeError("ERROR: Old log file '{0}' has a false header!".format(self.log_pool_path))
        if not log_content[-1].startswith('#@ END'):
            raise AttributeError("ERROR: Old log file '{0}' has an incomplete last generation! '#@ END' is missing!".format(self.log_pool_path))
        # read pop_size
        pop_size = int(log_content[1].strip().split(' = ')[1])
        iterations = int(log_content[2].strip().split(' = ')[1])
        optimization_type = log_content[3].strip().split(' = ')[1]
        number_of_genes = int(log_content[4].strip().split(' = ')[1])
        number_of_scores = int(log_content[5].strip().split(' = ')[1])
        
        # reverse the log content -> children are the first ones
        log_content.reverse()
        
        # set starting iteration
        self.starting_iteration = int(log_content[1].split(' ')[0])
        self.optimization_type = optimization_type
        
        # check if it is already finished
        if iterations <= self.starting_iteration:
            sys.exit('Calculations are already finished!')
        
        # read them with a dedicated method
        pop = Population()
        pareto_pop = Population()

        pop.revive_from_log(self.log_pool_path)
        pareto_pop.revive_from_log(self.log_pareto_path)

        # set pareto front 
        self.init_pareto_front = pareto_pop
        
        # perform a selection of the new offspring
        self.init_pop = pop
        print('Finished processing of existing log files')
    

class Worker_Bee(multiprocessing.Process):

    def __init__(self, task_queue, result_queue):
        multiprocessing.Process.__init__(self)
        self.task_queue = task_queue
        self.result_queue = result_queue

    def run(self):
        proc_name = self.name
        while True:
            try:
                next_task = self.task_queue.get(timeout = 5)
                if next_task is None:
                    # shutdown
                    #print('{0}: Exiting!'.format(proc_name))
                    self.task_queue.task_done()
                    break
                answer = next_task()
                self.task_queue.task_done()
                self.result_queue.put(answer)
            except Empty:
                time.sleep(5)
        return

class Task(object):
    
    def __init__(self, indi, score_func, data_dict):
        self.indi = indi
        self.score_func = score_func
        self.data_dict = data_dict

    def __call__(self):
        return self.score_func(self.indi, self.data_dict)
   


class GA_Optimizer(object):
    '''
    This class implements a genetic algorithm optimizer. 
    '''
    minimize = 'min'
    maximize = 'max'

    def __init__(self, population_size = None, generations = None,
                 mutate_pop = None, score_pop = None, select_pop = None,
                 analyze_pop = None, optimization_type = 'min', init_data = None,
                 crossover = False, crossover_frac = None, data_dict = None,
                 run_id = 0, log_data = False, print_generation = False):
        '''
        population_size -> # individuals
        generations -> # generations
        mutate_pop -> supply a function, which copies the data of each
                        individual, mutates it and creates a new individual
                        with the new data and returns a new population
                        afterwards. If None is given, it will assume that one
                        has a derivative of this class where it is implemented
                        as a class method!
        score_pop -> supply a function, which scores a given population.
                        If None is given, it will assume that one has a 
                        derivative of this class where it is implemented as 
                        a class method!
        select_pop -> supply a function, which selects a new pool from the
                        three given arguments: parents, children, pareto_front
                        a further argument is 'self.minimization' because the 
                        selection is sensitive to minimization or maximization
        analyze_pop -> supply a function, which analyzes a given population 
                        from the six given arguments: parents, children
                        pareto_front, generation, optimization_type and
                        data_dict. The result should be a dictionary e.g.
                        {'eval_step': eval_step, 'hyper_vol_pop' : x,
                         hyper_vol_pareto':y}
        optimization_type -> Is this a minimization or a maximization?
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
        self.pop_size = population_size
        self.generations = generations
        # If None is given, it is assumed, there is a derivative of this class
        # which implements this functions as class methods
        if mutate_pop is not None:
            self.mutate_pop = mutate_pop
        # If None is given, it is assumed, there is a derivative of this class
        # which implements this functions as class methods
        if score_pop is not None:
            self.score_pop = score_pop
        self.select_pop = select_pop
        self.analyze_pop = analyze_pop
        self.optimization_type = optimization_type
        # check if there is enough data for the initialization
        if init_data is None and len(init_data) < self.pop_size:
            raise ValueError('Not enough data for initialization!')
        else:
            self.init_data = init_data
        self.crossover = crossover
        self.crossover_frac = crossover_frac
        self.data_dict = data_dict
        self.run_id = run_id
        
        self.print_generation = print_generation
        
        self.log_data = log_data
        self.log_file_name = 'log_pool.LOG.LOG'
        self.pareto_log_file_name = 'log_pareto.LOG'
        self.analysis_log_file_name = 'log_analysis.LOG'
        self.working_directory = 'log_pool_{0}_{1}_{2}'.format(self.pop_size,
                                                               self.generations,
                                                               self.run_id)
        self.working_directory = os.path.join(os.getcwd(),
                                              self.working_directory)
        self.current_dir = os.getcwd()
        
        # check if there is already a run with the given settings, this might
        # indicate that there was an error and it can be resumed
        if self.log_data is True:
            # check if working dir exists
            if not os.path.exists(self.working_directory):
                os.mkdir(self.working_directory)
            # switch into it
            os.chdir(self.working_directory)
            # is there already some logging?
            if(os.path.exists(self.log_file_name)
               and os.path.exists(self.pareto_log_file_name)):
                # yes? then read the results
                # this method sets:
                # self.init_pop 
                # self.init_pareto_front
                # self.starting_generation
                self._read_existing_log_files()
            else:
                self.starting_generation = 0
                self.init_pop = Population()
        else:
            self.starting_generation = 0
            self.init_pop = Population()  
            
            
    def start_ga(self):
        '''
        Start the ga.
        '''
        start_time = time.time()
        
        # only initialize a new population, if there is none left
        if self.init_pop.size() > 0:
            parents = self.init_pop
            children = Population()
            pareto_front = self.init_pareto_front
        else:
            # init population
            init_pop = Population()
            for creation_id, new_data in enumerate(self.init_data):
                if creation_id < self.pop_size:
                    new_individual = Individual(new_data, 0, creation_id)
                    init_pop.add(new_individual)
                else:
                    break
            
            parents = self.score_pop(init_pop, self.data_dict)
                
            # children are empty
            children = Population()
            
            # get first pareto_front
            pareto_front = parents.get_pareto_frontier(self.optimization_type)
            

        for current_generation in range(self.starting_generation,
                                        self.generations):
            if self.print_generation is True:
                print("@Generation: {0}/{1}".format(current_generation,
                                                    self.generations))
            
            if self.crossover is True:
                intermediate_parents = self._ga_crossover(parents)
            elif self.crossover is False:
                intermediate_parents = parents
            elif self.crossover is not True and self.crossover is not False:
                intermediate_parents = self.crossover(parents, self.data_dict)
            else:
                print("I do not know what case I did not cover for the crossover?")
           
            # get mutants, add one to generation, as they are part of the next
            # generation
            children = self.mutate_pop(intermediate_parents,
                                       current_generation + 1, self.data_dict)
            
            # score the children
            children = self.score_pop(children, self.data_dict)
            
            # update pareto_front, combine all Population first
            united_population = Population()
            united_population.add_population(parents)
            united_population.add_population(children)
            united_population.add_population(pareto_front)
            
            pareto_front = united_population.get_pareto_frontier(self.optimization_type)
            
            if self.log_data is True:
                # log parents, children and pareto front
                self._log_population(parents, children, pareto_front,
                                     current_generation)
                
            # analyze if a function has been supplied
            if self.analyze_pop is not None:
                self._write_analysis(current_generation,
                                     self.analyze_pop(parents, children,
                                                      pareto_front,
                                                      current_generation,
                                                      self.optimization_type,
                                                      self.data_dict))
                
            # select offspring
            parents = self.select_pop(parents, children, pareto_front,
                                      self.optimization_type)
        
        end_time = time.time()
        
        os.chdir(self.current_dir)
        print('GA finished after {0}s.'.format(end_time - start_time))
        
        self.final_pop = parents
        self.final_pareto_front = pareto_front

    def _ga_crossover(self, crossover_parents):
        '''
        This method implements a crossover operator. It performs a single 
        crossover and uses a probability of 0.7, that it is 0.3 of the
        population 'survive' the operator without beeing changed.
        '''
        # these are a copy of crossover_parents
        new_parents = Population()
        new_parents.individuals = crossover_parents.individuals[:]
        # Population which contains the generated offspring
        offspring = Population()
        
        while not new_parents.is_empty():
            indi_number_1 = np.random.randint(0, len(new_parents.individuals))
            indi_1 = new_parents.individuals.pop(indi_number_1)
            indi_number_2 = np.random.randint(0, len(new_parents.individuals))
            indi_2 = new_parents.individuals.pop(indi_number_2)
            
            if np.random.random() <= self.crossover_frac:
                # perform crossover
                # get the data
                first_data = indi_1.data
                second_data = indi_2.data
                
                # point to perform the crossover
                crossover_point = np.random.randint(0, len(first_data))
                
                new_first_data = (first_data[:crossover_point]
                                  + second_data[crossover_point:])
                new_second_data = (second_data[:crossover_point]
                                   + first_data[crossover_point:])
                                   
                # make new instances, otherwise the parents data would change
                # which we don't want -> source of error!
                new_indi_1 = Individual(new_first_data, indi_1.generation,
                                       indi_1.creation_id)
                new_indi_2 = Individual(new_second_data, indi_2.generation,
                                        indi_2.creation_id)
                
                offspring.add(new_indi_1)
                offspring.add(new_indi_2)
                
            else:
                # keep a copy of the original ones (it is safer this way)
                offspring.add(Individual(indi_1.data, indi_1.generation,
                                         indi_1.creation_id))
                offspring.add(Individual(indi_2.data, indi_2.generation,
                                         indi_2.creation_id))
        
        # return the new offspring, which now only contains copies of the
        # original ones and the crossover individuals
        return offspring

    def _log_population(self, parents, children, pareto_front, current_generation):
        '''
        Log populations.
        '''
        if current_generation == 0:
            # if the gene is a sequence, it counts as one
            if isinstance(parents.individuals[0].data, str):
                number_of_genes = 1
            else:
                number_of_genes = len(parents.individuals[0].data)
            # line for r-headers
            number_of_scores = len(parents.individuals[0].score)
            
            with open(self.log_file_name, 'w') as f:
                f.write('#************************************CONFIG************************************\n')
                f.write('#* population_size = {0}\n'.format(self.pop_size))
                f.write('#* generations = {0}\n'.format(self.generations))
                f.write('#* number_of_genes = {0}\n'.format(number_of_genes))
                f.write('#* number_of_scores = {0}\n'.format(number_of_scores))
                f.write('#**************************************END*************************************\n')
                f.write('{0}\n'.format(parents.individuals[0].get_r_log_string()))
               
            with open(self.pareto_log_file_name, 'w') as f:
                f.write('#************************************CONFIG************************************\n')
                f.write('#* population_size = {0}\n'.format(self.pop_size))
                f.write('#* generations = {0}\n'.format(self.generations))
                f.write('#* number_of_genes = {0}\n'.format(number_of_genes))
                f.write('#* number_of_scores = {0}\n'.format(number_of_scores))
                f.write('#**************************************END*************************************\n')
                f.write('{0}\n'.format(parents.individuals[0].get_r_log_string()))
        
        # log parents and children
        with open(self.log_file_name, 'a') as f:
            f.write('#@ Generation: {0}\n'.format(current_generation))
        
        self._write_population(parents, self.log_file_name, current_generation)
        self._write_population(children, self.log_file_name, current_generation)
        
        with open(self.log_file_name, 'a') as f:
            f.write('#@ END\n')
        
        # log pareto front
        with open(self.pareto_log_file_name, 'a') as f:
            f.write('#@ Generation: {0}\n'.format(current_generation))
        
        self._write_population(pareto_front, self.pareto_log_file_name,
                               current_generation)
        
        with open(self.pareto_log_file_name, 'a') as f:
            f.write('#@ END\n')
    
    def _write_population(self, population, filename, current_generation = None):
        '''
        Write one population.
        '''
        with open(filename, 'a') as f:
            for indi in population.individuals:
                f.write(indi.get_log_string(current_generation))
    
    def _read_existing_log_files(self):
        '''
        This method reads the existing log files and sets 'self.init_pop'
        and 'self.init_pareto_front'.
        '''
        # read log file of all generations
        with open(self.log_file_name) as f:
            log_content = f.readlines()
            
        # read pareto log file
        with open(self.pareto_log_file_name) as f:
            pareto_log_content = f.readlines()
        
        if not log_content[0].startswith('#************************************CONFIG'):
            raise AttributeError("ERROR: Old log file '{0}' has a false header!".format(self.log_file_name))
        if not log_content[-1].startswith('#@ END'):
            raise AttributeError("ERROR: Old log file '{0}' has an incomplete last generation! '#@ END' is missing!".format(self.log_file_name))
        # read pop_size
        pop_size = int(log_content[1].strip().split(' = ')[1])
        generations = int(log_content[2].strip().split(' = ')[1])
        number_of_genes = int(log_content[3].strip().split(' = ')[1])
        number_of_scores = int(log_content[4].strip().split(' = ')[1])
        
        # reverse the log content -> children are the first ones
        log_content.reverse()
        pareto_log_content.reverse()
        
        # set starting generation
        self.starting_generation = int(log_content[1].split(' ')[0]) + 1
        
        # check if it is already finished
        if generations <= self.starting_generation:
            sys.exit('Calculations are already finished!')
        
        children = Population()
        parents = Population()
        pareto_front = Population()
        
        # read children
        for line in log_content[1:1 + pop_size]:
            children.add(self._revive_individual(line, number_of_genes,
                                                 number_of_scores))
            
        # read parents
        for line in log_content[1 + pop_size:1 + 2 * pop_size]:
            parents.add(self._revive_individual(line, number_of_genes,
                                                number_of_scores))
            
        # read last pareto front
        for line in pareto_log_content[1:]:
            if line.startswith('#@'):
                break
            else:
                pareto_front.add(self._revive_individual(line, number_of_genes,
                                                         number_of_scores))
        # set pareto front 
        self.init_pareto_front = pareto_front
        
        # perform a selection of the new offspring
        self.init_pop = self.select_pop(parents, children, pareto_front,
                                        self.optimization_type)
    
    def _revive_individual(self, line, number_of_genes, number_of_scores):
        '''
        This method takes a line from the log file and returns and individual.
        
        1. unique_id
        2. data
        ...
        2. + number_of_genes
        ...
        N- number_of_scores = score
        '''
        line = line.split(' ')
        unique_id = line[1]
        generation = line[0]
        creation_id = unique_id.split('_')[-1]
        
        # get genes
        data = []
        for i in line[2:2 + number_of_genes]:
            ## check if it is a sequence or a number
            #try:
            #    data_item = float(i)
            #except ValueError:
            #    data_item = i
            #data.append(data_item)
            ## here should be no conversion, because there could be binaey
            # binary encoded genes
            data.append(i)
       
       # data ist always a list!!!
        
        # get scores
        scores = []
        for i in line[-number_of_scores:]:
            scores.append(float(i))
            
        new_individual = Individual(data, generation, creation_id)
        new_individual.score = scores
        
        return new_individual
    
    def _write_analysis(self, current_generation, analysis_dict):
        '''
        This function writes out the supplied dictionary. 
        '''
        if current_generation == 0:
            with open(self.analysis_log_file_name, 'w') as f:
                # keys should be strings!
                str2write = ' '.join(analysis_dict.keys())
                f.write(str2write + '\n')
        
        with open(self.analysis_log_file_name, 'a') as f:
            cache_list = []
            for value in analysis_dict.values():
                cache_list.append(str(value))
            cache_str = ' '.join(cache_list)
            f.write(cache_str + '\n')



class Individual(object):
    
    minimize = 'min'
    maximize = 'max'
    
    def __init__(self, data, iteration, creation_id, parent_path = 'x',
                 mutation_path = 'x'):
        '''
        This initializes an individual with the given data. The data has to be
        a list!
        '''
        

        self.parent_path = parent_path
        self.mutation_path = mutation_path
        self.score = None
        self.iteration = iteration
        self.creation_id = creation_id   
        self.update_unique_id(self.iteration, self.creation_id)
        if not isinstance(data, list):
            raise AttributeError("Given data '{0}' of individual '{1}' is not a list!".format(data, self.unique_id))
        self.data = data
        self.mutations = None     
        
        self.counted = False
        
        #---Reda modifications
        self.n_p = None
        self.S_p = []
        self.mean_crowding_distance = None
        
    def get_dominant_number(self, population, optimization_type):
        '''
        This method returns the dominant number of this individual for the 
        given population. The dominant number counts the individuals in the 
        given population which dominate this individual (dominates means 
        x_i < y_i for at least one i).
        
        Arguments:
         *- population            
         *- optimization_type     either 'min' or 'max'
        '''
        if self.score is None:
            raise AttributeError("This individual '{0}' does not have a score!!!".format(self.unique_id))
        
        counter = 0
        for indi in population.individuals:
            if indi.score != self.score:
                #if dominates(indi.score, self.score, optimization_type):
                #    counter += 1
                if is_absolutely_fitter(indi.score, self.score, optimization_type):
                    counter += 1

        
        return counter
        
    def get_log_string(self, iteration):
        '''
        This method returns a logging string for this individual.
        Data always has to be a list! Even when it is only a protein sequence.
        '''
        data_str = ' '.join(str(i) for i in self.data)
        score_str = ' '.join(str(i) for i in self.score)
        
        return '{0} {1} {2} {3}\n'.format(iteration,
                                          self.unique_id,
                                          data_str, score_str)
        
    def get_r_log_string(self):
        '''
        This method returns an r log string line that can be put at the top
        of the logging file.
        '''
        if isinstance(self.data, str):
            number_of_genes = 1
        else:
            number_of_genes = len(self.data)
        # line for r-headers
        number_of_scores = len(self.score)
        r_str = 'iteration unique_id'
        for i in range(number_of_genes):
            r_str = r_str + ' data{0}'.format(i + 1)
        for i in range(number_of_scores):
            r_str = r_str + ' score{0}'.format(i + 1)
        
        return r_str
    
    def update_unique_id(self, iteration, creation_id):
        self.iteration = iteration
        self.creation_id = creation_id
        self.unique_id = '{0}_{1}'.format(self.iteration, self.creation_id)
        
    def score_with_population(self, population, competitor_id):
        '''
        This method returns a distance score to a given population. It
        calculates the distance to each individual in the population and
        returns the mean.
        '''
        if self.score is None:
            raise AttributeError("This individual '{0}' does not have a score!!!".format(self.unique_id))
        

        # array of this individuals score
        my_score = np.array(self.score)
        
        # make a list of all scores
        pop_scores = []
        pop_scores_all = []
        for individual in population.individuals:
            # we wouldn't want to include this individual, as it would 
            # influence the result
            pop_scores_all.append(individual.score)
            if(individual.unique_id != self.unique_id
               and individual.unique_id != competitor_id):
                pop_scores.append(individual.score)
        
        # there is the possibility, that the population only contains the two
        # competing individuals and then pop_scores is []
        if len(pop_scores) == 0:
            return 0
        
        # transform to an array
        pop_scores = np.array(pop_scores)

        # append my_score for all scores, there is the possibility that one
        # calculates the distance to the pareto front and then the individual
        # is not included and an exception could be raised
        pop_scores_all.append(list(my_score))
        pop_scores_all = np.array(pop_scores_all)
        
        # normalize 
        energy_min = np.amin(pop_scores_all, 0)
        energy_max = np.amax(pop_scores_all, 0)
        normalize = energy_max - energy_min
        normalize[np.nonzero(normalize == 0)] = 1
        pop_scores = pop_scores / normalize
        my_score = my_score / normalize
        
        # square subtracted vector
        for i, item in enumerate(pop_scores):
            pop_scores[i] = (item - my_score) ** 2
            
        # add x,y,z,... elements along '1' axis
        # take sqrt
        # calculate mean
        # return
        
        #---Reda modification
        self.mean_crowding_distance = np.sqrt(pop_scores.sum(1)).mean()
        
        return np.sqrt(pop_scores.sum(1)).mean()
    
    def score_energy_with_population(self, population, competitor_id, charge = 1):
        '''
        This method returns the potential of the individual.
        '''
        if self.score is None:
            raise AttributeError('This individual does not have a score!!!')
        
        # no individual means no energy
        if len(population.individuals) == 0:
            return 0
        
        # get the sign of the charge
        charge = np.sign(charge)
        
        # array of this individuals score
        my_score = np.array(self.score)
        
        # make a list of all scores
        pop_scores = []
        pop_scores_all = [] 
        for individual in population.individuals:
            # we wouldn't want to include this individual, as it would 
            # influence the result
            pop_scores_all.append(individual.score)
            if(individual.unique_id != self.unique_id
               and individual.unique_id != competitor_id):
                pop_scores.append(individual.score)
        
        # there is the possibility, that the population only contains the two
        # competing individuals and then pop_scores is []
        if len(pop_scores) == 0:
            return 0
        
        # transform to an array
        pop_scores = np.array(pop_scores)
        
        # append my_score for all scores, there is the possibility that one
        # calculates the distance to the pareto front and then the individual
        # is not included and an exception could be raised
        pop_scores_all.append(my_score)

        np.seterr('raise')
        
        # normalize 
        energy_min = np.amin(pop_scores_all, 0)
        energy_max = np.amax(pop_scores_all, 0)
        normalize = (energy_max - energy_min)
        # if there is only one point or two points, which share one coordinate
        # then there can be a 0, set it to 1, otherwise there would be a 
        # division with 0
        normalize[np.nonzero(normalize == 0)] = 1
        pop_scores = pop_scores / normalize
        my_score = my_score / normalize
        
        # square subtracted vector
        for i, item in enumerate(pop_scores):
            pop_scores[i] = (item - my_score) ** 2
        
        distance = np.sqrt(pop_scores.sum(1))
        
        distance = distance[np.nonzero(distance != 0)]
        
        potential = (charge / distance).sum()
        
        return potential

def revive_population(log_path, iteration = None):
    '''
    This functions reads a log file and returns a Population with the desired
    population. The 'log_path' can be either a 'log_pool' file or a 
    'log_pareto' file.
    If no iteration is given, it returns the last iteration population, if the
    file is empty it raises an error.
    '''
    with open(log_path) as f:
        log_content = f.readlines()

    pop_size = int(log_content[1].strip().split(' = ')[1])
    iterations = int(log_content[2].strip().split(' = ')[1])
    optimization_type = log_content[3].strip().split(' = ')[1]
    number_of_genes = int(log_content[4].strip().split(' = ')[1])
    number_of_scores = int(log_content[5].strip().split(' = ')[1])

    def revive_individual(line, number_of_genes, number_of_scores):
        '''
        This method takes a line from the log file and returns and individual.
        
        1. unique_id
        2. data
        ...
        2. + number_of_genes
        ...
        N- number_of_scores = score
        '''
        line = line.split(' ')
        unique_id = line[1]
        iteration = unique_id.split('_')[0]
        creation_id = unique_id.split('_')[-1]
        
        # get genes
        data = []
        for i in line[2:2 + number_of_genes]:
            # we only work on strings
            data_item = i
            data.append(data_item)
        
        # get scores
        scores = []
        for i in line[-number_of_scores:]:
            scores.append(float(i))
            
        new_individual = Individual(data, iteration, creation_id)
        new_individual.score = scores
        
        return new_individual
    
    pop_start = None
    if iteration is None:
        # find the last iteration
        for line_num,line in enumerate(log_content):
            if line.startswith('#@ Iteration: '):
                pop_start = line_num
                iteration = int(line.split(':')[1])
    else:
        # search for the start of the population
        # it is seperated, so that there are fewer errors
        for pop_start, line in enumerate(log_content):
            if line.startswith('#@ Iteration: {0}'.format(iteration)):
                break
    
    if pop_start is None:
        raise AttributeError("No population has been logged in '{0}'!".format(log_path))

    # start reviving
    pop_counter = 0
    pop = Population()
    for line in log_content[pop_start:]:
        if line.startswith('{0}'.format(iteration)):
            indi = revive_individual(line, number_of_genes, number_of_scores) 
            pop.add(indi)
            pop_counter += 1
        
        if pop_counter >= pop_size or line.startswith('#@ END'):
            break

    # save iteration step
    pop.iteration = iteration
    return pop


class Population(object):
    
    minimize = 'min'
    maximize = 'max'
    
    def __init__(self):
        '''
        Constructor
        
        Creates a new Population object.
        '''
        self.individuals = []
        # can be used to store information about the iteration step
        self.iteration = None
    
    
    def remove_individual(self, unique_id):
        '''
        Removes the first occurance of the individual from the population. 
        In general it should not happen, that there are two individuals with
        the same unique_id in the population.
        
        The method returns True, if the individual could be deleted and
        False if the individual could not be found!
        '''
        list_id = None
        for i, individual in enumerate(self.individuals):
            if individual.unique_id == unique_id:
                if list_id is None:
                    list_id = i
                else:
                    print("Found multiple occurences of '{0}', deleting only the first one!".format(unique_id))
            
        if list_id is None:
            return False
        else:
            self.individuals.pop(list_id)
            return True
    
    
    def is_empty(self):
        '''
        Check if the population contains individuals.
        '''
        if len(self.individuals) == 0:
            return True
        else:
            return False
        
    def add(self, i):
        '''
        Add an individual to the population (Manuel returned a boolean value
        from the add call of vector).
        '''
        self.individuals.append(i)
        
    def add_population(self, population):
        '''
        Add another population to this one (that is a copy of the list of 
        individuals, the individuals are not copied!).
        '''
        self.individuals.extend(population.individuals[:])
        
    def size(self):
        '''
        Returns the length of the individual list.
        '''
        return len(self.individuals)

    def clone(self):
        '''Clone this population. Notice that the individuals will not be
        cloned!
        '''
        pop = Population()
        pop.add_population(self)
        return pop
        
    def get_score_list_2_plot(self):
        '''
        Returns a list with all scores for plotting.
        '''
        score_list = []
        for score in self.individuals[0].score:
            score_list.append([])
            
        for indi in self.individuals:
            if indi.score is None:
                raise ValueError('Score of individual {0} is None!'.format(indi.unique_id))
            for i, score in  enumerate(indi.score):
                score_list[i].append(score)
        
        return score_list
    
    def get_score_list(self):
        '''
        This method returns a list with all score items.
        '''
        score_list = []
        for indi in self.individuals:
            if indi.score is None:
                raise ValueError('Score of individual {0} is None!'.format(indi.unique_id))
            score_list.append(indi.score)
        
        return score_list
    
    def get_data_list_2_plot(self):
        '''
        Returns a list with the data of each individual.
        '''
        data_list = []
        for data in self.individuals[0].data:
            if isinstance(data, str):
                break
            else:
                data_list.append([])
        
        for indi in self.individuals:
            if isinstance(indi.data, str):
                data_list.append(indi.data)
            else:
                for i, data in enumerate(indi.data):
                    data_list[i].append(data)
        
        return data_list
    
    def get_hypervolume(self, optimization_type, reference_point = None):
        '''
        Arguments:
         *- optimization_type     either 'min' or 'max'
        '''
        if self.size() == 0:
            return 0
        
        if optimization_type == self.minimize:
            return hypervolume(self.get_score_list(), optimization_type,
                               reference_point)
        elif optimization_type == self.maximize:
            return hypervolume(self.get_score_list(), optimization_type,
                               reference_point)
        else:
            raise ValueError("Unkown optimization type: '{0}'".format(optimization_type))
    
    def get_pareto_frontier(self, optimization_type):
        '''
        This method returns the pareto front of the population as a new
        population that consists of all pareto members.
        
        Arguments:
         *- optimization_type     either 'min' or 'max'
        '''
        if optimization_type == self.minimize:
            #return self._get_minimized_pareto_frontier()
            return optimizer_cython.get_minimized_pareto_frontier(self)
        elif optimization_type == self.maximize:
            #return self._get_maximized_pareto_frontier()
            return optimizer_cython.get_maximized_pareto_frontier(self)
        else:
            raise ValueError("Unkown optimization type: '{0}'".format(optimization_type))
    
    def get_pareto_shells(self, optimization_type):
        '''
        This method returns a list of populations, which make up a pareto
        shell. 
        First item in the list is the first pareto front and the last is 
        of course the shell with the most dominated individuals.
        
        Arguments:
         *- optimization_type     either 'min' or 'max'
        '''
        pareto_shells = []
        new_population = Population()
        new_population.add_population(self)
        
        # in the worst case there are n-shells
        loop_stopper = new_population.size() + 1
        while new_population.size() > 0 and loop_stopper > 0:
            loop_stopper -= 1
            frontier = new_population.get_pareto_frontier(optimization_type)
            pareto_shells.append(frontier)
            for indi in frontier.individuals:
                new_population.remove_individual(indi.unique_id)
        
        return pareto_shells
    
    def get_dominant_number_list(self, optimization_type):
        '''
        This method returns a list, with the dominant number for each
        individual. The dominant number describes, by how many other
        individuals a given individual is dominated. For a member of the pareto
        front this value is 0.
        '''
        if optimization_type != self.minimize and optimization_type != self.maximize:
            raise ValueError("Unkown optimization type: '{0}'".format(optimization_type))

        #dominant_number_list = []
        #for indi in self.individuals:
        #    dominant_number_list.append(indi.get_dominant_number(self,
        #                                                         optimization_type))
        scores = np.array(self.get_score_list())
        dominant_number_list = optimizer_cython.get_dominant_number_list(scores, optimization_type)
        return dominant_number_list

    def get_hypervolume_contribution(self, optimization_type):
        '''
        This method returns the hypervolume contribution of each individual.
        The method uses the worst fitness score in all dimensions with respect
        to the type of optimization as the reference point.
        
        Notice that the calculation of the hypervolume contribution only 
        makes sense for a pareto front. Notice further that there are points
        in a pareto front that have 0.0 contribution to the hypervolume. These
        elements are the "corners" of the pareto front and should be prefered!
        If there is just one individual in the population, this method returns
        [0.0]:browse confirm wa
        as well!
        
        Arguments:
         *- optimization_type     either 'min' or 'max'
        '''
        # if there is only one individual, we return 0
        if self.size() == 1:
            return [0.0]
    
        ref = self.get_reference_point(optimization_type)
        ps = self.get_score_list()

        # calculate the volume for all members
        hypervol_all = self.get_hypervolume(optimization_type, ref)
        
        hypervolume_list = []        
        for i in range(self.size()):
            ps_i = list(ps[0:i]) + list(ps[i + 1:])
            hypervol_i = hypervolume(ps_i, optimization_type, ref)
            hypervolume_list.append(hypervol_all - hypervol_i)
            
        return hypervolume_list 

    def get_reference_point(self, optimization_type):
        '''
        Calculate the reference point of the population.
        '''
        ps = self.get_score_list()
        n = min([len(p) for p in ps])
        if optimization_type == self.minimize:
            ref = [max(ps, key = lambda x: x[o])[o] for o in range(n)]
        elif optimization_type == self.maximize:
            ref = [min(ps, key = lambda x: x[o])[o] for o in range(n)]
        else:
            raise ValueError("Unkown input for optimization_type: '{0}', should be 'min' or 'max'!".format(optimization_type))
        
        return ref

    def _get_minimized_pareto_frontier(self):
        '''
        This method finds the pareto front of the population and returns a new
        population that consists of all pareto members.
        '''
        if self.is_empty():
            return Population()
        
        results = []
        # put first individual in result list
        results.append(self.individuals[0])
        # iterate over all individuals
        for i, item in enumerate(self.individuals):
            # first individual is already inside result list
            if i != 0:
                # list to check, if individual dominates all result members
                non_dominated_list = []
                pop_list = []
                # compare with result list
                for j, jtem in enumerate(results):
                    # list to check, if solution dominates members of result list
                    dominated_list = []
                    # iterate over fitness values
                    for k, ktem in enumerate(item.score):
                        if item.score[k] < jtem.score[k]:
                            # 1 -> solution is smaller
                            dominated_list.append(1)
                        elif item.score[k] == jtem.score[k]:
                            # 2 -> solution is equal
                            dominated_list.append(2)
                        else:
                            # 0 -> solution is bigger
                            dominated_list.append(0)
                    if 0 not in dominated_list and 1 in dominated_list:
                        # solution dominates result list
                        # -> remove this one, mark it in pop_list
                        pop_list.append(j)
                        non_dominated_list.append(1)
                    elif 0 in dominated_list and 1  in dominated_list:
                        # solution has bigger and smaller/equal elements
                        non_dominated_list.append(1)
                    else:
                        non_dominated_list.append(0)
                # reverse list, so the last elements get popped first
                pop_list.reverse()
                for element in pop_list:
                    results.pop(element)
                # if solution is not dominated by results
                if 1 in non_dominated_list and 0 not in non_dominated_list:
                    # append solution
                    results.append(item)
        p = self.clone()
        p.individuals = results
        return p
    
    def _get_maximized_pareto_frontier(self):
        '''
        This method finds the pareto front of the population and returns a new
        population that consists of all pareto members.
        '''
        if self.is_empty():
            return Population()
        
        results = []
        # put first individual in result list
        results.append(self.individuals[0])
        # iterate over all individuals
        for i, item in enumerate(self.individuals):
            # first individual is already inside result list
            if i != 0:
                # list to check, if individual dominates all result members
                non_dominated_list = []
                pop_list = []
                # compare with result list
                for j, jtem in enumerate(results):
                    # list to check, if solution dominates members of result list
                    dominated_list = []
                    # iterate over fitness values
                    for k, ktem in enumerate(item.score):
                        if item.score[k] > jtem.score[k]:
                            # 1 -> solution is bigger
                            dominated_list.append(1)
                        elif item.score[k] == jtem.score[k]:
                            # 2 -> solution is equal
                            dominated_list.append(2)
                        else:
                            # 0 -> solution is smaller
                            dominated_list.append(0)
                    if 0 not in dominated_list and 1 in dominated_list:
                        # solution dominates result list
                        # -> remove this one, mark it in pop_list
                        pop_list.append(j)
                        non_dominated_list.append(1)
                    elif 0 in dominated_list and 1  in dominated_list:
                        # solution has bigger and smaller/equal elements
                        non_dominated_list.append(1)
                    else:
                        non_dominated_list.append(0)
                # reverse list, so the last elements get popped first
                pop_list.reverse()
                for element in pop_list:
                    results.pop(element)
                # if solution is not dominated by results
                if 1 in non_dominated_list and 0 not in non_dominated_list:
                    # append solution
                    results.append(item)
        p = self.clone()
        p.individuals = results
        return p

    def revive_from_log(self, log_path, iteration = None):
        '''
        This method revives a population from a given log file path. If this
        population is not empty, it will raise an error, so that it will not
        override any data.
        '''
        if self.size() > 0:
            raise AttributeError("Can not revive log file '{0}', because this population is not empty!".format(log_path))
        pop = revive_population(log_path, iteration)
        self.individuals = pop.individuals

