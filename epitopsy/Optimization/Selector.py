'''
Created on Jan 17, 2012

@author: chris
'''

import numpy as np

from epitopsy.Optimization.Optimizer import Population, Individual


#===============================================================================
# Tournaments
#===============================================================================

def random_tournament(parents, children, pareto_front, optimization_type):
    """
    This selector chooses two individuals randomly from the combined 
    population of parents and their mutated offspring and compares their
    scores and the fitter individual is selected. In case that both are
    equally fit, one gets randomly chosen. The individuals compete until 
    no individuals are left. 
    """
    p = Population()
    result = Population()
    p.add_population(parents)
    p.add_population(children)
    #=======================================================================
    # iterate over half of the combined population
    #=======================================================================
    num_elements = len(parents.individuals)
    for i in range(num_elements):
        choose_1 = np.random.randint(0, len(p.individuals))
        c1 = p.individuals.pop(choose_1)
        choose_2 = np.random.randint(0, len(p.individuals))
        c2 = p.individuals.pop(choose_2)
        # find the fitter individual, choose randomly, if none dominates
        # the other
        fitter_individual = who_is_fitter_random(c1, c2, optimization_type)
        # add the fitter one to the returned population
        result.add(fitter_individual)
        
    return result

def nsga_tournament(parents, children, pareto_front, optimization_type):
    '''
    Replaces population using the non-dominated sorting technique from
    NSGA-II. Fill new parent population according to the best front
    respectively crowding distance within a front 
    '''
    pool = Population()
    pool.add_population(parents)
    pool.add_population(children)
    survivors = Population()
    
    # give each individual their attributes:
    for indi in pool.individuals:
        indi.crowding_distance = 0
        indi.pareto_rank = -1

    # get a list filled with the pareto fronts
    pareto_shells = pool.get_pareto_shells(optimization_type)

    # calculate the crowding distance for each individual (yes, this can be
    # optimized,but this way it can be used by other methods as well
    # -> copy + paste)
    for pareto_rank, shell in enumerate(pareto_shells):
        for objective, score in enumerate(shell.individuals[0].score):
            # sort by this objective
            shell.individuals.sort(key = lambda x: x.score[objective])
            for indi_index, indi in enumerate(shell.individuals):
                if indi_index == 0:
                    indi.crowding_distance = np.inf
                    indi.pareto_rank = pareto_rank + 1
                elif indi_index == shell.size() - 1:
                    indi.crowding_distance = np.inf
                    indi.pareto_rank = pareto_rank + 1
                else:
                    # individuals are sorted by the current objective and now we
                    # get the min and max value to normalize the fitness
                    diff_min = shell.individuals[0].score[objective]
                    diff_max = shell.individuals[-1].score[objective]
                    diff = diff_max - diff_min
                    if diff == 0:
                        diff = 1
                    distance_l = shell.individuals[indi_index - 1].score[objective]
                    distance_r = shell.individuals[indi_index + 1].score[objective]
                    indi.crowding_distance = (indi.crowding_distance
                                              + distance_l + distance_r)
    
    # fill survivors
    for shell in pareto_shells:
        if len(survivors.individuals) == len(parents.individuals):
            break

        elif len(survivors.individuals) + shell.size() > len(parents.individuals):
            # we can not add all elements, so we sort them by their crowding
            # distance and use the best until survivors is full
            shell.individuals.sort(key = lambda x: x.crowding_distance)

            i = 0
            while len(survivors.individuals) < len(parents.individuals):
                survivors.add(shell.individuals[i])
                i += 1

        else:
            for indi in shell.individuals:
                survivors.add(indi)

    return survivors





def nsga_random_tournament(parents, children, pareto_front, optimization_type):
    """
    Replaces population using the non-dominated sorting technique from NSGA-II.
    Fill new parent population according to the best front respectively crowding distance within a front 
    """
    pool = Population()
    pool.add_population(parents)
    pool.add_population(children)
    survivors = Population()
    
    # give each individual their attributes:
    for indi in pool.individuals:
        indi.crowding_distance = 0
        indi.pareto_rank = -1

    # get a list filled with the pareto fronts
    pareto_shells = pool.get_pareto_shells(optimization_type)

    # calculate the crowding distance for each individual (yes, this can be
    # optimized,but this way it can be used by other methods as well
    # -> copy + paste)
    for pareto_rank, shell in enumerate(pareto_shells):
        for objective, score in enumerate(shell[0].score):
            # sort by this objective
            shell.sort(key = lambda x: x.score[objective])
            for indi_index, indi in enumerate(shell):
                if indi_index == 0:
                    indi.crowding_distance = np.inf
                    indi.pareto_rank = pareto_rank + 1
                elif indi_index == len(shell) - 1:
                    indi.crowding_distance = np.inf
                    indi.pareto_rank = pareto_rank + 1
                else:
                    # individuals are sorted by the current objective and now we
                    # get the min and max value to normalize the fitness
                    diff_min = shell[0].score[objective]
                    diff_max = shell[-1].score[objective]
                    diff = diff_max - diff_min
                    if diff == 0:
                        diff = 1
                    distance_l = shell[indi_index - 1].score[objective]
                    distance_r = shell[indi_index + 1].score[objective]
                    indi.crowding_distance = (indi.crowding_distance
                                              + distance_l + distance_r)


    for i in range(len(parents.individuals)):
        choose_1 = np.random.randint(0, len(pool.individuals))
        c1 = pool.individuals.pop(choose_1)
        choose_2 = np.random.randint(0, len(pool.individuals))
        c2 = pool.individuals.pop(choose_2)
        
        if c1.pareto_rank < c2.pareto_rank:
            # c1 is better
            survivors.add(c1)
        elif c1.pareto_rank > c2.pareto_rank:
            # c2 is better
            survivors.add(c2)
        else:
            # both are in the same shell
            if c1.crowding_distance > c2.crowding_distance:
                # c1 should be prefered because it is lonley, there are 
                # no friends nearby
                survivors.add(c1)
            elif c1.crowding_distance < c2.crowding_distance:
                # c2 should be prefered
                survivors.add(c2)
            else:
                # very unlikely, both are equal -> choose random
                if np.random.randint(0, 2) == 1:
                    # c1 won the lottery
                    survivors.add(c1)
                else:
                    # c2 has won
                    survivors.add(c2)
                    
        
    
    return survivors

def pareto_energy_tournament(parents, children, pareto_front, optimization_type):
    """
    This selector chooses two individuals randomly from the combined 
    population of parents and their mutated offspring and compares their
    scores and the fitter individual is selected. In case that both are
    equally fit, their energy towards the pareto front will be caluclated and
    the one with the lower energy will be chosen. The individuals compete until 
    no individuals are left. 
    """
    p = Population()
    result = Population()
    p.add_population(parents)
    p.add_population(children)
    #=======================================================================
    # iterate over half of the combined population
    #=======================================================================
    num_elements = len(parents.individuals)
    for i in range(num_elements):
        choose_1 = np.random.randint(0, len(p.individuals))
        c1 = p.individuals.pop(choose_1)
        choose_2 = np.random.randint(0, len(p.individuals))
        c2 = p.individuals.pop(choose_2)
        # find the fitter individual, choose randomly, if none dominates
        # the other
        fitter_individual = who_is_fitter_min_energy_all(c1, c2,
                                                              pareto_front,
                                                              optimization_type)
        # add the fitter one to the returned population
        result.add(fitter_individual)
        
    return result


def pool_pareto_distance_tournament(parents, children, pareto_front,
                                    optimization_type):
    """
    This selector chooses two individuals randomly from the combined 
    population of parents and their mutated offspring and compares their
    scores and the fitter individual is selected. In case that both are
    equally fit, the one with the greater distance to . The individuals compete until 
    no individuals are left. 
    """
    p = Population()
    result = Population()
    p.add_population(parents)
    p.add_population(children)
    
    reference_pop = Population()
    reference_pop.add_population(parents)
    reference_pop.add_population(children)
    reference_pop.add_population(pareto_front)
    #=======================================================================
    # iterate over half of the combined population
    #=======================================================================
    num_elements = len(parents.individuals)
    for i in range(num_elements):
        choose_1 = np.random.randint(0, len(p.individuals))
        c1 = p.individuals.pop(choose_1)
        choose_2 = np.random.randint(0, len(p.individuals))
        c2 = p.individuals.pop(choose_2)
        # find the fitter individual, choose randomly, if none dominates
        # the other
        fitter_individual = who_is_fitter_max_distance_all(c1, c2,
                                                           reference_pop,
                                                           optimization_type)
        # add the fitter one to the returned population
        result.add(fitter_individual)
        
    return result



#===============================================================================
# def dominate(c1, c2, minimize):
#    
#    if minimize is True:
#        relation_case_0 = 1
#        relation_case_1 = 2
#        relation_case_2 = 0
#    else:
#        relation_case_0 = 0
#        relation_case_1 = 2
#        relation_case_2 = 1
#    # list to find the fitter individual
#    compare_list = []
#    for j in range(len(c1.score)):
#        if c1.score[j] < c2.score[j]:
#            compare_list.append(relation_case_0)
#        elif c1.score[j] == c2.score[j]:
#            compare_list.append(relation_case_1)
#        elif c1.score[j] > c2.score[j]:
#            compare_list.append(relation_case_2)
#            
#    #===================================================================
#    # check the relation of the individuals 
#    #===================================================================
#    if 0 not in compare_list and 1 in compare_list:
#        # c1 is smaller than c2
#        return 1
#    elif 1 not in compare_list and 0 in compare_list:
#        # c1 is bigger than c2
#        return 2
#    elif 2 in compare_list and 0 not in compare_list and 1 not in compare_list:
#        # c1 equal c2
#        return 0
#===============================================================================


#===============================================================================
# These functions determine which individual in a tournament is fitter
#===============================================================================

def who_is_fitter_random(c1, c2, optimization_type):
        '''
        This method returns the fitter individual form the two candidates
        c1 and c2. In case none dominates the other, one is randomly chosen.
        '''
        #===================================================================
        # For minimization:
        # append 1 if c1.value is smaller than c2's
        # append 2 if c1.value is equal to c2's
        # append 0 if C1.value is bigger than c2's
        # for maximization 1 is switched with 0
        #===================================================================
        p = Population()
        if optimization_type == p.minimize:
            relation_case_0 = 1
            relation_case_1 = 2
            relation_case_2 = 0
        elif optimization_type == p.maximize:
            relation_case_0 = 0
            relation_case_1 = 2
            relation_case_2 = 1
        else:
            raise ValueError("Unkown input for optimization_type: '{0}', should be 'min' or 'max'!".format(optimization_type))
        # list to find the fitter individual
        compare_list = []
        for j in range(len(c1.score)):
            if c1.score[j] < c2.score[j]:
                compare_list.append(relation_case_0)
            elif c1.score[j] == c2.score[j]:
                compare_list.append(relation_case_1)
            elif c1.score[j] > c2.score[j]:
                compare_list.append(relation_case_2)
        
        #===================================================================
        # check the relation of the individuals 
        #===================================================================
        if 0 not in compare_list and 1 in compare_list:
            # c1 is smaller than c2
            return c1
        elif 1 not in compare_list and 0 in compare_list:
            # c1 is bigger than c2
            return c2
        else:
            if np.random.randint(0, 2) == 1:
                return c1
            else:
                return c2


def who_is_fitter_max_distance_all(c1, c2, population, optimization_type):
    '''
    This method returns the fitter individual form the two candidates
    c1 and c2. In case none dominates the other, the distance to all other
    individual is calculated and the one with the greater distance is 
    chosen to increase the diversity.
    '''
    #===================================================================
    # For minimization:
    # append 1 if c1.value is smaller than c2's
    # append 2 if c1.value is equal to c2's
    # append 0 if C1.value is bigger than c2's
    # for maximization 1 is switched with 0
    #===================================================================
    p = Population()
    if optimization_type == p.minimize:
        relation_case_0 = 1
        relation_case_1 = 2
        relation_case_2 = 0
    elif optimization_type == p.maximize:
        relation_case_0 = 0
        relation_case_1 = 2
        relation_case_2 = 1
    else:
        raise ValueError("Unkown input for optimization_type: '{0}', should be 'min' or 'max'!".format(optimization_type))
    # list to find the fitter individual
    compare_list = []
    for j in range(len(c1.score)):
        if c1.score[j] < c2.score[j]:
            compare_list.append(relation_case_0)
        elif c1.score[j] == c2.score[j]:
            compare_list.append(relation_case_1)
        elif c1.score[j] > c2.score[j]:
            compare_list.append(relation_case_2)
    
    #===================================================================
    # check the relation of the individuals 
    #===================================================================
    if 0 not in compare_list and 1 in compare_list:
        # c1 is smaller than c2
        return c1
    elif 1 not in compare_list and 0 in compare_list:
        # c1 is bigger than c2
        return c2
    else:
        distance_c1 = c1.score_with_population(population, c2.unique_id)
        distance_c2 = c2.score_with_population(population, c1.unique_id)
        if distance_c1 > distance_c2:
            return c1
        elif distance_c1 < distance_c2:
            return c2
        # very unlikely but you never know
        else:
            if np.random.randint(0, 2) == 1:
                return c1
            else:
                return c2

def who_is_fitter_min_energy_all(c1, c2, population, optimization_type):
    '''
    This method returns the fitter individual form the two candidates
    c1 and c2. In case none dominates the other, the distance to all other
    individual is calculated and the one with the greater distance is 
    chosen to increase the diversity.
    '''
    #===================================================================
    # For minimization:
    # append 1 if c1.value is smaller than c2's
    # append 2 if c1.value is equal to c2's
    # append 0 if C1.value is bigger than c2's
    # for maximization 1 is switched with 0
    #===================================================================
    p = Population()
    if optimization_type == p.minimize:
        relation_case_0 = 1
        relation_case_1 = 2
        relation_case_2 = 0
    elif optimization_type == p.maximize:
        relation_case_0 = 0
        relation_case_1 = 2
        relation_case_2 = 1
    else:
            raise ValueError("Unkown input for optimization_type: '{0}', should be 'min' or 'max'!".format(optimization_type))
    # list to find the fitter individual
    compare_list = []
    for j in range(len(c1.score)):
        if c1.score[j] < c2.score[j]:
            compare_list.append(relation_case_0)
        elif c1.score[j] == c2.score[j]:
            compare_list.append(relation_case_1)
        elif c1.score[j] > c2.score[j]:
            compare_list.append(relation_case_2)
    
    #===================================================================
    # check the relation of the individuals 
    #===================================================================
    if 0 not in compare_list and 1 in compare_list:
        # c1 is smaller than c2
        return c1
    elif 1 not in compare_list and 0 in compare_list:
        # c1 is bigger than c2
        return c2
    else:
        
        
        energy_c1 = c1.score_energy_with_population(population, c2.unique_id)
        energy_c2 = c2.score_energy_with_population(population, c1.unique_id)
        if energy_c1 < energy_c2:
            return c1
        elif energy_c1 > energy_c2:
            return c2
        # very unlikely but you never know
        else:
            if np.random.randint(0, 2) == 1:
                return c1
            else:
                return c2
