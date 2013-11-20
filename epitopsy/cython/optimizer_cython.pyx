"""
@author: Christoph Wilms
"""
import numpy as np

# "cimport" is used to import special compile-time information
# about the numpy module (this is stored in a file numpy.pxd which is
# currently part of the Cython distribution).
cimport numpy as np
cimport cython

from libc.math cimport sqrt, round

# We now need to fix a datatype for our arrays. I've used the variable
# DTYPE for this, which is assigned to the usual NumPy runtime
# type info object.
DTYPE_vdw = np.int
DTYPE_esp = np.float
DTYPE_cmplx = np.complex

# "ctypedef" assigns a corresponding compile-time type to DTYPE_t. For
# every type in the numpy module there's a corresponding compile-time
# type with a _t-suffix.
ctypedef np.int_t DTYPE_vdw_t
ctypedef np.float_t DTYPE_esp_t
ctypedef np.complex_t DTYPE_cmplx_t

def get_dominant_number_list(np.ndarray[DTYPE_esp_t, ndim=2] score_list, o_type):
    cdef int list_len = score_list.shape[0]
    cdef int score_len = score_list.shape[1]
    cdef np.ndarray[DTYPE_vdw_t, ndim= 1] dom_list = np.zeros([list_len],DTYPE_vdw)
    cdef int dominating = 0
    cdef int counter
    for i from 0 <= i < list_len:
        counter = 0
        for j from 0 <= j < list_len:
            if i != j:
                dominating = 0
                for k from 0 <= k < score_len:
                    if o_type == 'min':
                        if score_list[j,k] > score_list[i,k]:
                            # j does not dominate i
                            dominating = 1
                    elif  o_type == 'max':
                        if score_list[j,k] < score_list[i,k]:
                            # j does not dominate i
                            dominating = 1

                if dominating == 0:
                    counter = counter + 1
                
        dom_list[i] = counter

    return dom_list


def get_minimized_pareto_frontier(pop):
    if pop.is_empty():
        return pop.clone()
    
    cdef int list_len = len(pop.individuals)
    cdef int score_len = len(pop.individuals[0].score)
    cdef int i = 0
    cdef int j = 0
    cdef int k = 0
    
    results = []
    # put first individual in result list
    results.append(pop.individuals[0])
    # iterate over all individuals
    for i in range(list_len):
        item = pop.individuals[i]
        # first individual is already inside result list
        if i != 0:
            # list to check, if individual dominates all result members
            non_dominated_list = []
            pop_list = []
            # compare with result list
            for j in range(len(results)):
                jtem = results[j]
                # list to check, if solution dominates members of result list
                dominated_list = []
                # iterate over fitness values
                for k in range(score_len):
                    ktem = jtem.score[k]
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
    p = pop.clone()
    p.individuals = results
    return p

def get_maximized_pareto_frontier(pop):
    '''
    This method finds the pareto front of the population and returns a new
    population that consists of all pareto members.
    '''
    if pop.is_empty():
        return pop.clone()
    
    cdef int list_len = len(pop.individuals)
    cdef int score_len = len(pop.individuals[0].score)
    cdef int i = 0
    cdef int j = 0
    cdef int k = 0
    
    results = []
    # put first individual in result list
    results.append(pop.individuals[0])
    # iterate over all individuals
    for i in range(list_len):
        item = pop.individuals[i]
        # first individual is already inside result list
        if i != 0:
            # list to check, if individual dominates all result members
            non_dominated_list = []
            pop_list = []
            # compare with result list
            for j in range(len(results)):
                jtem = results[j]
                # list to check, if solution dominates members of result list
                dominated_list = []
                # iterate over fitness values
                for k in range(score_len):
                    ktem = jtem.score[k]
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
    p = pop.clone()
    p.individuals = results
    return p

