'''
Created on May 8, 2012

@author: chris
'''

import math

def hypervolume(pareto_set, optimization_type, reference_point = None):
    '''Calculates the hypervolume by slicing objectives (HSO).
    
    The calculation has been 'stolen' from inspyred, but it is open source
    and so is this code :-)
    
    Arguments:
    
    - *pareto_set* -- the list or lists of objective values comprising the Pareto front
    - *optimization_type* -- 'min' for minimization or 'max' for maximization
    - *reference_point* -- the reference point to be used (default None)
    
    '''
    if optimization_type == 'min':
        return hypervolume_min(pareto_set, reference_point)
    elif optimization_type == 'max':
        return hypervolume_max(pareto_set, reference_point)
    else:
        raise ValueError("Unkown input for optimization_type: '{0}', should be 'min' or 'max'!".format(optimization_type))

def hypervolume_max(pareto_set, reference_point = None):
    """Calculates the hypervolume by slicing objectives (HSO).
    
    This function calculates the hypervolume (or S-measure) of a nondominated
    set using the Hypervolume by Slicing Objectives (HSO) procedure of `While, et al. 
    (IEEE CEC 2005) <http://www.lania.mx/~ccoello/EMOO/while05a.pdf.gz>`_.
    The *pareto_set* should be a list of lists of objective values.
    The *reference_point* may be specified or it may be left as the default 
    value of None. In that case, the reference point is calculated to be the
    maximum value in the set for all objectives (the ideal point). This function 
    assumes that objectives are to be maximized.
    
    Arguments:
    
    - *pareto_set* -- the list or lists of objective values comprising the Pareto front
    - *reference_point* -- the reference point to be used (default None)
    
    """
    def dominates(p, q, k = None):
        if k is None:
            k = len(p)
        d = True
        while d and k < len(p):
            d = not (q[k] > p[k])
            k += 1
        return d
        
    def insert(p, k, pl):
        ql = []
        while pl and pl[0][k] > p[k]:
            ql.append(pl[0])
            pl = pl[1:]
        ql.append(p)
        while pl:
            if not dominates(p, pl[0], k):
                ql.append(pl[0])
            pl = pl[1:]
        return ql

    def slice(pl, k, ref):
        p = pl[0]
        pl = pl[1:]
        ql = []
        s = []
        while pl:
            ql = insert(p, k + 1, ql)
            p_prime = pl[0]
            s.append((math.fabs(p[k] - p_prime[k]), ql))
            p = p_prime
            pl = pl[1:]
        ql = insert(p, k + 1, ql)
        s.append((math.fabs(p[k] - ref[k]), ql))
        return s

    ps = pareto_set
    ref = reference_point
    n = min([len(p) for p in ps])
    if ref is None:
        ref = [min(ps, key = lambda x: x[o])[o] for o in range(n)]
    pl = ps[:]
    pl.sort(key = lambda x: x[0], reverse = True)
    s = [(1, pl)]
    for k in range(n - 1):
        s_prime = []
        for x, ql in s:
            for x_prime, ql_prime in slice(ql, k, ref):
                s_prime.append((x * x_prime, ql_prime))
        s = s_prime
    vol = 0
    for x, ql in s:
        vol = vol + x * math.fabs(ql[0][n - 1] - ref[n - 1])
    return vol

def hypervolume_min(pareto_set, reference_point = None):
    """Calculates the hypervolume by slicing objectives (HSO).
    
    This function calculates the hypervolume (or S-measure) of a nondominated
    set using the Hypervolume by Slicing Objectives (HSO) procedure of `While, et al. 
    (IEEE CEC 2005) <http://www.lania.mx/~ccoello/EMOO/while05a.pdf.gz>`_.
    The *pareto_set* should be a list of lists of objective values.
    The *reference_point* may be specified or it may be left as the default 
    value of None. In that case, the reference point is calculated to be the
    maximum value in the set for all objectives (the ideal point). This function 
    assumes that objectives are to be minmized. 

    ###
    
    I have hijacked the function and now it calculates the correct hypervolume 
    for a minimization!

    ###
    
    Arguments:
    
    - *pareto_set* -- the list or lists of objective values comprising the Pareto front
    - *reference_point* -- the reference point to be used (default None)
    
    """
    def dominates(p, q, k = None):
        if k is None:
            k = len(p)
        d = True
        while d and k < len(p):
            d = not (q[k] < p[k])
            k += 1
        return d

        
    def insert(p, k, pl):
        ql = []
        #while pl and pl[0][k] > p[k]: -> changed for minimization
        while pl and pl[0][k] < p[k]:
            ql.append(pl[0])
            pl = pl[1:]
        ql.append(p)
        while pl:
            if not dominates(p, pl[0], k):
                ql.append(pl[0])
            pl = pl[1:]
        return ql

    def slice(pl, k, ref):
        p = pl[0]
        pl = pl[1:]
        ql = []
        s = []
        while pl:
            ql = insert(p, k + 1, ql)
            p_prime = pl[0]
            s.append((math.fabs(p[k] - p_prime[k]), ql))
            p = p_prime
            pl = pl[1:]
        ql = insert(p, k + 1, ql)
        s.append((math.fabs(p[k] - ref[k]), ql))
        return s

    ps = pareto_set
    ref = reference_point

    # dimension of the problem
    n = min([len(p) for p in ps])
    
    # find the reference point -> maximum in each dimension
    if ref is None:
        ref = [max(ps, key = lambda x: x[o])[o] for o in range(n)]
    
    # sort pareto
    pl = ps[:]
    pl.sort(key = lambda x: x[0])#, reverse = True)
    s = [(1, pl)]
    for k in range(n - 1):
        s_prime = []
        for x, ql in s:
            for x_prime, ql_prime in slice(ql, k, ref):
                s_prime.append((x * x_prime, ql_prime))
        s = s_prime
    vol = 0
    for x, ql in s:
        vol = vol + x * math.fabs(ql[0][n - 1] - ref[n - 1])
    return vol
