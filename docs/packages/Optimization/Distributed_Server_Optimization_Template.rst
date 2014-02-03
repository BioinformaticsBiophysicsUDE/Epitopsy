
:mod:`Distributed_Server_Optimization_Template` --- YYYYYYY
===========================================================

.. module:: Distributed_Server_Optimization_Template
   :synopsis: YYYYYYYY.
.. moduleauthor:: Christoph Wilms <christoph.wilms@uni-due.de>
.. sectionauthor:: Jean-Noel Grad <jean-noel.grad@uni-due.de>


This module provides YYYYY.


.. _Distributed_Server_Optimization_Template-syntax:

Module Syntax
-------------

How to use this optimizer:
    #. supply a config file
    #. the config is supposed to contain more parameters, please add them to 
       'the data_dict'
    #. modify the 'score_pop' class method of CMOptimizer
    #. modify the 'mutate_pop' class method of CMOptimizer
    #. supply a 'run_evaluator.py' script which will be copied to all clients
       and calculates the fitness for a given individual. The script gets the
       the following input: ``unique_id``, ``pdb_content``, ``gene_data``, ``data_dict``

.. _contents-of-module-Distributed_Server_Optimization_Template:

Module Contents
---------------

.. class:: CMOptimizer(GA_Optimizer, CMObject)

    Properties, which have to be specified in the config-file:

    .. attribute:: pop_size

            number of individuals

    .. attribute:: generations

            number of generations

    .. attribute:: mutate_pop

            supply a function, which copies the data of each
            individual, mutates it and creates a new individual
            with the new data and returns a new population
            afterwards

    .. attribute:: score_pop

            supply a function, which scores a given population

    .. attribute:: select_pop

            supply a function, which selects a new pool from the
            three given arguments: parents, children, pareto_front
            a further argument is 'self.minimization' because the 
            selection is sensitive to minimization or maximization

    .. attribute:: analyze_pop

            supply a function, which analyzes a given population 
            from the four given arguments: parents, children,
            pareto_front and current_generation. The result should
            be a dictionary e.g.::

                {'gen': current_gen, 'hyper_vol_parents' : x, hyper_vol_pareto':y}

    .. attribute:: minimization

            is this a minimization or a maximization?

    .. attribute:: init_data

            list of data that is used for the initialization

    .. attribute:: crossover

            can be ``True``, ``False`` or a function. The function works
            on the parent population, therefore it is necessary to 
            create new individuals, because the parents should not 
            be modified. After the crossover this population will 
            be mutated. If no function is supplied and crossover is
            ``True``, the implemented function will be used

    .. attribute:: crossover_frac

            fraction of the parent population, which will be
            mutated with the built in crossover

    .. attribute:: data_dict

            additional data, which can be used in any given function,
            a dictionary seems to be a good idea.

    .. attribute:: run_id

            can be used to create unique directory names for the same
            generation and population sizes

    .. attribute:: log_data

            log data or not

    .. attribute:: print_generation

            either print the current generation or not

    .. method:: score_pop(unscored_pop, data_dict)

        This method is used to score a population. Here the distributed computing
        is used! The actual scoring is performed by the supplied scoring script,
        which is run by the clients. For the proper function the client returns a
        list with the scores

    .. method:: mutate_pop(population, current_generation, data_dict)

        This function is used to mutate a given population.

    .. method:: secure_run()

        Docstring missing.

.. class:: CMOptimizerProcess(Process, CMOptimizer)

    Docstring missing.


