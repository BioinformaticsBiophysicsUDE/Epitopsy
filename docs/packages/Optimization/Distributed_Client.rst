
:mod:`Distributed_Client` --- YYYYYYY
======================================================

.. module:: Distributed_Client
   :synopsis: YYYYYYYY.
.. moduleauthor:: Christoph Wilms <christoph.wilms@uni-due.de>
.. sectionauthor:: Jean-Noel Grad <jean-noel.grad@uni-due.de>


This module provides YYYYY.


.. _Distributed_Client-syntax:

Module Syntax
-------------

Empty.

.. _contents-of-module-Distributed_Client:

Module Contents
---------------

.. class:: CMOptimizer_Client(CMParticipant)

    This class will be run on the clients. 
    
    Usage:

    #. always sync a script 'run_evaluation.py' with the client manager
    #. this script should contain a function:
       ``evaluate(path_to_dir, seq_data, pdb_path)``
       which returns a list with the calculated scores. If no pdb is used,
       pdb_path is set to ``None``. The function also has to make sure, that
       the modified pdb is also written to the pdd_path, because it will
       be read after the execution. Please switch the working directory
       always to the new directory, this way the source folder will not 
       be filled up with trash!


    .. method:: secure_run()

        Method that keeps the clients going.

.. class:: CMOptimizer_ClientProcess(Process, CMOptimizer_Client)

    Docstring missing.


