
:mod:`Wrapper` --- YYYYYYY
======================================================

.. module:: Wrapper
   :synopsis: YYYYYYYY.
.. moduleauthor:: Christoph Wilms <christoph.wilms@uni-due.de>
.. sectionauthor:: Jean-Noel Grad <jean-noel.grad@uni-due.de>


This module provides YYYYY.


.. _Wrapper-syntax:

Module Syntax
-------------

Empty.

.. _contents-of-module-Wrapper:

Module Contents
---------------

.. class:: Gromacs(object)

    This class performs a md simulation.

    .. method:: run()

        Set up the md, run equilibration and run the md.

    .. method:: _run_em_minimization()

        Docstring missing.

    .. method:: _run_production()

        Docstring missing.

    .. method:: _set_up_topology()

        Docstring missing.

    .. method:: _generate_box()

        Docstring missing.

    .. method:: _solvate_box()

        Docstring missing.

    .. method:: _add_counter_ions()

        Docstring missing.

    .. method:: _generate_tpr_input(mdp_file, gro_file, topol_file, tpr_output)

        Docstring missing.

    .. method:: _run_md(id_name)

        Docstring missing.

    .. method:: _write_mdp_files()

        Docstring missing.

.. class:: LatPack

    If a trajectory is calculated an already existing file of the same
    name is deleted.

    .. method:: run_fold()

        Run the program.

    .. method:: _run_fold(command)

        Docstring missing.

    .. method:: get_rmsd(fold_1, seq_1, fold_2, seq_2, pdb_path)

        This method converts the first structure and saves it as pdb_path,
        then it reads the data and converts the second and overrides
        pdb_path. After the calculation the rmsd is returned and the pdb_path
        will be removed.
        The second structure will be rotated

    .. method:: convert_fold(fold, seq, pdb_file)

        Docstring missing.

    .. method:: get_energy_and_fold_from_trajectory()

        :returns: A list with the energy as the first and the most populated
            fold as the second item.

    .. method:: get_min_energy_and_fold_from_trajectory()

        :returns: A list with the minimum energy in the trajectory as the first and the corresponding
            fold as the second item.

    .. method:: _get_fold_command(fold, step_size, moveSet)

        generate latFold run command

    .. method:: _get_conv_command(fold, seq)

        Docstring missing.

.. class:: Rosetta(object)

    Docstring missing.

    .. method:: score_pdb(pdb_path_list, energy_func = 'standard')

        This functions scores the given pdbs and returns their score as a
        list. The returned list should contain the energies in the same order
        as the given pdb_path_list

    .. method:: refine_pdb(pdb_path, new_path)

        This function calls rosetta to minimize the given pdb. The input and
        the output are superimposed afterwards.

    .. method:: calculate_ddG_monomer(refined_wt_path, new_sequence, new_path_prefix = None, iterations = 50, energy_func = 'standard')

        This function calls rosetta to mutate a given pdb. It is recommended to
        refine the structure before mutating it. The wild type should be copied
        to the 'new_dir' previously. The calculated ddG will be returned and is
        calculated like this:
        
        :math:`ddG = dG_{mut} - dG_{wt}`
        
        and is < 0 for stabilizing mutations and > 0 for destabilizing ones.
        
        If no 'new_path_prefix' is given, it will delete the mutant pdbs,
        otherwise they will be named like this:
        
        '{0}_{1}.pdb'.format(new_path_prefix, i) with i in range(iterations).

.. class:: FoldX(object)

    This class implements a method to run a FoldX simulation.

    .. attribute:: pdb_name

        refers to a structure that should be mutated and analyzed without '.pdb'!

    .. attribute:: sequence_parent

        the sequence of the pdb (pdb_name)

    .. attribute:: sequence_child

        a mutated version of sequence_parent

    .. attribute:: repair_flag

        indicates, if the pdb should be repaired at first

    .. attribute:: new_name

        the name of the scored structure

    .. attribute:: Temp

        the temperature of the experiment in K

    .. attribute:: pH

        the pH of the experiment

    .. attribute:: number_of_seeds

        the number of runs with different starting
        conditions (in the end the mean will be returned)

    .. attribute:: fixed_during_repair

        the residues which should be kept fixed.
        Format: (aminoacid)(chain)(number), e.g. LC7, RA17

    .. method:: _make_FoldX_run_Repair()

        Docstring missing.

    .. method:: _make_FoldX_run_Build()

        Docstring missing.

    .. method:: _make_FoldX_mut_List()

        Docstring missing.

    .. method:: _make_FoldX_pdb_List(pdb_name)

        Docstring missing.

    .. method:: _extract_FoldX_ddG()

        Docstring missing.

    .. method:: extract_FoldX_ddG()

        Docstring missing.

    .. method:: run_FoldX()

        Docstring missing.

.. class:: Eris(object)

    Docstring missing.

    .. method:: calculate_ddG()

        Run eris and return the score.

    .. method:: _get_mut_instructions(old_seq, new_seq)

        Get mutation instructions::

            old_seq: AYTP
            index:   1234
            new_seq: ARTN

            return '\"Y2R P4N\"'

