
:mod:`Modeller` --- YYYYYYY
======================================================

.. module:: Modeller
   :synopsis: YYYYYYYY.
.. moduleauthor:: Christoph Wilms <christoph.wilms@uni-due.de>
.. sectionauthor:: Jean-Noel Grad <jean-noel.grad@uni-due.de>


This module provides YYYYY.


.. _Modeller-syntax:

Module Syntax
-------------

Empty.

.. _contents-of-module-Modeller:

Module Contents
---------------

.. function:: align_structures_biopython(struct_path_ref, struct_path_query, new_query_path)

    Docstring missing.

.. function:: get_alignment(pdb_ref, pdb_query)

    Subroutine of :func:`align_structures_biopython`.

.. function:: get_sequence(pdb)

    Subroutine of :func:`align_structures_biopython`.

.. function:: built_homo_multimer_model(ref_struct_path, query_seq, n_models, new_query_name="query", pir_file="aln.pir", script_file="model.py", chain_id='A', new_folder="build_model", debug_mode=False, check_chainbreaks=True)

    This works only for homo multimers, where the sequence is the same
    for each chain! It renumbers the modelled structure, this may cause
    a difference for HETATM!

    :param ref_struct_path: path to the structure that contains one or more
      chains with the same sequence. For every non standard amino acid,
      a '.' is used (only the first chain is checked!).
    :param query_seq: string, which contains the modelled sequence
    :param new_query_name: names the final structures, i.e. "query_name_1.pdb"
    :param chain_id: this chain is used to extract the sequence information

    :returns: a list with the new structure paths.
            
.. function:: write_pir(ref_seq,ref_struct_path, query_seq, pir_file, chain_ids, n_chains, start_res_id, end_res_id)

    Docstring missing.

.. function:: write_homo_multimer_script(pir_file, ref_knowns, query_name, script_file, n_models)

    Docstring missing.

.. function:: built_monomer_model(ref_struct_path, query_seq, n_models, new_query_name="query", pir_file="aln.pir", script_file="model.py", chain_id='A', new_folder="build_model", debug_mode=False)

    This works only for monomers!

    :param ref_struct_path: path to the structure that contains chain X with
            the correct sequence, which is extracted during the
            calculation. For every non standard amino acid, a
            '.' is used.
    :param query_seq: string, which contains the modelled sequence
    :param new_query_name: names the final structures, i.e. "query_name_1.pdb"

    :returns: A list with the new structure paths.
            
.. function:: write_pir(ref_seq,ref_struct_path, query_seq, pir_file,chain_id)

    Docstring missing.

.. function:: write_monomer_script(pir_file, ref_knowns, query_name, script_file, n_models)

    Docstring missing.

.. function:: built_dimer_model(ref_struct_path, query_seq, n_models, new_query_name="query", pir_file="aln.pir", script_file="model.py", chain_id='A', new_folder="build_model", keep_both=False, debug_mode=False)

    This works only for homodimers, where one only uses the constraints from
    the homodimer as constraints for the conformational flexibility, e.g.
    constraints for the flexibility of the N-terminal loop of hedgehogs.
    The sequence has to be the same for both chains, of course!

    :param ref_struct_path: path to the structure that contains chain A+B with
      the correct sequence, which is extracted during the calculation.
      For every non standard amino acid, a '.' is used.
    :param query_seq: string, which contains the modelled sequence
    :param new_query_name: names the final structures, i.e. "query_name_1.pdb"

    :returns: a list with the new structure paths.
            
.. function:: write_pir(ref_seq,ref_struct_path, query_seq, pir_file)

    Docstring missing.

.. function:: write_dimer_script(pir_file, ref_knowns, query_name, script_file, n_models)

    Docstring missing.

.. function:: special_restraints(self, aln)

    Docstring missing.

.. function:: user_after_single_model(self)

    Docstring missing.

