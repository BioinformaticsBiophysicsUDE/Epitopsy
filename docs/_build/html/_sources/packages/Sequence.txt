
:mod:`Sequence` --- YYYYYYY
======================================================

.. module:: Sequence
   :synopsis: YYYYYYYY.
.. moduleauthor:: Christoph Wilms <christoph.wilms@uni-due.de>
.. sectionauthor:: Jean-Noel Grad <jean-noel.grad@uni-due.de>


This module provides YYYYY.


.. _Sequence-syntax:

Module Syntax
-------------

Empty.

.. _contents-of-module-Sequence:

Module Contents
---------------

.. function:: make_alignment_upper_case(orig_aln)

    Uses the alignment and transforms all characters to upper case.

    :param orig_aln: biopython alignment

    :returns: an upper case biopython alignment

.. function:: make_alignment_clean(orig_aln)

    Uses the alignment, removes sequences with non standard amino acids and
    reduces the alignment size if unnecessary gaps are there (this may happen,
    if a sequence has been removed).

    :param orig_aln: biopython alignment

    :returns: a clean biopython alignment

.. function:: format_stockholm_to_fasta(input_path, output_path)

    Change the format of an stockholm formatted alignment ot an fasta formatted
    alignment.

    :param input_path: path to the stockholm formatted format
    :param output_path: path to the new fasta alignment

.. function:: get_pairwise_alignment(ref_seq, query_seq)

    :param ref_seq: string of the reference sequence
    :param query_seq: string of the query sequence

    :returns: a dictionary with keys ("ref_seq", "query_seq") and the
      corresponding aligned sequences. It uses the first calculated
      alignment, equally scored alignments are not considered.

.. function:: get_almost_conserved_columns(aln, cutoff=0.95)

    :param aln: biopython alignment

    :returns: a list with the columns, where the amino acid with the highest
      frequency has a frequency larger than the cutoff. The columns start with 1!

