import os
import shutil
import numpy as np
import warnings
import subprocess

from Bio.SCOP.Raf import to_one_letter_code
from Bio.Alphabet import generic_protein
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.PDB import PDBParser, StructureAlignment,Superimposer,PDBIO
from Bio.Align import MultipleSeqAlignment
from Bio.PDB.PDBExceptions import PDBConstructionWarning

from epitopsy.Structure import PDBFile
from epitopsy.Sequence import get_pairwise_alignment
from epitopsy.tools.style import style

#### utility stuff
def align_structures_biopython(struct_path_ref, struct_path_query, new_query_path):
    def get_alignment(pdb_ref, pdb_query):
        seq_ref = get_sequence(pdb_ref)
        seq_query = get_sequence(pdb_query)

        aligned = get_pairwise_alignment(seq_ref, seq_query)
        aln_ref = aligned["ref_seq"]
        aln_query = aligned["query_seq"]

        aln = MultipleSeqAlignment([SeqRecord(Seq(aln_ref, generic_protein), id="ref"),
                SeqRecord(Seq(aln_query, generic_protein), id="query")])
        return aln


    def get_sequence(pdb):
        seq = ""
        if len(pdb) > 1:
            raise ValueError("Can not handle structures with more than one MODEL!\nThis structure has {0} MODELS!".format(len(pdb)))
        if len(pdb[0]) > 1:
            raise ValueError("Can not handle structures with more than one CHAIN!\nThis structure has {0} CHAINS!".format(len(pdb[0])))
        for model in pdb:
            for chain in model:
                for res in chain:
                    if res.resname in to_one_letter_code:
                        seq = "{0}{1}".format(seq, to_one_letter_code[res.resname])
        return seq

    struct_ref = struct_path_ref
    struct_query = struct_path_query

    parser = PDBParser()
    pdb_ref = parser.get_structure("ref", struct_ref)
    pdb_query = parser.get_structure("query", struct_query)

    aln = get_alignment(pdb_ref, pdb_query)

    coords_ref = []
    coords_query = []
    al=StructureAlignment(aln, pdb_ref, pdb_query)
    for (r1,r2) in al.get_iterator():
        if r1 is not None and r2 is not None:
            coords_ref.append(r1['CA'])
            coords_query.append(r2['CA'])

    coords_ref = np.array(coords_ref)
    coords_query = np.array(coords_query)

    super_imposer = Superimposer()
    super_imposer.set_atoms(coords_ref, coords_query)
    super_imposer.apply(pdb_query.get_atoms())
    io = PDBIO()
    io.set_structure(pdb_query)
    io.save(new_query_path)

#### homo multimer model
def built_homo_multimer_model(ref_struct_path, query_seq, n_models,
        new_query_name="query",
        pir_file="aln.pir", script_file="model.py", chain_id='A',
        new_folder="build_model", debug_mode=False,
        check_chainbreaks=True):
    '''
    This works only for homo multimeres, where the sequence is the same
    for each chain! Renumbers the modelled structure, this may cause
    a difference for HETATM!

    Args:
        ref_struct_path -> path to the structure that contains one or more
                        chains with the same sequence. For every non
                        standard amino acid, a '.' is used (only the first
                        chain is checked!).
        query_seq -> string, which contains the modelled sequence
        new_query_name -> names the final structures, i.e. "query_name_1.pdb"
        chain_id -> this chain is used to extract the sequence information

    Returns:
        A list with the new structure paths.
    '''
    wd = os.getcwd()
    if not os.path.exists(os.path.split(ref_struct_path)[-1]):
        shutil.copy(ref_struct_path, os.path.split(ref_struct_path)[-1])

    # make sure, it is not an absolute path
    ref_struct_path = os.path.split(ref_struct_path)[-1]

    # check chain break
    if check_chainbreaks:
        my_pdb = PDBFile(ref_struct_path)
        if my_pdb.contains_chain_break():
            raise AttributeError("Found chain breaks!")
    parser = PDBParser()
    pdb = parser.get_structure("dummy", ref_struct_path)
    if len(pdb) > 1:
        raise ValueError("Can not handle structures with more than one MODEL!\nThis structure has {0} MODELS!".format(len(pdb)))

    ref_seq = ""
    non_standard_seq = ""
    res_ids = []
    chain_ids = [ x.id for x in pdb[0].child_list[:]]
    n_chains = len(pdb[0])
    for res in pdb[0][chain_id]:
        res_ids.append(res.id[1])
        if res.resname in PDBFile.three2oneletter:
            ref_seq = "{0}{1}".format(ref_seq, PDBFile.three2oneletter[res.resname])
        else:
            # store the non standard amino acids in a variable
            non_standard_seq = "{0}.".format(non_standard_seq)

    start_res_id = min(res_ids)
    end_res_id = max(res_ids)

    aln_result = get_pairwise_alignment(ref_seq, query_seq)
    aln_ref = aln_result["ref_seq"]
    aln_query = aln_result["query_seq"]

    # add non standard amino acids to both string
    aln_ref = "{0}{1}".format(aln_ref, non_standard_seq)
    aln_query = "{0}{1}".format(aln_query, non_standard_seq)

    if os.path.exists(new_folder):
        raise ValueError("Folder already exists: {0}".format(new_folder))

    os.mkdir(new_folder)
    os.chdir(new_folder)
    shutil.copy(os.path.join(wd, ref_struct_path), ref_struct_path)

    ref_knowns = "ref"
    query_name = "query"

    def write_pir(ref_seq,ref_struct_path, query_seq, pir_file, chain_ids,
                  n_chains, start_res_id, end_res_id):
        all_ref_seq = []
        all_query_seq = []
        for i in range(n_chains):
            all_ref_seq.append(ref_seq)
            all_query_seq.append(query_seq)
        all_ref_seq = "/\n".join(all_ref_seq)
        all_query_seq = "/\n".join(all_query_seq)
        start_chain = sorted(chain_ids)[0]
        end_chain = sorted(chain_ids)[-1]
        pir_string = """>P1;ref
structureX:{struct_path}:{start_res}:{start_chain}:{end_res}:{end_chain}:::0.00:0.00
{ref_seq}*
>P1;query
sequence:query:::::::0.00:0.00
{query_seq}*"""
        if os.path.exists(pir_file):
            raise ValueError("Alingment file already exists: {0}".format(pir_file))
        with open(pir_file,"w") as f:
            f.write(pir_string.format(struct_path=ref_struct_path,
                                      ref_seq=all_ref_seq,
                                      query_seq=all_query_seq,
                                      start_res=start_res_id,
                                      end_res=end_res_id,
                                      start_chain=start_chain,
                                      end_chain=end_chain))
        return


    write_pir(aln_ref, ref_struct_path, aln_query,pir_file, chain_ids,
              n_chains, start_res_id, end_res_id)
    write_homo_multimer_script(pir_file, ref_knowns, query_name, script_file, n_models)
    p = subprocess.call(["python", script_file])
    if p == 0:
        # clean up
        file_list = os.listdir('.')
        file_list = filter(lambda x:x.startswith(query_name) and x.endswith("_fit.pdb"), file_list)
        file_list.sort()
        new_structures = []
        for i,item in enumerate(file_list):
            new_struct_path = os.path.join(wd,"{0}_{1}.pdb".format(new_query_name, i))
            parser = PDBParser()
            new_pdb = parser.get_structure("dummy", item)
            for chain in new_pdb[0]:
                for i,res in enumerate(chain.child_list):
                    old_id = list(res.id)
                     # we know that the structure has no chain breaks!
                    old_id[1] = start_res_id + i
                    res.id = tuple(old_id)

            if len(new_pdb[0].child_list) == 1:
                new_pdb[0].child_list[0].id = "A"
            io = PDBIO()
            io.set_structure(new_pdb)
            io.save(new_struct_path)
            #shutil.copy(item, new_struct_path)
            new_structures.append(new_struct_path)
        os.chdir(wd)

        ## align all structures
        cmd = ["python",style["pymol_align"],ref_struct_path]+new_structures
        #print(cmd)
        p = subprocess.call(cmd)

        if not debug_mode:
            shutil.rmtree(new_folder)

        return new_structures
    else:
        print("#### subprocess error")
        os.chdir(wd)



def write_homo_multimer_script(pir_file, ref_knowns, query_name, script_file, n_models):
    script_string = """from modeller import *
from modeller.automodel import *

log.verbose()

env = environ()
env.io.atom_files_directory = ['.', '../atom_files']
env.io.hetatm=True
a = automodel(env, alnfile='{pir_file}',
              knowns='{ref_knowns}', sequence='{query_name}',assess_methods=(assess.DOPE, assess.GA341))
a.final_malign3d=True
a.starting_model = 1
a.ending_model = {n_models}
a.make()
    """
    if os.path.exists(script_file):
        raise ValueError("Script file already exists: {0}".format(script_file))
    with open(script_file,"w") as f:
        f.write(script_string.format(pir_file=pir_file, ref_knowns=ref_knowns,
                query_name=query_name, n_models=n_models))
    return



#### monomer model
def built_monomer_model(ref_struct_path, query_seq, n_models,
        new_query_name="query",
        pir_file="aln.pir", script_file="model.py", chain_id='A',
        new_folder="build_model", debug_mode=False):
    '''
    This works only for monomers!

    Args:
        ref_struct_path -> path to the structure that contains chain X with
                        the correct sequence, which is extracted during the
                        calculation. For every non standard amino acid, a
                        '.' is used.
        query_seq -> string, which contains the modelled sequence
        new_query_name -> names the final structures, i.e. "query_name_1.pdb"

    Returns:
        A list with the new structure paths.
    '''
    wd = os.getcwd()
    if not os.path.exists(os.path.split(ref_struct_path)[-1]):
        shutil.copy(ref_struct_path, os.path.split(ref_struct_path)[-1])

    # make sure, it is not an absolute path
    ref_struct_path = os.path.split(ref_struct_path)[-1]

    parser = PDBParser()
    pdb = parser.get_structure("dummy", ref_struct_path)
    if len(pdb) > 1:
        raise ValueError("Can not handle structures with more than one MODEL!\nThis structure has {0} MODELS!".format(len(pdb)))

    ref_seq = ""
    non_standard_seq = ""
    res_ids = []
    for res in pdb[0][chain_id]:
        res_ids.append(res.id[1])
        if res.resname in PDBFile.three2oneletter:
            ref_seq = "{0}{1}".format(ref_seq, PDBFile.three2oneletter[res.resname])
        else:
            # store the non standard amino acids in a variable
            non_standard_seq = "{0}.".format(non_standard_seq)

    start_res_id = min(res_ids)

    aln_result = get_pairwise_alignment(ref_seq, query_seq)
    aln_ref = aln_result["ref_seq"]
    aln_query = aln_result["query_seq"]

    # add non standard amino acids to both string
    aln_ref = "{0}{1}".format(aln_ref, non_standard_seq)
    aln_query = "{0}{1}".format(aln_query, non_standard_seq)

    if os.path.exists(new_folder):
        raise ValueError("Folder already exists: {0}".format(new_folder))

    os.mkdir(new_folder)
    os.chdir(new_folder)
    shutil.copy(os.path.join(wd, ref_struct_path), ref_struct_path)

    ref_knowns = "ref"
    query_name = "query"

    def write_pir(ref_seq,ref_struct_path, query_seq, pir_file,chain_id):
        pir_string = """>P1;ref
structureX:{struct_path}::{chain_id}:::::0.00:0.00
{ref_seq}*
>P1;query
sequence:query:::::::0.00:0.00
{query_seq}*"""
        if os.path.exists(pir_file):
            raise ValueError("Alingment file already exists: {0}".format(pir_file))
        with open(pir_file,"w") as f:
            f.write(pir_string.format(struct_path=ref_struct_path, ref_seq=ref_seq,
                    query_seq=query_seq,chain_id=chain_id))
        return


    write_pir(aln_ref, ref_struct_path, aln_query,pir_file,chain_id)
    write_monomer_script(pir_file, ref_knowns, query_name, script_file, n_models)
    p = subprocess.call(["python", script_file])
    if p == 0:
        # clean up
        file_list = os.listdir('.')
        file_list = filter(lambda x:x.startswith(query_name) and x.endswith("_fit.pdb"), file_list)
        file_list.sort()
        new_structures = []
        for i,item in enumerate(file_list):
            new_struct_path = os.path.join(wd,"{0}_{1}.pdb".format(new_query_name, i))
            parser = PDBParser()
            new_pdb = parser.get_structure("dummy", item)
            ## debug
            for res in new_pdb[0][" "].child_list:
                old_id = list(res.id)
                old_id[1] += start_res_id - 1
                res.id = tuple(old_id)
            new_pdb[0][" "].id = "A"
            io = PDBIO()
            io.set_structure(new_pdb)
            io.save(new_struct_path)
            #shutil.copy(item, new_struct_path)
            new_structures.append(new_struct_path)
        os.chdir(wd)

        ## align all structures
        cmd = ["python",style["pymol_align"],ref_struct_path]+new_structures
        #print(cmd)
        p = subprocess.call(cmd)

        if not debug_mode:
            shutil.rmtree(new_folder)

        return new_structures
    else:
        print("#### subprocess error")
        os.chdir(wd)



def write_monomer_script(pir_file, ref_knowns, query_name, script_file, n_models):
    script_string = """from modeller import *
from modeller.automodel import *

log.verbose()

env = environ()
env.io.atom_files_directory = ['.', '../atom_files']
env.io.hetatm=True
a = automodel(env, alnfile='{pir_file}',
              knowns='{ref_knowns}', sequence='{query_name}',assess_methods=(assess.DOPE, assess.GA341))
a.final_malign3d=True
a.starting_model = 1
a.ending_model = {n_models}
a.make()
    """
    if os.path.exists(script_file):
        raise ValueError("Script file already exists: {0}".format(script_file))
    with open(script_file,"w") as f:
        f.write(script_string.format(pir_file=pir_file, ref_knowns=ref_knowns,
                query_name=query_name, n_models=n_models))
    return




#### homo dimer model for the sake of constraints
def built_dimer_model(ref_struct_path, query_seq, n_models, new_query_name="query",
        pir_file="aln.pir", script_file="model.py", chain_id='A',
        new_folder="build_model",
        keep_both=False,
        debug_mode=False):
    '''
    This works only for homodimers, where one only uses the constraints from
    the homodimer as constraints for the conformational flexibility, e.g.
    constraints for the flexibility of the N-terminal loop of hedgehogs.
    The sequence has to be the same for both chains, of course!

    Args:
        ref_struct_path -> path to the structure that contains chain A+B with
                        the correct sequence, which is extracted during the
                        calculation. For every non standard amino acid, a
                        '.' is used.
        query_seq -> string, which contains the modelled sequence
        new_query_name -> names the final structures, i.e. "query_name_1.pdb"

    Returns:
        A list with the new structure paths.
    '''
    wd = os.getcwd()
    if not os.path.exists(os.path.split(ref_struct_path)[-1]):
        shutil.copy(ref_struct_path, os.path.split(ref_struct_path)[-1])

    # make sure, it is not an absolute path
    ref_struct_path = os.path.split(ref_struct_path)[-1]

    # check if the pdb path is already in the working directory
    parser = PDBParser()
    pdb = parser.get_structure("dummy", ref_struct_path)
    if len(pdb) > 1:
        raise ValueError("Can not handle structures with more than one MODEL!\nThis structure has {0} MODELS!".format(len(pdb)))

    ref_seq = ""
    non_standard_seq = ""
    res_ids = []
    for res in pdb[0][chain_id]:
        res_ids.append(res.id[1])
        if res.resname in PDBFile.three2oneletter:
            ref_seq = "{0}{1}".format(ref_seq, PDBFile.three2oneletter[res.resname])
        else:
            # store the non standard amino acids in a variable
            non_standard_seq = "{0}.".format(non_standard_seq)

    start_res_id = min(res_ids)

    aln_result = get_pairwise_alignment(ref_seq, query_seq)
    aln_ref = aln_result["ref_seq"]
    aln_query = aln_result["query_seq"]

    # add non standard amino acids to both string
    aln_ref = "{0}{1}".format(aln_ref, non_standard_seq)
    aln_query = "{0}{1}".format(aln_query, non_standard_seq)

    if os.path.exists(new_folder):
        raise ValueError("Folder already exists: {0}".format(new_folder))

    os.mkdir(new_folder)
    os.chdir(new_folder)
    shutil.copy(os.path.join(wd, ref_struct_path), ref_struct_path)

    aln_ref = "".join(aln_ref)
    aln_query = "".join(aln_query)

    ref_knowns = "ref"
    query_name = "query"

    def write_pir(ref_seq,ref_struct_path, query_seq, pir_file):
        pir_string = """>P1;ref
structureX:{struct_path}::A::B:::0.00: 0.00
{ref_seq}/\n
{ref_seq}
*
>P1;query
sequence:query::A::B:::0.00: 0.00
{query_seq}/\n
{query_seq}
*"""

        if os.path.exists(pir_file):
            raise ValueError("Alingment file already exists: {0}".format(pir_file))
        with open(pir_file,"w") as f:
            f.write(pir_string.format(struct_path=ref_struct_path, ref_seq=ref_seq,
                    query_seq=query_seq))
        return


    write_pir(aln_ref, ref_struct_path, aln_query,pir_file)
    write_dimer_script(pir_file, ref_knowns, query_name, script_file, n_models)
    subprocess.call(["python", script_file])

    # clean up
    file_list = os.listdir('.')
    file_list = filter(lambda x:x.startswith(query_name) and x.endswith("_fit.pdb"), file_list)
    file_list.sort()
    new_structures = []
    for i,item in enumerate(file_list):
        new_struct_path = os.path.join(wd,"{0}_{1}.pdb".format(new_query_name, i))
        parser = PDBParser()
        new_pdb = parser.get_structure("dummy", item)
        if keep_both:
            pass
        else:
            new_pdb[0].detach_child("B")
        ## debug
        for res in new_pdb[0]["A"].child_list:
            old_id = list(res.id)
            old_id[1] += start_res_id - 1
            res.id = tuple(old_id)
        io = PDBIO()
        io.set_structure(new_pdb)
        io.save(new_struct_path)
        new_structures.append(new_struct_path)

    os.chdir(wd)

    ## align all structures
    cmd = ["python",style["pymol_align"],ref_struct_path]+new_structures
    #print(cmd)
    p = subprocess.call(cmd)

    if not debug_mode:
        shutil.rmtree(new_folder)

    return new_structures

def write_dimer_script(pir_file, ref_knowns, query_name, script_file, n_models):
    script_string = """from modeller import *
from modeller.automodel import *

log.verbose()

class MyModel(automodel):
    def special_restraints(self, aln):
        # Constrain the A and B chains to be identical (but only restrain
        # the C-alpha atoms, to reduce the number of interatomic distances
        # that need to be calculated):
        s1 = selection(self.chains['A']).only_atom_types('CA')
        s2 = selection(self.chains['B']).only_atom_types('CA')
        self.restraints.symmetry.append(symmetry(s1, s2, 1.0))
    def user_after_single_model(self):
        # Report on symmetry violations greater than 1A after building
        # each model:
        self.restraints.symmetry.report(1.0)

env = environ()
env.io.atom_files_directory = ['.', '../atom_files']
env.io.hetatm=True
a = MyModel(env, alnfile='{pir_file}',
              knowns='{ref_knowns}', sequence='{query_name}',assess_methods=(assess.DOPE, assess.GA341))
a.final_malign3d=True
a.starting_model = 1
a.ending_model = {n_models}
a.make()
    """
    if os.path.exists(script_file):
        raise ValueError("Script file already exists: {0}".format(script_file))
    with open(script_file,"w") as f:
        f.write(script_string.format(pir_file=pir_file, ref_knowns=ref_knowns,
                query_name=query_name, n_models=n_models))
    return

