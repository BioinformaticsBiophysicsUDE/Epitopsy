
:mod:`CalculateInteractionEnergyCython` --- Compute protein-ligand free binding energy
======================================================================================

.. module:: CalculateInteractionEnergyCython
   :synopsis: Compute protein-ligand free binding energy.
.. moduleauthor:: Christoph Wilms <christoph.wilms@uni-due.de>
.. sectionauthor:: Jean-Noel Grad <jean-noel.grad@uni-due.de>


This module provides tools to compute the Gibbs free energy of binding
between a ligand and a protein. It is meant to be used by
:func:`calculate_partition_function`.

.. _CalculateInteractionEnergyCython-syntax:

Module Syntax
-------------

Empty

.. _contents-of-module-CalculateInteractionEnergyCython:

Module Contents
---------------

.. function:: score_energy(potbox, ligand_coords, charges)

    Multiply the atomic charge by the magnitude of the protein electrostatic
    potential *V* for each atom in the molecular probe. Sum over all atoms.
    
    :math:`E = \sum_n^N q_n \left( i,j,k \right) \cdot V \left( i,j,k \right)`
    
    with i,j,k the coordinates on the grid.
    
    :param espbox: electrostatic box from an APBS calculation
    :type  espbox: np.ndarray
    :param ligand_coords: atomic positions of the molecular probe, rotated and
        snapped to the grid
    :type  ligand_coords: np.ndarray
    :param charges: all atomic charges from the PQR file
    :type  charges: list
    
    :returns: The sum of the electrostatic energy between the ligand and the
        protein, adimensional.


.. function:: transform_box_to_real_space(grid_coord, box_mesh_size, box_offset)

    Calculate the X,Y,Z coordinates in Angstroms of a particular grid point
    whose indices in the Numpy matrix are given by **grid_coord**.
    
    :math:`\displaystyle \text{new\_grid\_coord}[i] = \text{grid\_coord}[i] \cdot \text{box\_mesh\_size}[i] + \text{box\_offset}[i]`
    
    See :func:`transform_real_to_box_space` for the reverse operation.
    
    :param grid_coord: indices of a grid point from the Numpy matrix
        :attr:`DXFile.DXBox.box`
    :type  grid_coord: list
    :param box_mesh_size: grid resolution in Angstroms, as stored in
        :attr:`DXFile.DXBox.box_mesh_size`
    :type  box_mesh_size: list
    :param box_offset: vector between the grid origin and the DXBox geometrical
        center in Angstroms, as stored in :attr:`DXFile.DXBox.box_offset`
    :type  box_offset: list
    
    :returns: (X,Y,Z) coordinates of a grid point in Angstroms.


.. function:: transform_real_to_box_space(grid_coord, box_mesh_size, box_offset)

    Do the opposite of :func:`transform_box_to_real_space`: from X,Y,Z atomic
    coordinates in Angstroms, find the closest grid point.

    :math:`\displaystyle \text{new\_grid\_coord}[i] = \frac{\text{grid\_coord}[i] - \text{box\_offset}[i]}{\text{box\_mesh\_size}[i]}`
    
    :param grid_coord: X,Y,Z atomic coordinates in Angstroms
    :type  grid_coord: list
    :param box_mesh_size: grid resolution in Angstroms, as stored in
        :attr:`DXFile.DXBox.box_mesh_size`
    :type  box_mesh_size: list
    :param box_offset: vector between the grid origin and the DXBox geometrical
        center in Angstroms, as stored in :attr:`DXFile.DXBox.box_offset`
    :type  box_offset: list
    
    :returns: Grid point indices for a single atom.


.. function:: move_ligand(ligand_real_coords, coord)

    Compute the geometrical center of the molecular probe and move it to the
    coordinates of the grid point **coord**.
    
    :param ligand_real_coords: coordinates of every atoms in the molecular probe
    :type  ligand_real_coords: list
    :param coord: coordinates of the grid point, as returned by
        :func:`transform_box_to_real_space`
    :type  coord: list
    
    :returns: The coordinates of the molecular probe whose geometrical center
        coincides with the grid point given as parameter.


.. function:: transform_ligand_to_box(ligand_real_coords, box_mesh_size, box_offset)

    Snap each atom in the molecular probe to the closest grid points. Make use
    of :func:`transform_real_to_box_space`.
    
    :param ligand_real_coords: atomic coordinates of the molecular probe
    :type  ligand_real_coords: np.ndarray
    :param box_mesh_size: grid resolution in Angstroms
    :type  box_mesh_size: list
    :param box_offset: vector between the grid origin and the DXBox geometrical
        center in Angstroms
    :type  box_offset: list
    
    :returns: A list of grid points (discretized molecular probe).


.. function:: out_of_box(ligand_coords, box_dim)

    Check if the molecular probe has left the box.
    
    :param ligand_coords: atomic positions of a discretized molecular probe,
        as obtained from :func:`transform_ligand_to_box`
    :type  ligand_coords: np.ndarray
    :param box_dim: number of grid points in the X,Y,Z directions
    :type  box_dim: list
    
    :returns: ``True`` if the molecular probe left the box, ``False`` otherwise


.. function:: is_overlapp(vdwbox, ligand_coords, vdw_radii, solvent_score, protein_score)

    Check if the molecular probe overlaps with the protein.
    
    Method: for every atom in the ligand, find all grid points within the van
    der Waals radius of that atom and return ``True`` if any of them is equal
    to **protein_score**
    
    :param vdwbox: van der Waals box from an APBS calculation
    :type  vdwbox: np.ndarray
    :param ligand_coords: atomic positions of a discretized molecular probe,
        as obtained from :func:`transform_ligand_to_box`
    :type  ligand_coords: np.array
    :param vdw_radii: atomic radii snapped to the grid (*i.e.* the quantity of
        grid point that can be aligned in the X,Y, and Z directions
        respectively, typically between 1 and 6)
    :type  vdw_radii: np.ndarray
    :param solvent_score: score for the solvent, usually 1.0
    :type  solvent_score: float
    :param protein_score: score for the protein, usually 0.0
    :type  protein_score: float
    
    :returns: ``True`` if there is a least one grid point overlapping,
        ``False`` otherwise


.. automodule:: epitopsy.cython.CalculateInteractionEnergyCython
    :members:
    :undoc-members:
    :show-inheritance:

