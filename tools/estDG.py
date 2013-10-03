import numpy as np

from epitopsy.DXFile import VDWBox

# conc = 0.01-0.1 mol/l

def findVolLASPro(counter_box):

    # read the counter matrix
    #counter_matrix = read_dxfile(counter_matrix_path,'esp')

    # find the maximum of performed rotations in this matrix and every point
    # which has the this value in the matrix
    max_v = np.amax(counter_box.box)
#    solv_pos = np.nonzero(esp.box == max_v)
    # find the real surface of the protein
    solv_pos = np.nonzero(counter_box.box > 0)

    # create a new box of the same dimension of esp and write a 1 everywhere
    # where the rotation of the maximum value was found
    new_box = np.zeros(counter_box.box.shape)
    new_box[solv_pos] = 1

    # create a vdw object and find the surface indices
    vdw = VDWBox(new_box,counter_box.box_mesh_size,counter_box.box_offset)
    vdw.flood()
    vdw.find_solvent_surface()

    # find solvent, LAS , Protein

    solv = np.nonzero(vdw.box == 1)
    LAS = np.nonzero(vdw.box == 2)
    prot = np.nonzero(vdw.box == 0)

    return solv, LAS, prot

def calcDG(energy_box,counter_box,Temp=310,conc=None):
    '''
    Calculates the binding free energy of protein to ligand in kJ/mol and
    the k_D.

    Args:
        energy_box -> DXBox containing the energies
        counter_box -> DXBox containing the number of possible rotations at
                        each grid point
        Temp -> Temperature, default 310
        conc -> concentration

    Returns:
        A list containing the binding free energy of protein to ligand and
        the k_D.
    '''
    # gas constant
    R = 8.3144621

    # get solvent, LAS and protein indices
    SoLAPr = findVolLASPro(counter_box)

    # read energy matrix
    #EnMat = read_dxfile(energy_matrix_path, 'esp')

    # calculate e^(-E/(k_B T)) and sum up
    # energies are already in units of k_B T
    #solvAr = EnMat.box[SoLAPr[0]]
    #solvE = np.sum(np.exp(-solvAr))

    LASAr = energy_box.box[SoLAPr[1]]
    LASE = np.sum(np.exp(-1*LASAr))

    if conc is None:
        void = len(LASAr)
    else:
        void = calVoidVol(conc,SoLAPr[1],SoLAPr[2])

    #print(str(calcConc(void,SoLAPr[1],SoLAPr[2]))+' mol/l')

    # return DG estimate
    # DG = - R T ln( K )
    # K = sum( E_LAS ) / sum( E_nonLAS )
    # E_LAS = exp( -E_{onLAS} / (k_B T) )
    # E_nonLAS = exp( -E_{notonLAS} / (k_B T) )
    return [(-1*R*Temp*(np.log(LASE/void))/1000), LASE]

def calcConc(void,LASpoints,ProteinPoints,grid=[0.5,0.5,0.5]):

    '''
    calculates the concentration of protein and ligand for a given LAS layer
    and and the additional volume of the protein and a given volume which
    constitutes the rest of the water in which all is soluted

    all has to be given in grid points, the grid defines afterwards the volume
    of each gridpoint
    '''

    mol = 6.02214129 * 10**23
    in_liter = 10**-27
    water = ((len(LASpoints)*grid[0]*grid[1]*grid[2] * in_liter)
            + (void*grid[0]*grid[1]*grid[2] * in_liter))
    protein = len(ProteinPoints)*grid[0]*grid[1]*grid[2] * in_liter
    volume = protein+water
    conc = 1/(mol*volume)
    return conc

def calVoidVol(conc,LASpoints,ProteinPoints,grid=[0.5,0.5,0.5]):

    '''
    calculates a void volume, from a given LAS-layer, a protein volume (both in
    grid points) and a concentration in which the protein should be soluted'''

    mol = 6.02214129 * 10**23
    in_liter = 10**-27
    vol = 1.0/(mol*float(conc))
    waterLAS = len(LASpoints)*grid[0]*grid[1]*grid[2] * in_liter
    protein = len(ProteinPoints)*grid[0]*grid[1]*grid[2]* in_liter
    voidVol = vol-waterLAS-protein
    voidpoints = voidVol/(grid[0]*grid[1]*grid[2] * in_liter)
    return voidpoints
