import numpy as np
import os
import fortran_reader as fr
import parameters

def read_grids(it, path='', parameters_path=None, filename='grids', 
                digits=5, nparray=True, output_dictionary=True):
    """
    Reads the gridsXXXXX file, containing the description of the AMR 
     grid structure.

    Args:
        it: iteration number (int)
        path: path to the grids file (str)
        parameters_path: path to the parameters file (str). If None, 
         the parameters file is assumed to be in the same directory as
         the grids file.
        filename: name of the grids file (str)
        digits: number of digits in the iteration number (int)
        nparray: whether the output should be numpy arrays (bool)
        output_dictionary: whether the output should be a dictionary
            (bool)
        
    Returns:
        A dictionary with the grid information, if output_dictionary is True.
        Else, a tuple with the grid information.
    """
    if parameters_path is None:
        parameters_path = path

    nmax, nmay, nmaz, size = parameters.read_parameters(load_nma=True, load_nlevels=False,
                                                         load_size=True, path=parameters_path)
    rx = - size / 2 + size / nmax

    filename = filename + str(it).zfill(digits)
    with open(os.path.join(path, filename), 'r') as grids:
        # first, we load some general parameters
        irr, t, nl, mass_dmpart, _ = tuple(float(i) for i in grids.readline().split())
        irr = int(irr)
        # assert (it == irr)
        nl = int(nl)
        zeta = float(grids.readline().split()[0])
        # l=0
        _, ndxyz, _ = tuple(float(i) for i in grids.readline().split())[0:3]
        ndxyz = int(ndxyz)

        # vectors where the data will be stored
        npatch = [0]  # number of patches per level, starting with l=0
        npart = [ndxyz]  # number of dm particles per level, starting with l=0
        patchnx = [nmax]
        patchny = [nmay]
        patchnz = [nmaz]
        patchx = [0]
        patchy = [0]
        patchz = [0]
        patchrx = [rx]
        patchry = [rx]
        patchrz = [rx]
        pare = [0]

        for ir in range(1, nl + 1):
            level, npatchtemp, nparttemp = tuple(int(i) for i in grids.readline().split())[0:3]
            npatch.append(npatchtemp)
            npart.append(nparttemp)

            # ignoring a blank line
            grids.readline()

            # loading all values
            for i in range(sum(npatch[0:ir]) + 1, sum(npatch[0:ir + 1]) + 1):
                this_nx, this_ny, this_nz = tuple(int(i) for i in grids.readline().split())
                this_x, this_y, this_z = tuple(int(i) for i in grids.readline().split())
                this_rx, this_ry, this_rz = tuple(float(i) for i in grids.readline().split())
                this_pare = int(grids.readline())
                patchnx.append(this_nx)
                patchny.append(this_ny)
                patchnz.append(this_nz)
                patchx.append(this_x - 1)
                patchy.append(this_y - 1)
                patchz.append(this_z - 1)
                patchrx.append(this_rx)
                patchry.append(this_ry)
                patchrz.append(this_rz)
                pare.append(this_pare)

        # converts everything into numpy arrays if nparray set to True
        if nparray:
            npatch = np.array(npatch)
            npart = np.array(npart)
            patchnx = np.array(patchnx)
            patchny = np.array(patchny)
            patchnz = np.array(patchnz)
            patchx = np.array(patchx)
            patchy = np.array(patchy)
            patchz = np.array(patchz)
            patchrx = np.array(patchrx)
            patchry = np.array(patchry)
            patchrz = np.array(patchrz)
            pare = np.array(pare)

    if output_dictionary:
        return {'npatch': npatch, 
                'patchnx': patchnx, 'patchny': patchny, 'patchnz': patchnz,
                'patchx': patchx, 'patchy': patchy, 'patchz': patchz,
                'patchrx': patchrx, 'patchry': patchry, 'patchrz': patchrz,
                'pare': pare}
    else:
        return npatch, patchnx, patchny, patchnz, patchx, patchy, patchz, \
            patchrx, patchry, patchrz, pare

def read_grid_overlaps(it, path='', parameters_path=None, filename='grid_overlaps',
                       digits=5, grids_path=None, grids_filename='grids'):
    """
    Reads the grid_overlapsXXXXX file, containing the information of 
     the overlaps between patches at different and at the same level.

    Args:
        it: iteration number (int)
        path: path to the grid_overlaps file (str)
        parameters_path: path to the parameters file (str). If None, 
         the parameters file is assumed to be in the same directory as
         the grid_overlaps file.
        filename: name of the grid_overlaps file (str)
        digits: number of digits in the iteration number (int)
        grids_path: path to the grids file (str). If None, the grids file
         is assumed to be in the same directory as the grid_overlaps file.
        grids_filename: name of the grids file (str)


    Returns:
        
    """

    if parameters_path is None:
        parameters_path = path
    
    nmax, nmay, nmaz, size = parameters.read_parameters(load_nma=True, load_nlevels=False,
                                                        load_size=True, path=parameters_path)

    grids = read_grids(it, path=path, parameters_path=parameters_path,
                       filename=grids_filename)
    npatch = grids['npatch']
    patchnx = grids['patchnx']
    patchny = grids['patchny']
    patchnz = grids['patchnz']

    filename = filename + str(it).zfill(digits)

    with open(os.path.join(path, filename), 'rb') as f:
        cr0amr = [np.reshape(fr.read_record(f, 'i4'), (nmax,nmax,nmax), 'F')]
        solap = [np.ones(cr0amr[0].shape, dtype='i4')]
        for i in range(1,npatch.sum()+1):
            cr0amr.append(np.reshape(fr.read_record(f, 'i4'), (patchnx[i],patchny[i],patchnz[i]), 'F'))
            solap.append(np.reshape(fr.read_record(f, 'i4'), (patchnx[i],patchny[i],patchnz[i]), 'F'))

    return cr0amr, solap