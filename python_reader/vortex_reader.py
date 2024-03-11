import numpy as np
import os
import fortran_reader as fr
import parameters

###########################################################
### Grid building and particle-to-grid assignments
###########################################################

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
        cr0amr: at all patches and levels, =0 if the patch is refined 
                by a higher resolution level, =1 otherwise (list of
                numpy arrays)
        solap: at all patches and levels, =0 if there is another patch 
               at the same level overlapping with the current patch, =1
               otherwise (list of numpy arrays)
        
    """

    if parameters_path is None:
        parameters_path = path
    if grids_path is None:
        grids_path = path
    
    nmax, nmay, nmaz, size = parameters.read_parameters(load_nma=True, load_nlevels=False,
                                                        load_size=True, path=parameters_path)

    grids = read_grids(it, path=grids_path, parameters_path=parameters_path,
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


def read_gridded_density(it, path='', parameters_path=None, filename='gridded_density',
                       digits=5, grids_path=None, grids_filename='grids'):
    """
    Reads the gridded_densityXXXXX file, containing the (unnormalised)
     density field according to the same grid assignment procedure in
     vortex-p.

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
        density: at all patches and levels, the density field (list of
                numpy arrays), in arbitrary units.
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
        density = [np.reshape(fr.read_record(f, 'f4'), (nmax,nmax,nmax), 'F')]
        for i in range(1,npatch.sum()+1):
            density.append(np.reshape(fr.read_record(f, 'f4'), (patchnx[i],patchny[i],patchnz[i]), 'F'))

    return density


def read_gridded_kernel_length(it, path='', parameters_path=None, 
                               filename='gridded_kernel_length',
                               digits=5, grids_path=None, grids_filename='grids'):
    """
    Reads the gridded_kernel_lengthXXXXX file, which contains, in each 
     location, the kernel length used to assign the velocity field.

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
        kernel_length: at all patches and levels, the kernel length
                       numpy arrays), in input length units.
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
        kernel_length = [np.reshape(fr.read_record(f, 'f4'), (nmax,nmax,nmax), 'F')]
        for i in range(1,npatch.sum()+1):
            kernel_length.append(np.reshape(fr.read_record(f, 'f4'), (patchnx[i],patchny[i],patchnz[i]), 'F'))

    return kernel_length


def read_gridded_velocity(it, path='', parameters_path=None, 
                          filename='gridded_velocity',
                          digits=5, grids_path=None, grids_filename='grids'):
    """
    Reads the gridded_velocityXXXXX file, which contains, the gridded 
     version of the input velocity field

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
        vx, vy, vz: input velocity field, prior to any decomposition
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
        vx = [np.reshape(fr.read_record(f, 'f4'), (nmax,nmax,nmax), 'F')]
        vy = [np.reshape(fr.read_record(f, 'f4'), (nmax,nmax,nmax), 'F')]
        vz = [np.reshape(fr.read_record(f, 'f4'), (nmax,nmax,nmax), 'F')]
        for i in range(1,npatch.sum()+1):
            vx.append(np.reshape(fr.read_record(f, 'f4'), (patchnx[i],patchny[i],patchnz[i]), 'F'))
            vy.append(np.reshape(fr.read_record(f, 'f4'), (patchnx[i],patchny[i],patchnz[i]), 'F'))
            vz.append(np.reshape(fr.read_record(f, 'f4'), (patchnx[i],patchny[i],patchnz[i]), 'F'))

    return vx, vy, vz


def read_gridded_mach(it, path='', parameters_path=None, 
                          filename='gridded_mach',
                          digits=5, grids_path=None, grids_filename='grids'):
    """
    Reads the gridded_machXXXXX file, which contains, the gridded 
     version of the input Mach field (only used for the multiscale filter)

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
        mach: gridded input Mach number field
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
        mach = [np.reshape(fr.read_record(f, 'f4'), (nmax,nmax,nmax), 'F')]
        for i in range(1,npatch.sum()+1):
            mach.append(np.reshape(fr.read_record(f, 'f4'), (patchnx[i],patchny[i],patchnz[i]), 'F'))


    return mach


###########################################################
### Grid results
###########################################################

def read_vcomp(it, path='', parameters_path=None, 
               filename='vcomp',
               digits=5, grids_path=None, grids_filename='grids'):
    """
    Reads the vcompXXXXX file, which contains the compressive part of 
     the velocity field.

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
        vcompx, vcompy, vcompz: gridded compressive velocity field
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
        vcompx = [np.reshape(fr.read_record(f, 'f4'), (nmax,nmax,nmax), 'F')]
        vcompy = [np.reshape(fr.read_record(f, 'f4'), (nmax,nmax,nmax), 'F')]
        vcompz = [np.reshape(fr.read_record(f, 'f4'), (nmax,nmax,nmax), 'F')]
        for i in range(1,npatch.sum()+1):
            vcompx.append(np.reshape(fr.read_record(f, 'f4'), (patchnx[i],patchny[i],patchnz[i]), 'F'))
            vcompy.append(np.reshape(fr.read_record(f, 'f4'), (patchnx[i],patchny[i],patchnz[i]), 'F'))
            vcompz.append(np.reshape(fr.read_record(f, 'f4'), (patchnx[i],patchny[i],patchnz[i]), 'F'))

    return vcompx, vcompy, vcompz


def read_vsol(it, path='', parameters_path=None, 
               filename='vcomp',
               digits=5, grids_path=None, grids_filename='grids'):
    """
    Reads the vsolXXXXX file, which contains the solenoidal part of 
     the velocity field.

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
        vsolx, vsoly, vsolz: gridded compressive velocity field
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
        vsolx = [np.reshape(fr.read_record(f, 'f4'), (nmax,nmax,nmax), 'F')]
        vsoly = [np.reshape(fr.read_record(f, 'f4'), (nmax,nmax,nmax), 'F')]
        vsolz = [np.reshape(fr.read_record(f, 'f4'), (nmax,nmax,nmax), 'F')]
        for i in range(1,npatch.sum()+1):
            vsolx.append(np.reshape(fr.read_record(f, 'f4'), (patchnx[i],patchny[i],patchnz[i]), 'F'))
            vsoly.append(np.reshape(fr.read_record(f, 'f4'), (patchnx[i],patchny[i],patchnz[i]), 'F'))
            vsolz.append(np.reshape(fr.read_record(f, 'f4'), (patchnx[i],patchny[i],patchnz[i]), 'F'))

    return vsolx, vsoly, vsolz


def read_spot(it, path='', parameters_path=None, 
              filename='spot',
              digits=5, grids_path=None, grids_filename='grids'):
    """
    Reads the spotXXXXX file, which contains the scalar potential of 
     the Helmholtz-Hodge decomposition, whose gradient yields the
     compressive component.

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
        spot: gridded scalar potential
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
        spot = [np.reshape(fr.read_record(f, 'f4'), (nmax,nmax,nmax), 'F')]
        for i in range(1,npatch.sum()+1):
            spot.append(np.reshape(fr.read_record(f, 'f4'), (patchnx[i],patchny[i],patchnz[i]), 'F'))

    return spot


def read_vpot(it, path='', parameters_path=None, 
            filename='vpot',
            digits=5, grids_path=None, grids_filename='grids'):
    """
    Reads the vpotXXXXX file, which contains the vector potential of 
     the Helmholtz-Hodge decomposition, whose curl yields the
     solenoidal component.

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
        vsolx, vsoly, vsolz: gridded compressive velocity field
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
        vpotx = [np.reshape(fr.read_record(f, 'f4'), (nmax,nmax,nmax), 'F')]
        vpoty = [np.reshape(fr.read_record(f, 'f4'), (nmax,nmax,nmax), 'F')]
        vpotz = [np.reshape(fr.read_record(f, 'f4'), (nmax,nmax,nmax), 'F')]
        for i in range(1,npatch.sum()+1):
            vpotx.append(np.reshape(fr.read_record(f, 'f4'), (patchnx[i],patchny[i],patchnz[i]), 'F'))
            vpoty.append(np.reshape(fr.read_record(f, 'f4'), (patchnx[i],patchny[i],patchnz[i]), 'F'))
            vpotz.append(np.reshape(fr.read_record(f, 'f4'), (patchnx[i],patchny[i],patchnz[i]), 'F'))

    return vpotx, vpoty, vpotz


def read_divv(it, path='', parameters_path=None, 
              filename='divv',
              digits=5, grids_path=None, grids_filename='grids'):
    """
    Reads the divvXXXXX file, which contains the divergence of the
     velocity field.

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
        divv: gridded div(v) field
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
        divv = [np.reshape(fr.read_record(f, 'f4'), (nmax,nmax,nmax), 'F')]
        for i in range(1,npatch.sum()+1):
            divv.append(np.reshape(fr.read_record(f, 'f4'), (patchnx[i],patchny[i],patchnz[i]), 'F'))

    return divv


def read_curlv(it, path='', parameters_path=None, 
            filename='curlv',
            digits=5, grids_path=None, grids_filename='grids'):
    """
    Reads the curlvXXXXX file, which contains the curl of the velocity 
     field.

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
        curlvx, curlvy, curlvz: gridded curl(v) field
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
        curlvx = [np.reshape(fr.read_record(f, 'f4'), (nmax,nmax,nmax), 'F')]
        curlvy = [np.reshape(fr.read_record(f, 'f4'), (nmax,nmax,nmax), 'F')]
        curlvz = [np.reshape(fr.read_record(f, 'f4'), (nmax,nmax,nmax), 'F')]
        for i in range(1,npatch.sum()+1):
            curlvx.append(np.reshape(fr.read_record(f, 'f4'), (patchnx[i],patchny[i],patchnz[i]), 'F'))
            curlvy.append(np.reshape(fr.read_record(f, 'f4'), (patchnx[i],patchny[i],patchnz[i]), 'F'))
            curlvz.append(np.reshape(fr.read_record(f, 'f4'), (patchnx[i],patchny[i],patchnz[i]), 'F'))

    return curlvx, curlvy, curlvz


###########################################################
### Filter results
###########################################################

def read_gridded_filtlen(it, path='', parameters_path=None, 
                         filename='gridded_filtlen',
                         digits=5, grids_path=None, grids_filename='grids'):
    """
    Reads the gridded_filtlenXXXXX file, which contains the filter  
     length used in the multiscale filter.

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
        filtlen: gridded filter length field
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
        filtlen = [np.reshape(fr.read_record(f, 'f4'), (nmax,nmax,nmax), 'F')]
        for i in range(1,npatch.sum()+1):
            filtlen.append(np.reshape(fr.read_record(f, 'f4'), (patchnx[i],patchny[i],patchnz[i]), 'F'))

    return filtlen


def read_gridded_vturb(it, path='', parameters_path=None, 
                       filename='gridded_vturb',
                       digits=5, grids_path=None, grids_filename='grids'):
    """
    Reads the gridded_vturbXXXXX file, which contains the turbulent  
     velocity field, as computed by the multiscale filter.

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
        vturbx, vturby, vturbz: gridded turbulent velocity field
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
        vturbx = [np.reshape(fr.read_record(f, 'f4'), (nmax,nmax,nmax), 'F')]
        vturby = [np.reshape(fr.read_record(f, 'f4'), (nmax,nmax,nmax), 'F')]
        vturbz = [np.reshape(fr.read_record(f, 'f4'), (nmax,nmax,nmax), 'F')]
        for i in range(1,npatch.sum()+1):
            vturbx.append(np.reshape(fr.read_record(f, 'f4'), (patchnx[i],patchny[i],patchnz[i]), 'F'))
            vturby.append(np.reshape(fr.read_record(f, 'f4'), (patchnx[i],patchny[i],patchnz[i]), 'F'))
            vturbz.append(np.reshape(fr.read_record(f, 'f4'), (patchnx[i],patchny[i],patchnz[i]), 'F'))

    return vturbx, vturby, vturbz


def read_shocked(it, path='', parameters_path=None, 
                 filename='shocked',
                 digits=5, grids_path=None, grids_filename='grids'):
    """
    Reads the shockedXXXXX file, which marks the shocked regions 
     used by the multiscale filter as an stopping criterion.

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
        shocked: gridded shocked regions
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
        shocked = [np.reshape(fr.read_record(f, 'i1'), (nmax,nmax,nmax), 'F')]
        for i in range(1,npatch.sum()+1):
            shocked.append(np.reshape(fr.read_record(f, 'i1'), (patchnx[i],patchny[i],patchnz[i]), 'F'))

    return shocked


###########################################################
### Particle-wise results
###########################################################

def read_errorparticles(it, path='', parameters_path=None,
                         filename='error-particles',  digits=5):
    """
    Reads the error-particlesXXXXX file, which contains the error
     (difference between the input velocity and the smoothed velocity)
     in the velocity assignment.

    Args:
        it: iteration number (int)
        path: path to the error-particles file (str)
        parameters_path: path to the parameters file (str). If None, 
         the parameters file is assumed to be in the same directory as
         the error-particles file.
        filename: name of the error-particles file (str)
        digits: number of digits in the iteration number (int)

    Returns:
        error: error in the velocity assignment
    """

    filename = filename + str(it).zfill(digits)

    with open(os.path.join(path, filename), 'rb') as f:
        error = fr.read_record(f, 'f4')

    return error


def read_velocityparticles(it, path='', parameters_path=None,
                           filename='velocity-particles',  digits=5,
                           return_original_position=False,
                           return_original_velocity=False,
                           return_original_mass=False,
                           return_smoothed_velocity=True,
                           return_smoothed_compressive_velocity=True,
                           return_smoothed_solenoidal_velocity=True):
    """
    Reads the velocity-particlesXXXXX file, which contains the
     decomposition results reinterpolated to the particle positions.

    Args:
        it: iteration number (int)
        path: path to the velocity-particles file (str)
        parameters_path: path to the parameters file (str). If None, 
         the parameters file is assumed to be in the same directory as
         the velocity-particles file.
        filename: name of the velocity-particles file (str)
        digits: number of digits in the iteration number (int)

    Returns:

    """

    filename = filename + str(it).zfill(digits)

    results = {}

    with open(os.path.join(path, filename), 'rb') as f:
        npart = fr.read_record(f, 'i4')[0]
        
        if return_original_position:
            results['orig_xpart'] = fr.read_record(f, 'f4')
            results['orig_ypart'] = fr.read_record(f, 'f4')
            results['orig_zpart'] = fr.read_record(f, 'f4')

        else:
            for i in range(3):
                fr.skip_record(f)

        if return_original_velocity:
            results['orig_vxpart'] = fr.read_record(f, 'f4')
            results['orig_vypart'] = fr.read_record(f, 'f4')
            results['orig_vzpart'] = fr.read_record(f, 'f4')
        else:
            for i in range(3):
                fr.skip_record(f)

        if return_original_mass:
            results['orig_mpart'] = fr.read_record(f, 'f4')
        else:
            fr.skip_record(f)

        if return_smoothed_velocity:
            results['smooth_vxpart'] = fr.read_record(f, 'f4')
            results['smooth_vypart'] = fr.read_record(f, 'f4')
            results['smooth_vzpart'] = fr.read_record(f, 'f4')
        else:
            for i in range(3):
                fr.skip_record(f)

        if return_smoothed_compressive_velocity:
            results['smooth_vxcomppart'] = fr.read_record(f, 'f4')
            results['smooth_vycomppart'] = fr.read_record(f, 'f4')
            results['smooth_vzcomppart'] = fr.read_record(f, 'f4')
        else:
            for i in range(3):
                fr.skip_record(f)

        if return_smoothed_solenoidal_velocity:
            results['smooth_vxsolpart'] = fr.read_record(f, 'f4')
            results['smooth_vysolpart'] = fr.read_record(f, 'f4')
            results['smooth_vzsolpart'] = fr.read_record(f, 'f4')
        else:
            for i in range(3):
                fr.skip_record(f)

    return results
    
