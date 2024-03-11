"""
from MASCLET framework
https://github.com/dvallesp/masclet_framework

parameters module
Functions for reading and creating parameters JSON files

Created by David Vallés
"""

#  Last update on 17/3/20 19:42

# GENERAL PURPOSE AND SPECIFIC LIBRARIES USED IN THIS MODULE

import json
import os

# FUNCTIONS DEFINED IN THIS MODULE
def read_parameters_file(filename='vortex_parameters.json', path=''):
    """
    Returns dictionary containing the parameters file, that have been previously written with the
    write_parameters() function in this same module.

    Args:
        filename: name of the vortex JSON parameters file (str)
        path: path to the file (str)

    Returns:
        dictionary containing the parameters (and their names), namely:
        NMAX, NMAY, NMAZ: number of l=0 cells along each direction (int)
        NPALEV: maximum number of refinement cells per level (int)
        NLEVELS: maximum number of refinement level (int)
        NAMRX (NAMRY, NAMRZ): maximum size of refinement patches (in l-1 cell units) (int)
        SIZE: side of the simulation box in the chosen units (typically Mpc or kpc) (float)

    """
    filepath = os.path.join(path, filename)
    with open(filepath) as json_file:
        data = json.load(json_file)
    return data


def read_parameters(filename='vortex_parameters.json', path='', load_nma=True, load_nlevels=True,
                    load_size=True):
    """
    Returns MASCLET parameters in the old-fashioned way (as a tuple).

    Args:
        filename: name of the vortex JSON parameters file (str)
        path: path to the file (str)
        load_nma: whether NMAX, NMAY, NMAZ are read (bool)
        load_nlevels: whether NLEVELS is read (bool)
        load_size: whether SIZE is read (bool)

    Returns:
        tuple containing, in this exact order, the chosen parameters from:
        NMAX, NMAY, NMAZ: number of l=0 cells along each direction (int)
        NPALEV: maximum number of refinement cells per level (int)
        NLEVELS: maximum number of refinement level (int)
        NAMRX (NAMRY, NAMRZ): maximum size of refinement patches (in l-1 cell units) (int)
        SIZE: side of the simulation box in the chosen units (typically Mpc or kpc) (float)

    """
    parameters = read_parameters_file(filename=filename, path=path)
    returnvariables = []
    if load_nma:
        returnvariables.extend([parameters[i] for i in ['NMAX', 'NMAY', 'NMAZ']])
    if load_nlevels:
        returnvariables.append(parameters['NLEVELS'])
    if load_size:
        returnvariables.append(parameters['SIZE'])
    return tuple(returnvariables)


def write_parameters(nmax, nmay, nmaz, nlevels, size, 
                     filename='vortex_parameters.json', path=''):
    """
    Creates a JSON file containing the parameters of a certain simulation

    Args:
        nmax: number of l=0 cells along the X-direction (int)
        nmay: number of l=0 cells along the Y-direction (int)
        nmaz: number of l=0 cells along the Z-direction (int)
        nlevels: maximum number of refinement level (int)
        size: side of the simulation box in the chosen units (typically Mpc or kpc) (float)
        filename: name of the vortex JSON parameters file to be saved (str)
        path: path to the file to be saved (str)

    Returns: nothing. A file is created in the specified path
    """
    parameters = {'NMAX': nmax, 'NMAY': nmay, 'NMAZ': nmaz,
                  'NLEVELS': nlevels,
                  'SIZE': size}

    with open(os.path.join(path,filename), 'w') as json_file:
        json.dump(parameters, json_file)
