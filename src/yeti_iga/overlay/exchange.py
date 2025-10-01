"""
Module for import/export of data
"""

import json
import numpy as np

from . import Patch, ElasticMaterial

def extract_keys(d, keys):
    return [d[k] for k in keys]

def load_geomdl_json(filename):
    """
    Load a set of geomdl entities (volumes, surfaces, ...) and return a list
    of patchs.

    Parameters
    ----------
    filename : string
        path to the file to read

    Return
    ------
    patchs : list
        a list of patchs
    """


    with open(filename, "r", encoding="utf-8") as f:
        data = json.load(f)
        patchs = []
        # TODO handle exceptions
        if data['shape']['type'] == 'volume':
            patch_type = 'volume'
            for ipatch in range(data['shape']['count']):
                patchs.append(read_one_patch_json_geomdl(data['shape']['data'][ipatch],
                                                         patch_type))
        elif data['shape']['type'] == 'surface':
            patch_type = 'surface'
            for ipatch in range(data['shape']['count']):
                patchs.append(read_one_patch_json_geomdl(data['shape']['data'][ipatch],
                                                         patch_type))
        else:
            raise Exception(f"Shape data type {data['shape']['type']} is not handled")


        return patchs

def read_one_patch_json_geomdl(patch_dict, patch_type):
    """
    Read data for one patch from a dict contained in a json file generated
    by geomdl

    WARNING : does not take into accountpaâtchs with repeated inner knots

    Parameters
    ----------
    patch_dict : dict
        Dictionnary containing data for one patch
    patch_type : string
        String indicatng patch type : 'volume', ...

    Return
    ------
    patch : Patch
        patch data for yeti
    """

    assert patch_dict['type'] == 'spline'

    if patch_type == 'volume':
        # for a solid volume, dimension of parametric space is dimension of phys space
        dimension = patch_dict['dimension']

        degrees = extract_keys(patch_dict, ['degree_u', 'degree_v', 'degree_w'])
        knotvectors = extract_keys(patch_dict, ['knotvector_u', 'knotvector_v', 'knotvector_w'])
        sizes_cp = extract_keys(patch_dict, ['size_u', 'size_v', 'size_w'])

        cps = np.array(patch_dict['control_points']['points'])          #.reshape(*sizes_cp, 3)

        # Points are given with geomdl connectivity convention : v increasing, u increasing, w increasing
        connectivity = np.arange(np.prod(sizes_cp))
        connectivity = connectivity.reshape(sizes_cp[1], sizes_cp[0], sizes_cp[2], order="F")

        nb_spans = np.subtract(sizes_cp, degrees)

        yeti_connectivity = np.empty((np.prod(nb_spans), np.prod(np.array(degrees)+1)), dtype=int)
        spans = np.empty((np.prod(nb_spans), dimension), dtype=int)

        ispan = 0
        for i in range(nb_spans[0]):
            for j in range(nb_spans[1]):
                for k in range(nb_spans[2]):
                    icp = 0
                    for kk in reversed(range(degrees[2]+1)):
                        for jj in reversed(range(degrees[1]+1)):
                            for ii in reversed(range(degrees[0]+1)):
                                yeti_connectivity[ispan, icp] = connectivity[j+jj, i+ii, k+kk]
                                icp += 1
                    spans[ispan, :] = np.add(degrees, [i, j, k])
                    ispan += 1

        return Patch(element_type='U1',
                     degrees=degrees,
                     knot_vectors = knotvectors,
                     control_points=cps,
                     weights=np.ones(cps.shape[0]),
                     connectivity=yeti_connectivity,
                     spans=spans,
                     material=ElasticMaterial(1., 1.))

    elif patch_type == 'surface':
        dimension = patch_dict['dimension']
        if dimension != 3:
            raise Exception('2D solid geometry is not handled')
        degrees = extract_keys(patch_dict, ['degree_u', 'degree_v'])
        knotvectors = extract_keys(patch_dict, ['knotvector_u', 'knotvector_v'])
        sizes_cp = extract_keys(patch_dict, ['size_u', 'size_v'])

        cps = np.array(patch_dict['control_points']['points'])

        connectivity = np.arange(np.prod(sizes_cp))
        # TODO à vérifier
        connectivity = connectivity.reshape(sizes_cp[1], sizes_cp[0], order="F")

        nb_spans = np.subtract(sizes_cp, degrees)

        yeti_connectivity = np.empty((np.prod(nb_spans), np.prod(np.array(degrees)+1)), dtype=int)
        # spans = np.empty((np.prod(nb_spans), dimension))
        # For a surface, dimension of parametric space is 2
        spans = np.empty((np.prod(nb_spans), 2), dtype=int)

        ispan = 0
        for i in range(nb_spans[0]):
            for j in range(nb_spans[1]):
                    icp = 0
                    for jj in reversed(range(degrees[1]+1)):
                        for ii in reversed(range(degrees[0]+1)):
                            yeti_connectivity[ispan, icp] = connectivity[j+jj, i+ii]
                            icp += 1
                    spans[ispan, :] = np.add(degrees, [i, j])
                    ispan += 1

        return Patch(element_type='U3',
                     degrees=degrees,
                     knot_vectors = knotvectors,
                     control_points=cps,
                     weights=np.ones(cps.shape[0]),
                     connectivity=yeti_connectivity,
                     spans=spans,
                     material=ElasticMaterial(1., 1.),
                     properties=np.array([1.]))
