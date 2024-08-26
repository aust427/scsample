import h5py
import numpy as np
<<<<<<< Updated upstream
=======
from tqdm import tqdm 
>>>>>>> Stashed changes

from .groupcat import load_subvolume, file_path, load_galprop, load_snapshot_subhalos
from .utility import load_header 


def load_subvolume(base_path, subvolume, group, fields, flag=False):
    """Returns queried results for a specific subvolume

    Similar to loadObjects() from illustris_python, modified in that does not iterate over chunks and can return whole
    arrays without having to allocate space and insert into specific positions since we are loading one file.
    :param base_path: base path to data repository
    :param subvolume: what subvolume to to load
    :param group: what catalog to query
    :param fields: list of fields to retrieve
    :param flag: if fields need to be checked
    :return: dictionary of results
    """
    result = {}

    with h5py.File(file_path(base_path, subvolume, 'sfh'), 'r') as f:
        if flag:
            if not fields:
                fields = list(f[group].keys())

            for field in fields:
                if field not in f[group].keys():
                    raise Exception("Catalog does not have requested field [{}]!".format(field))

        for field in fields:
            result[field] = f[group][field][:]

    return result
    
    
<<<<<<< Updated upstream
def load_snapshot(base_path, snap_num, subvolumes, group, fields, flag=False):
=======
def load_snapshot(base_path, snap_num, subvolumes, group, fields, flag=False, verbose=True):
>>>>>>> Stashed changes
    n_init = []

    for subvolume in subvolumes: 
        n_init.append(load_header(base_path, subvolume, 'sfh')['Nsubgroups_ThisFile_Redshift'][snap_num])
  
    # initialize objects structure
    result = {}
    
    with h5py.File(file_path(base_path, subvolumes[0], 'sfh'), 'r') as f:
        if not fields:
            fields = list(f[group].keys())

        for field in fields:
            if field not in f[group].keys():
                raise Exception("Catalog does not have requested field [{}]!".format(field))

            shape = list(f[group][field].shape)
            shape[0] = np.sum(n_init)

            # allocate within return dict
            result[field] = np.zeros(shape, dtype=f[group][field].dtype)
    
    offset = 0

<<<<<<< Updated upstream
    for subvolume in subvolumes:
=======
    for subvolume in tqdm(subvolumes, disable = not verbose):
>>>>>>> Stashed changes
        filter_fields = load_subvolume(base_path, subvolume, 'Linkprop', fields=None, flag=True)
        
        subvol_result = load_subvolume(base_path, subvolume, group, fields, flag=False)

        idx = filter_fields['LinkpropSnapNum'] == snap_num # filter_condition
        
        for field in subvol_result.keys():
            if len(subvol_result[field].shape) != 1:
                result[field][offset:offset+n_init[0], :] = subvol_result[field][idx]
            else:
                result[field][offset:offset+n_init[0]] = subvol_result[field][idx]
                
        offset += n_init[0]
        del n_init[0]
        
    return result


def load_linkprop(base_path, subvolume, fields=None):        
    return load_subvolume(base_path, subvolume, 'Linkprop', fields, flag=True)

def load_histprop(base_path, subvolume, fields=None):
    return load_subvolume(base_path, subvolume, 'Histprop', fields, flag=True)

def load_snapshot_linkprop(base_path, snap_num, subvolume, fields=None):
    return load_snapshot(base_path, snap_num, subvolume, 'Linkprop', fields, flag=True)

def load_snapshot_histprop(base_path, snap_num, subvolume, fields=None):
    return load_snapshot(base_path, snap_num, subvolume, 'Histprop', fields, flag=True)

def crossmatch_galprop(base_path, subvolume, galprop_fields=None): 
    galprop_idx = load_linkprop(base_path, subvolume, fields=['LinkproptoGalprop'])['LinkproptoGalprop']
    galprop = load_galprop(base_path, subvolume, fields=galprop_fields)

    for key in galprop.keys(): 
        galprop[key] = galprop[key][galprop_idx]
    
    return galprop

<<<<<<< Updated upstream
def crossmatch_snapshot_galprop(base_path, snap_num, subvolume, galprop_fields=None):
    galprop_idx = load_snapshot_linkprop(base_path, snap_num, subvolume, 
                                fields=['LinkproptoGalprop_Snapshot'])['LinkproptoGalprop_Snapshot']
    galprop = load_snapshot_subhalos(base_path, snap_num, subvolume, fields=galprop_fields)
=======
def crossmatch_snapshot_galprop(base_path, snap_num, subvolume, galprop_fields=None, verbose=True):
    galprop_idx = load_snapshot_linkprop(base_path, snap_num, subvolume, 
                                fields=['LinkproptoGalprop_Snapshot'])['LinkproptoGalprop_Snapshot']
    galprop = load_snapshot_subhalos(base_path, snap_num, subvolume, fields=galprop_fields, verbose=verbose)
>>>>>>> Stashed changes
    
    for key in galprop.keys(): 
        galprop[key] = galprop[key][galprop_idx]
        
    return galprop