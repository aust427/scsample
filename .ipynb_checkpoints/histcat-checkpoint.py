import h5py
import numpy as np
from tqdm import tqdm 

from .groupcat import load_subvolume, file_path, load_galprop, load_snapshot_subhalos
from .utility import load_header 


def load_subvolume(base_path, subvolume, group, fields, flag=False, file_name='volume'):
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

    with h5py.File(file_path(base_path, file_name=file_name), 'r') as f:
        subvol = f['{}_{}_{}'.format(*subvolume)]    
        if flag:
            if not fields:
                fields = list(subvol[group].keys())

            for field in fields:
                if field not in subvol[group].keys():
                    raise Exception("Catalog does not have requested field [{}]!".format(field))

        for field in fields:
            result[field] = subvol[group][field][:]

    return result
    
    
def load_snapshot(base_path, snap_num, subvolumes, group, fields, flag=False, verbose=True, file_name='volume'):
    n_init = []

    for subvolume in subvolumes: 
        n_init.append(load_header(base_path, subvolume, file_name)['Nsubgroups_ThisSubvol_Redshift_SFH'][snap_num])
  
    # initialize objects structure
    result = {}
    
    with h5py.File(file_path(base_path, file_name), 'r') as f:
        subvol = f['{}_{}_{}'.format(*subvolumes[0])]
        if not fields:
            fields = list(subvol[group].keys())

        for field in fields:
            if field not in subvol[group].keys():
                raise Exception("Catalog does not have requested field [{}]!".format(field))

            shape = list(subvol[group][field].shape)
            shape[0] = np.sum(n_init)

            # allocate within return dict
            result[field] = np.zeros(shape, dtype=subvol[group][field].dtype)
    
    offset = 0

    for subvolume in tqdm(subvolumes, disable = not verbose):
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


def load_linkprop(base_path, subvolume, fields=None, file_name='volume'):        
    return load_subvolume(base_path, subvolume, 'Linkprop', fields, flag=True, file_name=file_name)

def load_histprop(base_path, subvolume, fields=None, file_name='volume'):
    return load_subvolume(base_path, subvolume, 'Histprop', fields, flag=True, file_name=file_name)

def load_snapshot_linkprop(base_path, snap_num, subvolume, fields=None, file_name='volume'):
    return load_snapshot(base_path, snap_num, subvolume, 'Linkprop', fields, flag=True, file_name=file_name)

def load_snapshot_histprop(base_path, snap_num, subvolume, fields=None, file_name='volume'):
    return load_snapshot(base_path, snap_num, subvolume, 'Histprop', fields, flag=True, file_name=file_name)

def crossmatch_galprop(base_path, subvolume, galprop_fields=None, file_name='volume'): 
    galprop_idx = load_linkprop(base_path, subvolume, 
                                fields=['LinkproptoGalprop'], file_name=file_name)['LinkproptoGalprop']
    galprop = load_galprop(base_path, subvolume, fields=galprop_fields)

    for key in galprop.keys(): 
        galprop[key] = galprop[key][galprop_idx]
    
    return galprop

def crossmatch_snapshot_galprop(base_path, snap_num, subvolume, galprop_fields=None, verbose=True, file_name='volume'):
    galprop_idx = load_snapshot_linkprop(base_path, snap_num, subvolume, 
                                fields=['LinkproptoGalprop_Snapshot'], file_name=file_name)['LinkproptoGalprop_Snapshot']
    galprop = load_snapshot_subhalos(base_path, snap_num, subvolume, fields=galprop_fields, verbose=verbose)
    
    for key in galprop.keys(): 
        galprop[key] = galprop[key][galprop_idx]
        
    return galprop