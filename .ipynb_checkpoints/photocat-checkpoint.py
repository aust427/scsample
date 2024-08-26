import h5py
import numpy as np

from .groupcat import load_subvolume, file_path, load_header


def load_photometry_subvolume(base_path, subvolume, group, fields, flag=False):
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

    with h5py.File(file_path(base_path, subvolume, 'photometry'), 'r') as f:
        if flag:
            if not fields:
                fields = list(f[group].keys())

            for field in fields:
                if field not in f[group].keys():
                    raise Exception("Catalog does not have requested field [{}]!".format(field))

        for field in fields:
            result[field] = f[group][field][:]

    return result
    
    
def load_photometry_snapshot(base_path, snap_num, subvolumes, group, fields, flag=False):
    n_init = []

    for subvolume in subvolumes: 
        n_init.append(load_header(base_path, subvolume, 'photometry')['Nsubgroups_ThisFile_Redshift'][snap_num])
  
    # initialize objects structure
    result = {}
    
    with h5py.File(file_path(base_path, subvolumes[0], 'photometry'), 'r') as f:
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

    for subvolume in subvolumes:
        filter_fields = load_photometry_subvolume(base_path, subvolume, 'Linkprop', fields=None, flag=True)
        
        subvol_result = load_photometry_subvolume(base_path, subvolume, group, fields, flag=False)

        idx = filter_fields['LinkpropSnapNum'] == snap_num # filter_condition
        
        # result['LinkpropGalpropIndex'] = filter_fields['LinkpropGalpropIndex'][idx]
        
        # for field in subvol_result.keys():
        #     result[field] = subvol_result[field][idx]
        
        for field in subvol_result.keys():
            if len(subvol_result[field].shape) != 1:
                result[field][offset:offset+n_init[0], :] = subvol_result[field][idx]
            else:
                result[field][offset:offset+n_init[0]] = subvol_result[field][idx]
                
        offset += n_init[0]
        del n_init[0]
        
    return result


def load_linkprop(base_path, subvolume, fields=None):
    return load_photometry_subvolume(base_path, subvolume, 'Linkprop', fields, flag=True)

def load_histprop(base_path, subvolume, fields=None):
    return load_photometry_subvolume(base_path, subvolume, 'Histprop', fields, flag=True)

def load_photoprop(base_path, subvolume, fields=None):
    return load_photometry_subvolume(base_path, subvolume, 'Photoprop', fields, flag=True)

def load_snapshot_linkprop(base_path, snap_num, subvolume, fields=None):
    return load_photometry_snapshot(base_path, snap_num, subvolume, 'Linkprop', fields, flag=True)

def load_snapshot_histprop(base_path, snap_num, subvolume, fields=None):
    return load_photometry_snapshot(base_path, snap_num, subvolume, 'Histprop', fields, flag=True)

def load_snapshot_photoprop(base_path, snap_num, subvolume, fields=None):
    return load_photometry_snapshot(base_path, snap_num, subvolume, 'Photoprop', fields, flag=True)