import h5py
import numpy as np
from tqdm import tqdm 

from .utility import load_header, file_path


# def load_matches(base_path, subvolume, group):
#     """

#     :param base_path:
#     :param subvolume:
#     :param group:
#     :return:
#     """
#     result = {}

#     with h5py.File(file_path(base_path, file_name='matches'), 'r') as f:
#         fields = list(f[group].keys())

#         for field in fields:
#             result[field] = f[group][field][:]

#     return result


def load_subvolume(base_path, subvolume, group, fields, matches, flag, file_name='volume'):
    """Returns queried results for a specific subvolume

    Similar to loadObjects() from illustris_python, modified in that does not iterate over chunks and can return whole
    arrays without having to allocate space and insert into specific positions since we are loading one file.
    :param base_path: base path to data repository
    :param subvolume: what subvolume to to load
    :param group: what catalog to query
    :param fields: list of fields to retrieve
    :param matches:
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

        if matches:
            result = {**result, **load_matches(base_path, subvolume, group)}

    return result


def load_snapshot(base_path, snap_num, subvolumes, group, fields, matches, verbose=True, file_name='volume'):
    """Returns objects queried for all subvolumes

    :param base_path: base path to data repository
    :param snap_num: snapshot number you want to load
    :param subvolumes: list of subvolumes to query
    :param group: what group you want to query from the hdf5 files
    :param fields: list of fields to query
    :param matches:
    :return:
    """
    n_init = []

    snap_key = 'N{}_ThisSubvol_Redshift'.format('groups' if group == 'Haloprop' else 'subgroups')
    for subvolume in subvolumes: 
        n_init.append(load_header(base_path, subvolume, file_name)[snap_key][snap_num])
        
    # initialize objects structure
    result = {}
    
    with h5py.File(file_path(base_path, file_name), 'r') as f:
        subvol = f['{}_{}_{}'.format(*subvolumes[0])]
        # using snapshot in the new update
        filter_field = '{}SnapNum'.format(group)
        
        if not fields:
            fields = list(subvol[group].keys())

        # make sure the snapshot field is included in fields
        if filter_field not in fields:
            fields.append(filter_field)   
            
        for field in fields:
            if field not in subvol[group].keys():
                raise Exception("Catalog does not have requested field [{}]!".format(field))

            shape = list(subvol[group][field].shape)
            shape[0] = np.sum(n_init)

            # allocate within return dict
            result[field] = np.zeros(shape, dtype=subvol[group][field].dtype)

    # if matches:
    #     with h5py.File(file_path(base_path, subvolumes[0], 'matches'), 'r') as f:
    #         for field in f[group].keys():
    #             result[field] = np.zeros(shape, dtype=f[group][field].dtype)

    # header = load_header(base_path, subvolumes[0])
    # filter_condition = header['Redshifts'][snap_num]

    offset = 0

    for subvolume in tqdm(subvolumes, disable = not verbose):
        subvol_result = load_subvolume(base_path, subvolume, group, fields, matches, False, file_name=file_name)

        idx = subvol_result[filter_field][:] == snap_num # filter_condition

        for field in subvol_result.keys():
            if len(subvol_result[field].shape) != 1:
                result[field][offset:offset+n_init[0], :] = subvol_result[field][idx]
            else:
                result[field][offset:offset+n_init[0]] = subvol_result[field][idx]

        offset += n_init[0]
        del n_init[0]
                
    return result


def load_haloprop(base_path, subvolume, fields=None, matches=False, verbose=True, file_name='volume'):
    """Returns a specific subvolume's haloprop for all snapshots."""
    return load_subvolume(base_path, subvolume, 'Haloprop', fields, matches, True, file_name=file_name)


def load_galprop(base_path, subvolume, fields=None, matches=False, verbose=True, file_name='volume'):
    """Returns a specific subvolume's galprop for all snapshots."""
    return load_subvolume(base_path, subvolume, 'Galprop', fields, matches, True, file_name=file_name)

    
def load_snapshot_halos(base_path, snap_num, subvolumes, fields=None, matches=False, verbose=True, file_name='volume'):
    """Returns all halos from queried subvolumes at a specific snapshot."""
    return load_snapshot(base_path, snap_num, subvolumes, "Haloprop", fields, matches, verbose, file_name=file_name)


def load_snapshot_subhalos(base_path, snap_num, subvolumes, fields=None, matches=False, verbose=True, file_name='volume'):
    """Returns all subhalos from queried subvolumes at a specific snapshot."""
    return load_snapshot(base_path, snap_num, subvolumes, "Galprop", fields, matches, verbose, file_name=file_name)



