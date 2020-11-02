import h5py
import numpy as np


def subvolume_path(base_path, subvolume):
    """Returns the path to the subvolume hdf5 file

    Similar to gcPath() from illustris_python, modified to load specific subvolumes.
    :param base_path: base path to data repository
    :param subvolume: what subvolume to load
    :return: path to subvolume file
    """
    file_path = base_path + 'outputs/subvolume_%i_%i_%i.hdf5' % (subvolume[0], subvolume[1], subvolume[2])

    return file_path


def load_subvolume(base_path, subvolume, group, fields, flag):
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
    f = h5py.File(subvolume_path(base_path, subvolume), 'r')

    if flag:
        if not fields:
            fields = list(f[group].keys())

        for field in fields:
            if field not in f[group].keys():
                raise Exception("Catalog does not have requested field [" + field + "]!")

    result = {}

    for field in fields:
        result[field] = f[group][field][:]

    return result


def load_snapshot(base_path, snap_num, subvolumes, group, fields=None):
    """Returns objects queried for all subvolumes

    :param base_path: base path to data repository
    :param snap_num: snapshot number you want to load
    :param subvolumes: list of subvolumes to query
    :param group: what group you want to query from the hdf5 files
    :param fields: list of fields to query
    :return:
    """
    # subvolumes = [list(subvolume) for subvolume in set(tuple(subvolume) for subvolume in subvolumes)]

    # easy to check on 1st item, might have to check all of them
    if len(subvolumes) < 1:
        raise Exception("Subvolumes is empty!")
    if type(subvolumes[0]) is not list:
        raise Exception("Subvolume is not of correct format! (not a list)")
    if len(subvolumes[0]) != 3:
        raise Exception("Subvolume is not of correct length! (length != 3)")

    n_init = []

    snap_key = 'N%s_ThisFile_Redshift' % ('groups' if group == 'Haloprop' else 'subgroups')
    for subvolume in subvolumes: 
        n_init.append(load_header(base_path, subvolume)[snap_key][snap_num])
        
    # initialize objects structure
    result = {}
    
    with h5py.File(subvolume_path(base_path, subvolumes[0]), 'r') as f:
        # galprop and haloprop both have a redshift quantity so we can use that to query for the snapshot we want
        filter_field = group + 'Redshift'
        
        if not fields:
            fields = list(f[group].keys())

        # make sure the redshift field is included in fields
        if filter_field not in fields:
            fields.append(filter_field)   
            
        for field in fields:
            if field not in f[group].keys():
                raise Exception("Catalog does not have requested field [" + field + "]!")

            shape = list(f[group][field].shape)
            shape[0] = np.sum(n_init)

            # allocate within return dict
            result[field] = np.zeros(shape, dtype=f[group][field].dtype)

    header = load_header(base_path, subvolumes[0])
    filter_condition = header['Redshifts'][snap_num]

    offset = 0

    for subvolume in subvolumes:
        subvolume_result = load_subvolume(base_path, subvolume, group, fields, False)

        idx = subvolume_result[filter_field][:] == filter_condition
        
        for field in fields:
            if len(subvolume_result[field].shape) != 1:
                result[field][offset:offset+n_init[0], :] = subvolume_result[field][idx]
            else:
                result[field][offset:offset+n_init[0]] = subvolume_result[field][idx]

        offset += n_init[0]
        del n_init[0]
        
    return result


def load_haloprop(base_path, subvolume, fields=None):
    """Returns a specific subvolume's haloprop for all snapshots."""
    return load_subvolume(base_path, subvolume, 'Haloprop', fields, True)


def load_galprop(base_path, subvolume, fields=None):
    """Returns a specific subvolume's galprop for all snapshots."""
    return load_subvolume(base_path, subvolume, 'Galprop', fields, True)


def load_snapshot_halos(base_path, snap_num, subvolumes, fields=None):
    """Returns all halos from queried subvolumes at a specific snapshot."""
    return load_snapshot(base_path, snap_num, subvolumes, "Haloprop", fields)


def load_snapshot_subhalos(base_path, snap_num, subvolumes, fields=None):
    """Returns all subhalos from queried subvolumes at a specific snapshot."""
    return load_snapshot(base_path, snap_num, subvolumes, "Galprop", fields)


def load_header(base_path, subvolume):
    """Returns the header from a queried subvolume."""
    with h5py.File(subvolume_path(base_path, subvolume), 'r') as f:
        header = dict(f['Header'].attrs.items())
        header.update({key: f['Header'][key][:] for key in f['Header'].keys()})
        
    return header


def load_single(base_path, snap_num, subvolume, halo_id=-1, subhalo_id=-1):
    return

