import h5py
import numpy as np


def subvolume_path(base_path, subvolume):
    """Returns the path to the subvolume hdf5 file

    Similar to gcPath() from illustris_python, modified to load specific subvolumes.
    :param base_path: base path to data repository
    :param subvolume: what subvolume you want to load
    :return: path to subvolume file
    """
    catalog_path = base_path + '/%i_%i_%i/' % (subvolume[0], subvolume[1], subvolume[2])
    file_path = catalog_path + 'subvolume.hdf5'

    return file_path


def load_subvolume_objects(base_path, subvolume, group, fields, flag):
    """Returns queried results for a specific subvolume

    Similar to loadObjects() from illustris_python, modified in that does not iterate over chunks and can return whole
    arrays without having to allocate space and insert into specific positions since we are loading one file.
    :param base_path: base path to data repository
    :param subvolume: what subvolume you want to load
    :param group: what catalog to query
    :param fields: list of fields to retrieve
    :param flag: if fields need to be checked
    :return: dictionary of results
    """
    f = h5py.File(subvolume_path(base_path, subvolume), 'r')

    if flag:
        header = dict(f['Header'].attrs.items())

        if not fields:
            fields = list(f[group].keys())

        for field in fields:
            if field not in f[group].keys():
                raise Exception("Catalog does not have requested field [" + field + "]!")

    result = {}

    for field in fields:
        result[field] = f[group][field][:]

    return result


def load_haloprop(base_path, subvolume, fields=None):
    """Returns a specific subvolume's haloprop for all snapshots.

    Similar to loadHalos() from illustris_python except for a specific subvolume rather than snapshot.
    :param base_path: base path to data repository
    :param subvolume: what subvolume you want to load
    :param fields: fields that you want to query. If none, returns all possible fields
    :return: dictionary of results
    """
    return load_subvolume_objects(base_path, subvolume, 'Group', fields, True)


def load_galprop(base_path, subvolume, fields=None):
    """Returns a specific subvolume's galprop for all snapshots.

    Similar to loadSubhalos() from illustris_python except for a specific subvolume rather than snapshot.
    :param base_path: base path to data repository
    :param subvolume: what subvolume you want to load
    :param fields: fields that you want to query. If none, returns all possible fields
    :return: dictionary of results
    """
    return load_subvolume_objects(base_path, subvolume, 'Subhalo', fields, True)


def load_snap_objects(base_path, snap_num, subvolumes, group, fields=None):
    """Returns objects queried for all subvolumes

    :param base_path: base path to data repository
    :param snap_num: snapshot number you want to load
    :param subvolumes:
    :param group: what group you want to query from the hdf5 files
    :param fields:
    :return:
    """
    with h5py.File(subvolume_path(base_path, subvolumes[0])) as f:

        header = dict(f['Header'].attrs.items())

        n_partitions = np.cbrt(header['Nsubvolumes'])

        if not fields:
            fields = list(f[group].keys())

        for field in fields:
            if field not in f[group].keys():
                raise Exception("Catalog does not have requested field [" + field + "]!")

            # pop the exception?

    if group == 'galprop':
        filter_field = 'SubhaloRedshift'
        filter_condition = header['redshifts'][snap_num]
    else:
        filter_field = 'GroupSnapNum'
        filter_condition = snap_num

    if filter_field not in fields:
        fields.append(filter_field)

    # initialize
    objects = {}

    for field in fields:
        object[field] = np.array([])

    for i in subvolumes:
        subvolume_objects = load_subvolume_objects(base_path, subvolumes[i], group, fields, True)

        idx = subvolume_objects[filter_field] == filter_condition

        # speed test for np.append later ...
        for field in fields:
            objects[field] = np.append(object[field], subvolume_objects[field][idx])

    return


def load_snap_subhalos(base_path, snap_num, subvolumes, fields=None):
    """Returns all subhalos from all subvolumes at a given snapshot.


    :param base_path:
    :param snap_num:
    :param subvolumes:
    :param fields:
    :return:
    """
    return load_snap_objects(base_path, snap_num, subvolumes, "galprop", fields)


def load_snap_halos(base_path, snap_num, subvolumes, fields=None):
    """Returns all halos from all subvolumes at a given snapshot.

    :param base_path:
    :param snap_num:
    :param subvolumes:
    :param fields:
    :return:
    """
    return load_snap_objects(base_path, snap_num, subvolumes, "galprop", fields)


def load_header(base_path, subvolume):
    return


def load_single(base_path, snap_num, subvolume, halo_id=-1, subhalo_id=-1):
    return
