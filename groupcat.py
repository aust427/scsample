import h5py
import numpy as np


def subvolume_path(base_path, subvolume):
    """Returns the path to the subvolume hdf5 file

    Similar to gcPath() from illustris_python, modified to load specific subvolumes.
    :param base_path: base path to data repository
    :param subvolume: what subvolume to load
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


def load_snap_objects(base_path, snap_num, subvolumes, group, fields=None):
    """Returns objects queried for all subvolumes

    :param base_path: base path to data repository
    :param snap_num: snapshot number you want to load
    :param subvolumes: list of subvolumes to query
    :param group: what group you want to query from the hdf5 files
    :param fields: list of fields to query
    :return:
    """
    with h5py.File(subvolume_path(base_path, subvolumes[0]), 'r') as f:

        header = dict(f['Header'].attrs.items())
        header['redshifts'] = f['Header']['Redshifts'][:]

        if not fields:
            fields = list(f[group].keys())

        for field in fields:
            if field not in f[group].keys():
                raise Exception("Catalog does not have requested field [" + field + "]!")

            # pop the exception?

    # galprop and haloprop both have a redshift quantity so we can use that to query for the snapshot we want
    filter_condition = header['redshifts'][snap_num]
    filter_field = group + 'Redshift'

    # make sure the redshift field is included in fields
    if filter_field not in fields:
        fields.append(filter_field)

    # initialize objects structure
    object = {}

    # either have to post-process SAM outputs to include a table of # of halos and # of subhalos per snapshot
    # or can just begin with an empty list and append until ends
    for field in fields:
        object[field] = np.array([])

    # rather than loop over all n^3 subvolumes, user submits list of subvolumes (i.e. [[0, 0, 0], [0, 0, 1]])
    # that way it allows for usage of partial volumes if users want to download a specific subset of subvolumes
    for subvolume in subvolumes:
        subvolume_object = load_subvolume_objects(base_path, subvolume, group, fields, True)

        idx = subvolume_object[filter_field][:] == filter_condition
        for field in fields:
            object[field] = np.append(object[field], subvolume_object[field][idx])

    return object


def load_haloprop(base_path, subvolume, fields=None):
    """Returns a specific subvolume's haloprop for all snapshots."""
    return load_subvolume_objects(base_path, subvolume, 'Haloprop', fields, True)


def load_galprop(base_path, subvolume, fields=None):
    """Returns a specific subvolume's galprop for all snapshots."""
    return load_subvolume_objects(base_path, subvolume, 'Galprop', fields, True)


def load_snap_halos(base_path, snap_num, subvolumes, fields=None):
    """Returns all halos from queried subvolumes at a specific snapshot."""
    return load_snap_objects(base_path, snap_num, subvolumes, "Haloprop", fields)


def load_snap_subhalos(base_path, snap_num, subvolumes, fields=None):
    """Returns all subhalos from queried subvolumes at a specific snapshot."""
    return load_snap_objects(base_path, snap_num, subvolumes, "Galprop", fields)


def load_header(base_path, subvolume):
    return


def load_single(base_path, snap_num, subvolume, halo_id=-1, subhalo_id=-1):
    return

