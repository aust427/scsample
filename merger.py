import h5py
import numpy as np
from .groupcat import load_subvolume

TREE_PATH = 'postprocessing/tree_offsets'


def find_most_massive(forest, group):
    """function for loading the most massive progenitors from the trees in the provided forest

    :param forest: dictionary of trees to query
    :param group: what group was queried from the hdf5 files
    :return:
    """
    most_massive_idx = []

    for tree in forest:
        for z in np.sort(np.unique(forest[tree]['%sRedshift' % group])):
            mvir_idx = np.flatnonzero(forest[tree]['%sRedshift' % group] == z)
            most_massive_idx.append(mvir_idx[np.argmax(forest[tree]['%sMvir' % group][mvir_idx])])

        for field in forest[tree].keys():
            if len(forest[tree][field].shape) != 1:
                forest[tree][field] = forest[tree][field][most_massive_idx, :]
            else:
                forest[tree][field] = forest[tree][field][most_massive_idx]

    return forest


def field_check(fields, group):
    """Helper function for checking if the right fields are included when querying for most massive progenitor

    :param group: what group you want to query from the hdf5 files
    :param fields: list of fields to query
    :return:
    """
    if ('%sMvir' % group) not in fields:
        fields.append(group + 'Mvir')
    if ('%sRedshift' % group) not in fields:
        fields.append(group + 'Redshift')
    return fields


def load_tree(base_path, halo_id, group, fields=None, most_massive=False):
    """Returns merger queried merger tree from haloprop or galprop

    :param base_path: base path to data repository
    :param halo_id: root halo id of the tree you want to load
    :param group: what group you want to query from the hdf5 files
    :param fields: list of fields to query
    :param most_massive: if you want to load most massive progenitors only
    :return:
    """
    if halo_id == -1:
        raise Exception("No halo id specified!")

    # if fields is none then all fields are loaded regardless
    if fields is not None and most_massive:
        fields = field_check(fields, group)

    tree_dict = {halo_id: {}}

    with h5py.File('%s%s/offsets_lookup.hdf5' % (base_path, TREE_PATH), 'r') as f:
        tree_dict[halo_id] = {}
        loc = np.where(f['Lookup_Table']['GalpropRootHaloID'][:] == halo_id)[0][0]
        for key in f['Lookup_Table']:
            tree_dict[halo_id][key] = (f['Lookup_Table'][key][loc])

    result = load_subvolume(base_path, tree_dict[halo_id]['Subvolume'][:], group, fields, True)

    if fields is None:
        fields = result.keys()

    # load the offset file to do first round of index filtering
    with h5py.File('%s%s%s' % (base_path, TREE_PATH, '/offsets_%i_%i_%i.hdf5'
                                                     % tuple(tree_dict[halo_id]['Subvolume'][:])), 'r') as f:
        offset = f['Offsets']['%sOffsets' % group][tree_dict[halo_id]['%sOffsets' % group]]
        for field in fields:
            if len(result[field].shape) != 1:
                result[field] = result[field][offset[0]:offset[1], :]
            else:
                result[field] = result[field][offset[0]:offset[1]]

    tree_result = {halo_id: result}

    if most_massive:
        tree_result = find_most_massive(tree_result, group)

    return tree_result


def load_tree_catalog(base_path, subvolumes, group, fields=None, most_massive=False):
    """Returns merger queried forest of merger trees from haloprop or galprop

    :param base_path: base path to data repository
    :param subvolumes: subvolumes to load and query
    :param group: what group you want to query from the hdf5 files
    :param fields: list of fields to query
    :param most_massive: if you want to load most massive progenitors only
    :return:
    """
    if fields is not None and most_massive:
        fields = field_check(fields, group)

    forest_result = {}

    for subvolume in subvolumes:
        result = load_subvolume(base_path, subvolume, group, fields, True)

    if most_massive:
        forest_result = find_most_massive(forest_result, group)

    return forest_result


def load_tree_galprop(base_path, halo_id=-1, fields=None, most_massive=False):
    """Returns a specific subhalo merger tree."""
    return load_tree(base_path, halo_id, 'Galprop', fields, most_massive)


def load_tree_haloprop(base_path, halo_id=-1, fields=None, most_massive=False):
    """Returns a specific halo merger tree."""
    return load_tree(base_path, halo_id, 'Haloprop', fields, most_massive)


def load_galprop_trees(base_path, subvolumes, fields=None, most_massive=False):
    """Returns all subhalo trees from queried subvolumes."""
    return load_tree_catalog(base_path, subvolumes, 'Galprop', fields, most_massive)


def load_haloprop_trees(base_path, subvolumes, fields=None, most_massive=False):
    """Returns all halo trees from queried subvolumes."""
    return load_tree_catalog(base_path, subvolumes, 'Haloprop', fields, most_massive)
