import h5py
import numpy as np
from .groupcat import load_subvolume, file_path


def find_most_massive(forest, group):
    """function for loading the most massive progenitors from the trees in the provided forest

    :param forest: dictionary of trees to query
    :param group: what group was queried from the hdf5 files
    :return:
    """

    for tree in forest:
        most_massive_idx = []

        for snap in np.sort(np.unique(forest[tree]['%sSnapNum' % group]))[::-1]:
            mvir_idx = np.flatnonzero(forest[tree]['%sSnapNum' % group] == snap)
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
    if ('%sSnapNum' % group) not in fields:
        fields.append(group + 'SnapNum')
    if ('%sRedshift' % group) not in fields:
        fields.append(group + 'Redshift')
    return fields


def load_tree(base_path, halo_id, group, fields, matches, most_massive):
    """Returns merger queried merger tree from haloprop or galprop

    :param base_path: base path to data repository
    :param halo_id: root halo id of the tree you want to load
    :param group: what group you want to query from the hdf5 files
    :param fields: list of fields to query
    :param matches:
    :param most_massive: if you want to load most massive progenitors only
    :return:
    """
    if halo_id == -1:
        raise Exception("No halo id specified!")

    # if fields is none then all fields are loaded regardless
    if fields is not None and most_massive:
        fields = field_check(fields, group)

    tree_dict = {}

    with h5py.File('{}/lookup/tree_lookup.hdf5'.format(base_path), 'r') as f:
        tree_dict[halo_id] = {}
        loc = np.where(f['Lookup_Table']['RootHaloID'][:] == halo_id)[0][0]
        tree_dict['Subvolume'] = f['Lookup_Table']['Subvolume'][loc]
        tree_dict['%sOffsets' % group] = f['Lookup_Table']['%sOffsets' % group][loc]

    result = load_subvolume(base_path, tree_dict['Subvolume'][:], group, fields, matches, True)

    if fields is None:
        fields = result.keys()

    # load the offset file to do first round of index filtering # might need to encase tree_dict.. in tuple
    with h5py.File(file_path(base_path, tree_dict['Subvolume'][:], 'tree_offsets'), 'r') as f:
        idx = f['Offsets']['%sOffsets' % group][list(tree_dict['%sOffsets' % group])]

        for field in fields:
            if len(result[field].shape) != 1:
                result[field] = result[field][idx[0]:idx[1], :]
            else:
                result[field] = result[field][idx[0]:idx[1]]

    tree_result = {halo_id: result}

    if most_massive:
        tree_result = find_most_massive(tree_result, group)

    return tree_result


def load_forest(base_path, subvolumes, group, fields, matches, most_massive):
    """Returns merger queried forest of merger trees from haloprop or galprop

    :param base_path: base path to data repository
    :param subvolumes: subvolumes to load and query
    :param group: what group you want to query from the hdf5 files
    :param fields: list of fields to query
    :param matches:
    :param most_massive: if you want to load most massive progenitors only
    :return:
    """
    if fields is not None and most_massive:
        fields = field_check(fields, group)

    forest_result = {}

    with h5py.File('{}/lookup/tree_lookup.hdf5'.format(base_path), 'r') as f:
        for subvolume in subvolumes:
            result = load_subvolume(base_path, subvolume, group, fields, matches, True)

            with h5py.File(file_path(base_path, subvolume, 'tree_offsets'), 'r') as g:
                subvolume_loc = np.where((f['Lookup_Lookup']['Subvolume'][:] == subvolume).all(axis=1))[0][0]
                subvolume_offset = f['Lookup_Lookup']['Offsets'][subvolume_loc]

                trees = f['Lookup_Table']['RootHaloID'][subvolume_offset[0]:subvolume_offset[1]]
                offsets = f['Lookup_Table']['%sOffsets' % group][subvolume_offset[0]:subvolume_offset[1]]

                for i in range(0, len(trees)):
                    forest_result[trees[i]] = {}
                    idx = g['Offsets']['%sOffsets' % group][offsets[i][0]:offsets[i][1]]
                    for field in result.keys():
                        if len(result[field].shape) != 1:
                            forest_result[trees[i]][field] = result[field][idx, :]
                        else:
                            forest_result[trees[i]][field] = result[field][idx]

    if most_massive:
        forest_result = find_most_massive(forest_result, group)

    return forest_result


def load_tree_galprop(base_path, halo_id=-1, fields=None, matches=False, most_massive=False):
    """Returns a specific subhalo merger tree."""
    return load_tree(base_path, halo_id, 'Galprop', fields, matches, most_massive)


def load_tree_haloprop(base_path, halo_id=-1, fields=None, matches=False, most_massive=False):
    """Returns a specific halo merger tree."""
    return load_tree(base_path, halo_id, 'Haloprop', fields, matches, most_massive)


def load_galprop_forest(base_path, subvolumes, fields=None, matches=False, most_massive=False):
    """Returns all subhalo trees from queried subvolumes."""
    return load_forest(base_path, subvolumes, 'Galprop', fields, matches, most_massive)


def load_haloprop_forest(base_path, subvolumes, fields=None, matches=False, most_massive=False):
    """Returns all halo trees from queried subvolumes."""
    return load_forest(base_path, subvolumes, 'Haloprop', fields, matches, most_massive)
