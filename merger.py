import h5py
import numpy as np
from .groupcat import load_subvolume

TREE_PATH = '/postprocessing/tree_offsets/'


def find_most_massive(forest):
    for tree in forest:
        print(tree)
    return 0


def field_check(fields, group):
    if (group + 'Mvir') not in fields:
        fields.append(group + 'Mvir')
    if (group + 'Redshift') not in fields:
        fields.append(group + 'Redshift')
    return fields


def load_tree(base_path, halo_id, group, fields=None, most_massive=False):
    if halo_id == -1:
        raise Exception("No halo id specified!")

    if fields is not None and most_massive:
        fields = field_check(fields, group)

    tree_dict = {halo_id: {}}

    with h5py.File(base_path + TREE_PATH + 'offsets_lookup.hdf5', 'r') as f:
        tree_dict[halo_id] = {}
        loc = np.where(f['Lookup_Table']['GalpropRootHaloID'][:] == halo_id)[0][0]
        for key in f['Lookup_Table']:
            tree_dict[halo_id][key] = (f['Lookup_Table'][key][loc])

    result = load_subvolume(base_path, tree_dict[halo_id]['Subvolume'][:], group, fields, True)

    if fields is None:
        fields = result.keys()

    with h5py.File(base_path + TREE_PATH + 'offsets_%i_%i_%i.hdf5'
                   % tuple(tree_dict[halo_id]['Subvolume'][:]), 'r') as f:
        offset = f['Offsets'][group + 'Offsets'][tree_dict[halo_id][group + 'Offsets']]
        for field in fields:
            if len(result[field].shape) != 1:
                result[field] = result[field][offset[0]:offset[1], :]
            else:
                result[field] = result[field][offset[0]:offset[1]]

    tree_result = {halo_id: result}

    if most_massive:
        tree_result = find_most_massive(tree_result)

    return tree_result


def load_tree_catalog(base_path, subvolume, halo_id, group, fields=None, most_massive=False):
    return 0


def load_tree_galprop(base_path, halo_id=-1, fields=None, most_massive=False):
    """Returns a specific subhalo merger tree."""
    return load_tree(base_path, halo_id, 'Galprop', fields, most_massive)


def load_tree_haloprop(base_path, halo_id=-1, fields=None, most_massive=False):
    """Returns a specific halo merger tree."""
    return load_tree(base_path, halo_id, 'Haloprop', fields, most_massive)


def load_galprop_trees(base_path, subvolume, fields=None, most_massive=False):
    """Returns all subhalo trees from queried subvolumes."""
    return load_tree_catalog()


def load_haloprop_trees(base_path, subvolume, fields=None, most_massive=False):
    """Returns all halo trees from queried subvolumes."""
    return load_tree_catalog()
