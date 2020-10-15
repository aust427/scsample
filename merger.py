import h5py
import numpy as np
import groupcat


def load_roots(base_path):
    """Returns tree root ids and what subvolume they are in."""
    return 0


def find_most_massive():
    return 0


def load_tree(base_path, halo_id, group, fields=None, most_massive=False):
    if halo_id == -1:
        raise Exception("No halo id specified!")

    root_locations = load_roots(base_path)

    loc = '000'

    subvolume = groupcat.load_subvolume(base_path, [int(char) for char in loc], group, fields, True)

    return 0


def load_tree_galprop(base_path, halo_id=-1, fields=None, most_massive=False):
    """Returns a specific subhalo merger tree."""
    return load_tree(base_path, halo_id, 'Galprop', fields, most_massive)


def load_tree_haloprop(base_path, halo_id=-1, fields=None, most_massive=False):
    """Returns a specific halo merger tree."""
    return load_tree(base_path, halo_id, 'Haloprop', fields, most_massive)


def load_galprop_trees(base_path, subvolume, fields=None, most_massive=False):
    """Returns all subhalo trees from queried subvolumes."""
    return 0


def load_haloprop_trees(base_path, subvolume, fields=None, most_massive=False):
    """Returns all halo trees from queried subvolumes."""
    return 0
