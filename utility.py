import numpy as np
import h5py 



def gen_subvolume_list(nsub):
    n = int(np.round(nsub**(1/3))) 
    
    return [[i, j, k] for i in range(n) for j in range(n) for k in range(n)]


def load_header(base_path, subvolume, file_name):
    """Returns the header from a queried subvolume."""
    with h5py.File(file_path(base_path, subvolume, file_name), 'r') as f:
        header = dict(f['Header'].attrs.items())
        header.update({key: f['Header'][key][:] for key in f['Header'].keys()})
        
    return header


def file_path(base_path, subvolume, file_name):
    """Returns the path to the subvolume hdf5 file

    Similar to gcPath() from illustris_python, modified to load specific subvolumes.
    :param base_path: base path to data repository
    :param subvolume: what subvolume to load
    :param file_name:
    :return: path to file
    """
    return '{}/{}_{}_{}/{}.hdf5'.format(base_path, *subvolume, file_name)


def result_to_df(res):
    return 