import numpy as np
import h5py 



def gen_subvolume_list(nsub):
    n = int(np.round(nsub**(1/3))) 
    
    return [[i, j, k] for i in range(n) for j in range(n) for k in range(n)]


def load_header(base_path, subvolume, file_name='volume'):
    """Returns the header from az queried subvolume."""
    with h5py.File(file_path(base_path, file_name=file_name ), 'r') as f:
        subvol = f['{}_{}_{}'.format(*subvolume)]
        header = dict(subvol['Header'].attrs.items())
        header.update({key: subvol['Header'][key][:] for key in subvol['Header'].keys()})
        
    return header


def file_path(base_path, file_name):
    """Returns the path to the subvolume hdf5 file

    Similar to gcPath() from illustris_python, modified to load specific subvolumes.
    :param base_path: base path to data repository
    :param subvolume: what subvolume to load
    :param file_name:
    :return: path to file
    """
    return '{}/{}.hdf5'.format(base_path, file_name)



def get_metadata(base_path, keys): 
    g_keys = []
    ha_keys = []
    l_keys = []
    hi_keys = []
    
    for key in keys: 
        if key[0:3] == 'Gal': 
            g_keys.append(key)
        elif key[0:4] == 'Halo': 
            ha_keys.append(key)
        elif key[0:4] == 'Link': 
            l_keys.append(key)
        elif key[0:4] == 'Hist': 
            hi_keys.append(key)
        else: 
            print('Key not in any file!') 
    
    d = []
    
    with h5py.File(file_path(base_path, 'volume'), 'r') as f: 
        subvol = f['{}'.format(list(f.keys())[0])] # grab the first subvolume 

        for key in g_keys: 
            d.append([key] + [x for xs in list(subvol['Galprop'][key].attrs.items()) for x in list(xs)])
        
        for key in ha_keys: 
            d.append([key] + [x for xs in list(subvol['Haloprop'][key].attrs.items()) for x in list(xs)])
    
        for key in l_keys: 
            d.append([key] + [x for xs in list(subvol['Linkprop'][key].attrs.items()) for x in list(xs)])
            
        for key in hi_keys: 
            d.append([key] + [x for xs in list(subvol['Histprop'][key].attrs.items()) for x in list(xs)])
            
    return np.array(d)[:, [0, 4, 2]]


def print_units(base_path, keys): 
    d = get_metadata(base_path, keys)[:, [0, 1]]

    max_key = max([len(d[i, 0]) for i in range(d.shape[0])])
    
    print(" | ".join([str(l).ljust(s) for l, s in zip(['Key', 'Units'], [max_key, -1])]))
    _ = [print(" | ".join([str(l).ljust(s) for l, s in zip(i, [max_key, -1])])) for i in d]
        
    return

def print_desc(base_path, keys): 
    d = get_metadata(base_path, keys)[:, [0, 2]]
    
    max_key = max([len(d[i, 0]) for i in range(d.shape[0])])

    print(" | ".join([str(l).ljust(s) for l, s in zip(['Key', 'Desc.'], [max_key, -1])]))
    _ = [print(" | ".join([str(l).ljust(s) for l, s in zip(i, [max_key, -1])])) for i in d]

    return 


def print_metadata(base_path, keys):
    d = get_metadata(base_path, keys)
    
    max_key = max([len(d[i, 0]) for i in range(d.shape[0])])
    max_unit = max([len(d[i, 1]) for i in range(d.shape[0])])

    _ = [print(" | ".join([str(l).ljust(s) for l, s in zip(i, [max_key, max_unit, -1])])) for i in d]
        
    return 


def result_to_df(res):
    return 

