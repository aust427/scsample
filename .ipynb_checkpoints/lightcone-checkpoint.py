
# function to load galprop given lightcone index 
# so load subvolume, index on Q0 for ex., put to array 
# eventually sort by z 

# same-ish thing for haloprop 

# get histories for lightcone 
# need to chop sfh on redshifts  

def load_galprop(): 
    return 

def load_haloprop(): 
    return 

def load_linkprop(): 
    return 

def load_histprop(): 
    return 



# how to construct the lightcone file 
# it's probably each subvolume needs its own file, so lightcone.hdf5
# because you are going to pull out the information they need and reconstruct at the end
# code will load subvolume 1, index by Q1, load subvolume 2, index by Q1, etc. 
# organize result by redshift and return 

# similarly for hist, we just need to do some fancy masking but you can figure that out 


# this does not have to be fancy names cause it's just lookup info
# lightcone['Q1'] 
# lightcone['Q1']['pos']
# lightcone['Q1']['z_lightcone'] 
# lightcone['Q1']['subvolume_index'] 

# lightprop? 