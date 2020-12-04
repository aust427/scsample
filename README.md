# illustris_sam
Python module developed for loading and querying the Illustris-SAM catalogs. Heavily inspired by the illustris-python module created by Dylan Nelson. Created by Austen Gabrielpillai. 

## Background
The module is compatible with the following soon to be published datasets: 

```
Simulation   | n_subvolumes | Sizes (GB) | Snapshots | nhalos | nsubhalos 
TNG100-1-SAM | 125          |            | 0 - 99    |        |
TNG300-1-SAM | 343          |            | 1 - 99    |        |
TNG50-1-SAM  | 216          |            | 0 - 99    |        |
```

If using the above datasets, please cite ... For more information ... 

## Installing 


## Directory File Structure 
A brief explanation of the files created for this project that will be used by this module:
 
	subvolume.hdf5 - Post-processed output file from running the Santa-Cruz SAM on Consistent-Tree files created from running Rockstar on IllustrisTNG. Each subvolume is a partition of the IllustrisTNG volume. Contains two catalogs: Galprop (galaxy properties) and Haloprop (halo properties). For a list of fields, please see [link].
	tree-offsets.hdf5 - Pre-computed offsets used for querying merger trees. 
	matches.hdf5 - Bijective matches between our SAM subhalos / halos and TNG's particle subhalos / halos created by running the SubLink algorithm on Rockstar and Subfind catalogs. For a list of fields, please see [link].  
	tree_lookup.hdf5 - Lookup file for querying individual and sets of trees across subvolumes. 
	
To use this module, we follow a similar structure outline in the TNG Project's script examples, and it is assumed that you have downloaded the proper files and have set up your directories accordingly: 

``` 
TNG-SAM/
TNG-SAM/L75n1820TNG/
TNG-SAM/L75n1820TNG/output/
	subvolume catalog: TNG-SAM/L75n1820TNG/output/0_0_0/subvolume.hdf5
	merger tree offsets: TNG-SAM/L75n1820TNG/output/0_0_0/tree-offsets.hdf5
	TNG bijective matches: TNG-SAM/L75n1820TNG/output/0_0_0/matches.hdf5 (OPTIONAL)
	tree lookup file: TNG-SAM/L75n1820TNG/output/lookup/tree_lookup.hdf5
```
	 
From here, you would set basePathSAM = 'TNG-SAM/L75n1820TNG'. For more detailed instruction from here, please take a look at the provided Jupyter Notebook as an example of the capabilities of both the module and the catalog. 