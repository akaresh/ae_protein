# Learning representations for protein structural fragments

## Organization ##

```
|-- data_collection
|-- fragment_generators
|-- sample_data
|-- tests
|   |-- arch_tests
|   |-- frag_gen_tests
|   |-- vis_tests
|-- trainers
|-- visualization
```

* **data\_collection**
	Contains scripts to download cif files for protein structures in dataset. 
	To get started with constructing the fragment dataset, you need a 		
	pdbid\_list.txt file from PISCES server. 
	- pull\_pdbs.py: Downloads pdbs from PISCES output. Saves the chains 
	specified in PISCES list.
* **fragment\_generators**
	Methods to produce the structural fragment dataframes. 
* **sample\_data**
	Datasets for testing. Advised to only place single `.cif` files here, or a 
	samll dataframe to benchmark architectures. 
* **tests**
	Tests for the repository. Writing with `pytest`.
* **trainers**
	Location of training scripts, architecture definitions, and tools for 
	training. 
* **visualization**
	Library for visualization of protein structure fragments. 

## Priorities ##

2.  Add tests to architectures.py in main
	 - using pytest
5.  document all functions and classes
7.  more efficient fragment cosntruction method
8.  make method for c-beta extended
9.  add fragment method for bb+sc
13. Build more informative README
	 - include file hiearchy
	 - how to use

## Development Guide ##

The repository has a Makefile to run style-checks, unit/functional tests, and manage the conda environments. 

* `make install` if you are are cloning the repo for the first time. 
	- This creates the `fqa` environment from the checked in `fqa.yaml` 		
	configuration file. 
* `make setup` if you want to update your local `fqa` environment.
* `make con` if you have installed new packages in the `fqa` environment, update 		
	the `fqa.yaml` and commit the changes. 
* `make lint` to run `pycodestyle` with our specific style. 
* `make test` to run our tests. Not implemented. 

## Plan Overview ##

+ Fragment sizes = 3, 4, 5, 7, 8, 9, 11, 15, 21, 25, 30 (11)
+ fragment types = ca, bb, bb+cen, bb+sc (4)
+ selective fragment sets (5) [44*5]
	+ all alpha
	+ all beta
	+ all turns
	+ all coils
	+ helix-turn-helix
+ FC autoencoders 
+ CNN distance matrix autoencoders
+ find places of low/high quality
	+ use lDDt?

