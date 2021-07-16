# protein autoencoder

## Priorities ##

1.  architectures.py set up
	 - contain all PyTorch Model Class definitions we are using
	 - UPDATE: architectures.py set up, but not tested
2.  Add tests to architectures.py in main
	 - use unittests or pytest?
3.  add CNN definitions to architectures.py
4.  add method to create distance matrices from fragments
	 - test simple CNN with distance matrices
5.  document all functions and classes
6.  implement PEP8 style checking in repo
7.  more efficient fragment cosntruction method
8.  make method for c-beta extended
9.  add fragment method for bb+sc
10.  make a cif writer
11.  improve pymol visualizationn
	 - probably focus on getting pymol to work in a notebook
12. visualize nn input/outputs overlays
13. Build more informative README
	 - include file hiearchy
	 - how to use

## Plan Overview ##

+ Fragment sizes = 3, 5, 7, 9, 11, 15, 21, 25, 30
+ fragment types = ca, bb, bb+cen, bb+sc
+ selective fragment sets
	+ all alpha
	+ all beta
	+ all turns
	+ all coils
	+ helix-turn-helix
+ FC autoencoders 
+ CNN distance matrix autoencoders
+ find places of low/high quality

