# Learning representations for protein structural fragments

## Priorities ##

2.  Add tests to architectures.py in main
	 - using pytest
5.  document all functions and classes
7.  more efficient fragment cosntruction method
8.  make method for c-beta extended
9.  add fragment method for bb+sc
10. make a cif writer
11. improve pymol visualizationn
	 - probably focus on getting pymol to work in a notebook
12. visualize nn input/outputs overlays
13. Build more informative README
	 - include file hiearchy
	 - how to use




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

