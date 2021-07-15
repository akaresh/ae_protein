#!/usr/bin/python3

#http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec190 11.1.7

import sys
import argparse
import pickle

import pandas as pd
from Bio.PDB.PDBIO import PDBIO

# print(pd.__version__)
# sys.exit()

parser = argparse.ArgumentParser(description='making fragments')
parser.add_argument('--atoms', '-a', required=False, type=str, default=None,
	metavar='<path>', help='path to prebuilt atom data frame')
parser.add_argument('--cifs', '-c', required=False, type=str,
	metavar='<path>', help='path to file containing set of cifs to build')

arg = parser.parse_args()

assert(arg.atoms != None)



atm_frame = pd.read_pickle(arg.atoms, compression='xz')

#writing in consideration only for CA dataframe
for idx, row in atm_frame.iterrows():
# 	print(row['xyz_set'], len(row['xyz_set']))
# 	print(row['fragment_seq'], len(row['fragment_seq']))
# 	sys.exit()
	for r, c in zip(row['fragment_seq'], row['xyz_set']):
		print(r,c)
	sys.exit()

# Example: saving a structure


# >>> io = PDBIO()
# >>> io.set_structure(s)
# >>> io.save("out.pdb")

#https://fitzkee.chemistry.msstate.edu/sites/default/files/ch8990/pymol-tutorial.pdf
