#!/usr/bin/python3

#http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec190 11.1.7

import sys
import argparse
import pickle

import pandas as pd
import Bio.PDB.PDBIO as io

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

# with open(arg.atoms, 'wb') as f:
#     pickle.dump(data_frame_object, f)

# print(data_frame_object)
# sys.exit()

# for idx, row in atm_frame.iterrows():
# 	print(idx, row)

# Example: saving a structure


# >>> io = PDBIO()
# >>> io.set_structure(s)
# >>> io.save("out.pdb")

#https://fitzkee.chemistry.msstate.edu/sites/default/files/ch8990/pymol-tutorial.pdf
