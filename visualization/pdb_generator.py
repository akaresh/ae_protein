#!/usr/bin/python3

#http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec190 11.1.7
#http://rasbt.github.io/biopandas/tutorials/Working_with_PDB_Structures_in_DataFrames/

import sys
import argparse
import pickle

import pandas as pd
#from Bio.PDB.PDBIO import PDBIO
import Bio.Data.IUPACData as conv

# print(pd.__version__)
# sys.exit()

parser = argparse.ArgumentParser(description='making fragments')
parser.add_argument('--atoms', '-a', required=False, type=str, default=None,
	metavar='<path>', help='path to prebuilt atom data frame')
parser.add_argument('--cifs', '-c', required=False, type=str,
	metavar='<path>', help='path to file containing set of cifs to build')

arg = parser.parse_args()

assert(arg.atoms != None)

#converting to 3 letter code for pdb file
aa_dict = conv.protein_letters_1to3_extended


atm_frame = pd.read_pickle(arg.atoms, compression='xz')

#writing in consideration only for CA dataframe
for idx, row in atm_frame.iterrows():
	pdb_to_list = []
	
	for f_id, r, c in zip(row['fragment_ids'], row['fragment_seq'], row['xyz_set']):
							#atom, fragment_id, atom_type, residue_3_letter_code, chain name, 
							#residue number, coordinates
		pdb_to_list.append(['ATOM', f_id, row['fragment_type'], aa_dict[r].upper(), 
							row['chain_id'], f_id, c[0], c[1], c[2]])
		
		print(pdb_to_list)
	
	filepdb = open(f'f{idx}.pdb', 'w')
	for l in pdb_to_list:
		filepdb.write(str(l)+'\n')
		print(l)
	filepdb.close()
		
	sys.exit()

# Example: saving a structure


# >>> io = PDBIO()
# >>> io.set_structure(s)
# >>> io.save("out.pdb")

#https://fitzkee.chemistry.msstate.edu/sites/default/files/ch8990/pymol-tutorial.pdf
