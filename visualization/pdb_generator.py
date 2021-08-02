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

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='making fragments')
	parser.add_argument('--atoms', '-a', required=False, type=str, default=None,
		metavar='<path>', help='path to prebuilt atom data frame')
	parser.add_argument('--cifs', '-c', required=False, type=str,
		metavar='<path>', help='path to file containing set of cifs to build')

	arg = parser.parse_args()

	assert(arg.atoms != None)

	#converting to 3 letter code for pdb file
	aa_dict = conv.protein_letters_1to3_extended


	frag_frame = pd.read_pickle(arg.atoms, compression='xz')
	
	###pdb convertor
	#writing in consideration only for CA dataframe
	for idx, row in frag_frame.iterrows():
		# for i in row['xyz_set']:
# 			print(i)
# 		sys.exit()
		pdb_to_list = []
	
		for f_id, r, c in zip(row['fragment_ids'], row['fragment_seq'], row['xyz_set']):
								#atom, fragment_id, atom_type, residue_3_letter_code, chain name, 
								#residue number, coordinates
			pdb_to_list.append(['ATOM', (f_id), row['fragment_type'], aa_dict[r].upper(), 
								row['chain_id'], (f_id), (c[0]), (c[1]), (c[2])])
			
		for a in pdb_to_list:
			begin = f'ATOM  {a[1]:>5} {a[2]:<4} {a[3]:>3} {a[4]}{a[5]:>4}    '
			coords = f'{a[6]:>8.3f}{a[7]:>8.3f}{a[8]:>8.3f}'
			end = "  1.00  0.00           C  "
			print(begin+coords+end)
		print('new fragment\n')
		#sys.exit()
			
			
			
			#print(pdb_to_list)
	
	
	
