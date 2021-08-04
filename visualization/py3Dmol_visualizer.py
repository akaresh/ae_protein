#Note
#1. Alias for training_tools

###Plan###
'''
1. Visualisation before the model
2. Visualisation after the normalisation
3. Visualisation after the model
'''


def pdb_gen(row, dtype):
	#change dtype later
	#converting to 3 letter code for pdb file
	aa_dict = conv.protein_letters_1to3_extended

	pdb_to_list = []

	if dtype == 'xyz':
		print('initial')

		for f_id, r, c in zip(row['fragment_ids'], row['fragment_seq'], row['xyz_set']):
								#atom, fragment_id, atom_type, residue_3_letter_code, chain name,
								#residue number, coordinates
			#print(c)
			pdb_to_list.append(['ATOM', (f_id), row['fragment_type'], "XXX",
								row['chain_id'], (f_id), (c[0]), (c[1]), (c[2])])

	elif dtype == 'normalized':
		print('normalized')

		for f_id, r, c in zip(row['fragment_ids'], row['fragment_seq'], [row['norm_frag'][i:i+3] for i in
			                  range(0, len(row['norm_frag']), 3)]):
								#atom, fragment_id, atom_type, residue_3_letter_code, chain name,
								#residue number, coordinates
			pdb_to_list.append(['ATOM', (f_id), row['fragment_type'], aa_dict[r].upper(),
								row['chain_id'], (f_id), (c[0]), (c[1]), (c[2])])
	elif dtype == 'cnn':
		pass

	for a in pdb_to_list:
		begin = f'ATOM  {a[1]:>5} {a[2]:<4} {a[3]:>3} {a[4]}{a[5]:>4}    '
		coords = f'{a[6]:>8.3f}{a[7]:>8.3f}{a[8]:>8.3f}'
		end = "  1.00  0.00           C  "
		print(begin+coords+end)

	return pdb_to_list

def pdb_writer(atoms=None, seq=None, chain=None, coords=None):
	aa_dict = conv.protein_letters_1to3_extended

	###asserts and checks
	#check atoms, seq, chain, coords
	#if seq is type str, length < 4

	res_id = 0
	prev = None
	for i, (a, s, ch, coo) in enumerate(zip(atoms, seq, chain, coords)):
		if s != prev: res_id += 1
		x, y, z = coo

		begin = f'ATOM  {i+1:>5} {a:<4} {aa_dict[s].upper():>3} {ch}{res_id:>4}    '
		coordinates = f'{x:>8.3f}{y:>8.3f}{z:>8.3f}'
		end = "  1.00  0.00           C  "
		print(begin+coordinates+end)
	return True



#!/usr/bin/python3
if __name__ == '__main__':
	import sys
	import argparse


	#https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/bbPlane.py
	# import pymol
	# from pymol import cmd
	# from pymol import stored

	import pandas as pd
	import sys

	import Bio.Data.IUPACData as conv

	from training_tools import normalize_frag
	import architectures


	#https://biopython.org/docs/1.75/api/Bio.PDB.HSExposure.html?highlight=write%20pdb%20file

	#Bio.PDB.HSExposure.HSExposureCA


	parser = argparse.ArgumentParser(
		description='desc')
	parser.add_argument(
		'--df', '-d', required=True, type=str,
		metavar='<path>', help='path to xz compressed pandas dataframe')

	arg = parser.parse_args()

	df = pd.read_pickle(arg.df, compression='xz').sample(
		frac=1.0, random_state=42)

	df['norm_frag'] = df.xyz_set.apply(normalize_frag)

	# print(df.head(3))
	# print(df.columns)
	# print(df.shape)
	# print(df['xyz_set'])
	# print(df['norm_frag'])
	# print()
	#sys.exit()

	pdb_writer(atoms=['CA']*7, seq=df.fragment_seq[1], chain=['A']*7, coords=df.xyz_set[0])
	sys.exit()

	for idx, frag_row in df.iterrows():
		print(frag_row)
		initial = pdb_gen(frag_row, 'xyz')
		#normalized_frag = pdb_gen(frag_row, 'normalized')
		#after_cnn = pdb_gen(frag_row, 'cnn')
		sys.exit()

	###code for test examples




