import sys

import numpy as np
import pandas as pd
import Bio.Data.IUPACData as conv
import Bio.PDB.MMCIFParser as mmcifparser
from molmass import Formula
#import Bio.PDB.internal_coords as internal_coords

cifparser = mmcifparser()
aa_dict = (conv.protein_letters_3to1_extended)


"""
1. just ignore structures with any missing res and atoms
2. just grab turns, before i was making fragments through the whole structure
3. just grab xyz of CA cluster the vectors of positions in the turns
	CA coordinates of the turns
4. install dssp
5. make pdb
5. run dssp on the turns
6. MLP on turns

learning
CNN autoencoder
distance matrix between all the atoms
M, MxM 
CNN down to some latent space
convolution transpose 
"""


def make_atom_frame(files):
	assert(type(files) == list)

	seq_data = []

	for cif in files:
		psplit = cif.split('/')
		xyz  = cifparser.get_structure(psplit[-1][:-4], cif)
		index = 0

		# if ftype == 'bbcen' or ftype == 'bbsc':
		# 	build C-beta's for glycine
		# 	internal_coords.IC_Residue.gly_Cbeta = True
		# 	#xyz.internal_to_atom_coordinates()
		# 	xyz.atom_to_internal_coordinates()

		for i, m in enumerate(xyz.get_models()):
			for j, c in enumerate(m.get_chains()):
				for r in c.get_residues():
					index += 1
					res = r.get_resname()
					res = res.capitalize()
					#skipping non-residues
					if res not in aa_dict.keys(): continue
					for atom in r.get_atoms():
						seq_data.append((index,
										 atom.get_full_id()[0],
										 atom.get_full_id()[1],
										 atom.get_full_id()[2],
										 atom.get_full_id()[3][1],
										 aa_dict.get(res),
										 atom.get_id(),
										 atom.element,
										 atom.get_coord()[0],
										 atom.get_coord()[1],
										 atom.get_coord()[2]))
	
	df = pd.DataFrame(seq_data, columns =['Index', 'Molecule_Name','Model_ID',
										  'Chain_ID', 'Residue_ID', 'Residue',
										  'Atom', 'Type', 'X', 'Y', 'Z'])
	
	return df


def make_CAframe(atomdf, size):
	new = []
	start = True
	
	fidxs    = None
	pos      = None
	frag_seq = None
	
	for counter, (idx, row) in enumerate(atomdf.iterrows()):
		if counter % 10000 == 0: print(counter)
		if row.Atom != 'CA': continue
		
		dic = {'pdb_id':row.Molecule_Name, 'model_id':row.Model_ID,
			'chain_id':row.Chain_ID, 'fragment_ids':None,
			'fragment_seq':None, 'xyz_set':None, 'fragment_type':'CA'}
		
		pdbid = row.Molecule_Name
		resid = row.Index
		chid  = row.Chain_ID
		
		if fidxs is None and pos is None and frag_seq is None:
			fidxs = [resid]
			pos = [[row.X, row.Y, row.Z]]
			frag_seq = str(row.Residue)
		
		# print(fidxs, row.Residue)
		
		if start:
			skip  = False
			start = False
			for i in range(1,size):
				df_row = atomdf[atomdf['Molecule_Name'] == pdbid]
				df_row = df_row[df_row['Model_ID'] == dic['model_id']]
				df_row = df_row[df_row['Chain_ID'] == chid]
				df_row = df_row[df_row['Index'] == (fidxs[i-1]+1)]
				df_row = df_row[df_row['Atom'] == 'CA']
				
				if df_row.empty:
					skip = True
					break
				fidxs.append(df_row['Index'].values[0])
				pos.append([df_row.X.values[0],
							df_row.Y.values[0],
							df_row.Z.values[0]])
				frag_seq += str(df_row['Residue'].values[0])
			
			if skip:
				fidxs    = None
				pos      = None
				frag_seq = None
				start    = True
				continue
		else:
			fidxs    = fidxs[1:]
			pos      = pos[1:]
			frag_seq = frag_seq[1:]
			
			df_row = atomdf[atomdf['Molecule_Name'] == pdbid]
			df_row = df_row[df_row['Index'] == (fidxs[-1]+1)]
			df_row = df_row[df_row['Chain_ID'] == chid]
			df_row = df_row[df_row['Atom'] == 'CA']
			
			if df_row.empty:
				fidxs    = None
				pos      = None
				frag_seq = None
				start    = True			
				continue
			
			fidxs.append(df_row['Index'].values[0])
			pos.append([
				df_row.X.values[0],
				df_row.Y.values[0],
				df_row.Z.values[0]])
			frag_seq += str(df_row['Residue'].values[0])
		
		dic['fragment_ids'] = fidxs
		dic['fragment_seq'] = frag_seq
		dic['xyz_set'] = pos
		new.append(dic)
			
	df = pd.DataFrame(new, columns =['pdb_id', 'model_id',
									 'chain_id', 'fragment_ids',
									 'fragment_seq', 'xyz_set',
									 'fragment_type'])
	return df


def make_bbframe(atomdf, size):
	new = []
	backbone = ['N', 'CA', 'C', 'O']
	
	start    = True
	fidxs    = []
	pos      = []
	frag_seq = ''	
	
	for idx, row in atomdf.iterrows():
		if row.Atom != 'N': continue
		
		dic = {'pdb_id':row.Molecule_Name, 'model_id':row.Model_ID,
			'chain_id':row.Chain_ID, 'fragment_ids':None,
			'fragment_seq':None, 'xyz_set':None, 'fragment_type':'bb'}
		
		pdbid = row.Molecule_Name
		resid = row.Index
		chid  = row.Chain_ID
		
		if start:
			start = False
			skip  = False
			for i in range(0,size):
				df_row = atomdf[atomdf['Molecule_Name'] == pdbid]
				df_row = df_row[df_row['Index'] == (resid+i)]
				df_row = df_row[df_row['Chain_ID'] == chid]
				
				if df_row.empty:
					skip = True
					break
				
				fidxs.append(df_row['Index'].values[0])
				frag_seq += str(df_row['Residue'].values[0])
				
				for b in backbone:
					df_row1 = df_row[df_row['Atom'] == b]
					if df_row1.empty: raise(f'df_row1 empty at {b}')
					pos.append([df_row1.X.values[0],
								df_row1.Y.values[0],
								df_row1.Z.values[0]])
			
			if skip:
				fidxs    = []
				pos      = []
				frag_seq = ''
				start    = True
				continue
		else:
			fidxs    = fidxs[1:]
			pos      = pos[4:]
			frag_seq = frag_seq[1:]
			
			df_row = atomdf[atomdf['Molecule_Name'] == pdbid]
			df_row = df_row[df_row['Index'] == (fidxs[-1]+1)]
			df_row = df_row[df_row['Chain_ID'] == chid]
			
			if df_row.empty:
				fidxs    = []
				pos      = []
				frag_seq = ''
				start    = True			
				continue
			
			fidxs.append(df_row['Index'].values[0])
			frag_seq += str(df_row['Residue'].values[0])
			for b in backbone:
				df_row1 = df_row[df_row['Atom'] == b]
				if df_row1.empty: raise(f'df_row1 empty at {b}')
				pos.append([
					df_row1.X.values[0],
					df_row1.Y.values[0],
					df_row1.Z.values[0]])
		
		dic['fragment_ids'] = fidxs
		dic['fragment_seq'] = frag_seq
		dic['xyz_set'] = pos
		new.append(dic)
	
	df = pd.DataFrame(new, columns =['pdb_id','model_id',
									 'chain_id', 'fragment_ids',
									 'fragment_seq', 'xyz_set',
									 'fragment_type'])	
	return df


def res_cen(dict_positions):
	#print(dict_positions)
	total_mass = 0
	mx = 0
	my = 0
	mz = 0
	for keys, values in dict_positions.items():
		total_mass += (Formula(values[-1]).mass)
		mx += ((Formula(values[-1]).mass)*values[0])
		my += ((Formula(values[-1]).mass)*values[1])
		mz += ((Formula(values[-1]).mass)*values[2])
	return [mx/total_mass, my/total_mass, mz/total_mass]


# def make_bbcen(atomdf, size):
# 	backbone = ['N', 'CA', 'C', 'O']
# 	new = []
# 	
# 	start    = True
# 	fidxs    = []
# 	pos      = []
# 	frag_seq = ''
# 	
# 	for idx, row in atomdf.iterrows():
# 		if row.Atom != 'N': continue
# 		
# 		dic = {'pdb_id':row.Molecule_Name, 'model_id':row.Model_ID,
# 			'chain_id':row.Chain_ID, 'fragment_ids':None,
# 			'fragment_seq':None, 'xyz_set':None,
# 			'fragment_type':'bb+cen'}
# 		resid = row.Index
# 		chid  = row.Chain_ID
# 		
# 		fidxs = []
# 		pos = []
# 		frag_seq = ''
# 		atoms = []
# 		
# 		if start:
# 			start = False
# 			skip  = False
# 			for i in range(0, size):
# 				df_row = atomdf[atomdf['Molecule_Name'] == pdbid]
# 				df_row = atomdf[atomdf['Index'] == (resid + i)]
# 				df_row = df_row[df_row['Chain_ID'] == chid]
# 				
# 				if df_row.empty:
# 					skip = True
# 					break
# 				
# 				if df_row['Residue'] == 'G': 
# 					skip = True
# 					break
# 				
# 				fidxs.append(df_row['Index'].values[0])
# 				frag_seq += str(df_row['Residue'].values[0])
# 				
# 				res_atoms_pos = {}
# 				atoms_sublist = []
# 				for sub_id, sub_row in df_row.iterrows():
# 					if sub_row.Atom in backbone: continue
# #						atoms_sublist.append(sub_row.Atom)
# #						pos.append([sub_row.X, sub_row.Y, sub_row.Z])
# 					else:
# 						if sub_row.Atom in res_atoms_pos.keys():
# 							raise ValueError('Atom Duplication')
# 						res_atoms_pos[sub_row.Atom] = [sub_row.X,
# 													   sub_row.Y,
# 													   sub_row.Z,
# 													   sub_row.Type]
# 				
# 				centroid_pos = res_cen(res_atoms_pos)
# 				
# 				for b in backbone:
# 					df_row1 = df_row[df_row['Atom'] == b]
# 					if df_row1.empty: raise(f'df_row1 empty at {b}')
# 					pos.append([df_row1.X.values[0],
# 								df_row1.Y.values[0],
# 								df_row1.Z.values[0]])
# 				
# 				
# 				
# 			if skip:
# 				fidxs    = []
# 				pos      = []
# 				frag_seq = ''
# 				start    = True
# 				continue
# 			
# 			res_atoms_pos = {}
# 			atoms_sublist = []
# 			for sub_id, sub_row in df_row.iterrows():
# 				if sub_row.Atom in backbone: continue
# #					atoms_sublist.append(sub_row.Atom)
# #					pos.append([sub_row.X, sub_row.Y, sub_row.Z])
# 				else:
# 					if sub_row.Atom in res_atoms_pos.keys():
# 						raise ValueError('Atom Duplication')
# 					res_atoms_pos[sub_row.Atom] = [sub_row.X,
# 												   sub_row.Y,
# 												   sub_row.Z,
# 												   sub_row.Type]
# 			
# 			centroid_pos = res_cen(res_atoms_pos)
# 			
# 			df_row1 = df_row[df_row['Atom'] == b]
# 			if df_row1.empty: raise(f'df_row1 empty at {b}')
# 			pos.append([df_row1.X.values[0],
# 								df_row1.Y.values[0],
# 								df_row1.Z.values[0]])
# 			
# 			for b in backbone:
# 				df_row1 = df_row[df_row['Atom'] == b]
# 			
# 			atoms_sublist.append('Centroid')
# 			atoms.append(atoms_sublist)
# 
# 			if len(res_atoms_pos) > 0:
# 				#print(res_atoms_pos)
# 				pos.append(res_cen(res_atoms_pos))
# 				#print(res_cen(res_atoms_pos))
# 			else:
# 				#print(res_cen)
# 				#print('glycine') = glycine extention
# 				pass
# 
# 		if skip: continue
# 
# 		dic['fragment_ids'] = fidxs
# 		dic['fragment_seq'] = frag_seq
# 		dic['xyz_set'] = pos
# 		dic['atoms_list'] = atoms
# 		new.append(dic)
# 	
# 	df = pd.DataFrame(new, columns =['pdb_id', 'model_id',
# 									 'chain_id', 'fragment_ids',
# 									 'fragment_seq', 'xyz_set',
# 									 'fragment_type', 'atoms_list'])	
# 	return df


def make_fragment_frame(atomdf, size, ftype=None):
	assert(ftype != None)
	
	if   ftype == 'CA':    return make_CAframe(atomdf, size)
	elif ftype == 'bb':    return make_bbframe(atomdf, size)
	elif ftype == 'bbcen': return make_bbcen(atomdf, size)
	elif ftype == 'bbsc':  return make_bbsc(atomdf, size)
