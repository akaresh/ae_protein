import sys

import pandas as pd
import Bio.Data.IUPACData as conv
import Bio.PDB.MMCIFParser as mmcifparser
from molmass import Formula

parser = mmcifparser()
aa_dict = (conv.protein_letters_3to1_extended)

def make_atom_frame(files):
	assert(type(files) == list)
	
	seq_data = []
	
	for cif in files:
		psplit = cif.split('/')
		xyz  = parser.get_structure(psplit[-1][:-4], cif)
		index = 0
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
	for counter, (idx, row) in enumerate(atomdf.iterrows()):
		if counter % 10000 == 0: print(counter)
		if row.Atom != 'CA': continue
		
		dic = {'pdb_id':row.Molecule_Name, 'model_id':row.Model_ID,
			'chain_id':row.Chain_ID, 'fragment_ids':None,
			'fragment_seq':None, 'xyz_set':None, 'fragment_type':'CA'}
		
		pdbid = row.Molecule_Name
		resid = row.Index
		chid  = row.Chain_ID
		fidxs = [resid]
		#print(fidxs, row.Residue)
		pos = [[row.X, row.Y, row.Z]]
		frag_seq = str(row.Residue)
		skip = False

		for i in range(1,size):
			df_row = atomdf[atomdf['Molecule_Name'] == pdbid]
			df_row = df_row[df_row['Index'] == (fidxs[i-1]+1)]
			df_row = df_row[df_row['Chain_ID'] == chid]
			df_row = df_row[df_row['Atom'] == 'CA']

			if df_row.empty:
				skip = True
				break
			fidxs.append(df_row['Index'].values[0])
			pos.append([df_row.X.values[0], 
						df_row.Y.values[0], 
						df_row.Z.values[0]])
			frag_seq += str(df_row['Residue'].values[0])

		if skip: continue

		dic['fragment_ids'] = fidxs
		dic['fragment_seq'] = frag_seq
		dic['xyz_set'] = pos
		new.append(dic)
	
	df = pd.DataFrame(new, columns =['pdb_id', 'model_id',
									 'chain_id', 'fragment_ids',
									 'fragment_seq', 'xyz_set',
									 'fragment_type'])
	return df
#new = make_CAframe(df, k)


def make_bbframe(atomdf, size):
	new = []
	#size = 3
	backbone = ['N', 'CA', 'C', 'O']
	for idx, row in atomdf.iterrows():
		if row.Atom != 'N': continue

		dic = {'pdb_id':row.Molecule_Name, 'model_id':row.Model_ID,
			'chain_id':row.Chain_ID, 'fragment_ids':None,
			'fragment_seq':None, 'xyz_set':None, 'fragment_type':'bb'}
		resid = row.Index
		chid  = row.Chain_ID
		fidxs = []
		pos = []
		frag_seq = ''

		skip = False
		for i in range(0,size):
			df_row = atomdf[atomdf['Index'] == (resid + i)]
			df_row = df_row[df_row['Chain_ID'] == chid]

			if df_row.empty:
				skip = True
				break
			fidxs.append(df_row['Index'].values[0])
			frag_seq += str(df_row['Residue'].values[0])

			for b in backbone:
				df_row1 = df_row[df_row['Atom'] == b]
				pos.append([df_row1.X.values[0],
							df_row1.Y.values[0],
							df_row1.Z.values[0]])

		if skip: continue

		dic['fragment_ids'] = fidxs
		dic['fragment_seq'] = frag_seq
		dic['xyz_set'] = pos
		new.append(dic)
	
	df = pd.DataFrame(new, columns =['pdb_id','model_id',
									 'chain_id', 'fragment_ids',
									 'fragment_seq', 'xyz_set',
									 'fragment_type'])	
	return df

# new1 = make_bbframe(df, 5)
# print(new1)

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

def make_bbcen(atomdf, size):
	backbone = ['N', 'CA', 'C', 'O']
	new3 = []
	for idx, row in atomdf.iterrows():
		if row.Atom != 'N': continue

		dic = {'pdb_id':row.Molecule_Name, 'model_id':row.Model_ID,
			'chain_id':row.Chain_ID, 'fragment_ids':None,
			'fragment_seq':None, 'xyz_set':None, 'atoms_list':None,
			'fragment_type':'bb+cen'}
		resid = row.Index
		chid  = row.Chain_ID
		fidxs = []
		pos = []
		frag_seq = ''
		atoms = []

		skip = False

		for i in range(0,size):
			df_row = atomdf[atomdf['Index'] == (resid + i)]
			df_row = df_row[df_row['Chain_ID'] == chid]

			if df_row.empty:
				skip = True
				break
			fidxs.append(df_row['Index'].values[0])
			frag_seq += str(df_row['Residue'].values[0])

			res_atoms_pos = {}
			atoms_sublist = []
			for sub_id, sub_row in df_row.iterrows():
				#print(sub_row.Atom)
				if sub_row.Atom in backbone:
					atoms_sublist.append(sub_row.Atom)
					pos.append([sub_row.X, sub_row.Y, sub_row.Z])
				else:
					if sub_row.Atom in res_atoms_pos.keys():
						raise ValueError('Atom Duplication')
					res_atoms_pos[sub_row.Atom] = [sub_row.X,
												   sub_row.Y,
												   sub_row.Z,
												   sub_row.Type]
			atoms_sublist.append('Centroid')
			atoms.append(atoms_sublist)

			if len(res_atoms_pos) > 0:
				#print(res_atoms_pos)
				pos.append(res_cen(res_atoms_pos))
				#print(res_cen(res_atoms_pos))
			else:
				#print('glycine') = glycine extention
				pass

		if skip: continue

		dic['fragment_ids'] = fidxs
		dic['fragment_seq'] = frag_seq
		dic['xyz_set'] = pos
		dic['atoms_list'] = atoms
		new3.append(dic)
	
	df = pd.DataFrame(new, columns =['pdb_id', 'model_id',
									 'chain_id', 'fragment_ids',
									 'fragment_seq', 'xyz_set',
									 'fragment_type', 'atoms_list'])	
	return df

# a = (make_bbcen(df, 4))

def make_fragment_frame(atomdf, size, ftype=None):
	assert(ftype != None)
	
	if   ftype == 'CA':    return make_CAframe(atomdf, size)
	elif ftype == 'bb':    return make_bbframe(atomdf, size)
	elif ftype == 'bbcen': return make_bbcen(atomdf, size)
	elif ftype == 'bbsc':  return make_bbsc(atomdf, size)
