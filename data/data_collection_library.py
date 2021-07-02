import sys

import pandas as pd
import Bio.Data.IUPACData as conv
import Bio.PDB.MMCIFParser as mmcifparser
from molmass import Formula

parser = mmcifparser()
aa_dict = (conv.protein_letters_3to1_extended)


xyz  = parser.get_structure()

seq_data = []
index = 0
for i, m in enumerate(xyz.get_models()):
	for j, c in enumerate(m.get_chains()):
		for r in c.get_residues():
			index += 1
			res = r.get_resname()
			res = res.capitalize()

			#skipping non-residues
			if res not in aa_dict.keys():
				continue

			for atom in r.get_atoms():

				seq_data.append((index, atom.get_full_id()[0], atom.get_full_id()[1], atom.get_full_id()[2],
								 atom.get_full_id()[3][1], aa_dict.get(res), atom.get_id(), atom.element,
								 atom.get_coord()[0], atom.get_coord()[1], atom.get_coord()[2]))

df = pd.DataFrame(seq_data, columns =['Index', 'Molecule_Name','Model_ID', 'Chain_ID', 'Residue_ID',
									  'Residue', 'Atom', 'Type', 'X', 'Y', 'Z'])
pd.set_option('max_columns', None)
print(df.head(10))
