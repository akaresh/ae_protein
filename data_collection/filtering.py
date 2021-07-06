import sys
import argparse
import os
import urllib


from gemmi import cif #reading cif files
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

parser = argparse.ArgumentParser(description='choosing complete pdb')
parser.add_argument('--file', required=True, type=str,
	metavar='<path>', help='path to text file with all of the pdbs')
parser.add_argument('--n', required=False, type=int, default = 10,
	metavar='<number>', help='number of pdbs for the dataset')
parser.add_argument('--name', required=False, type=str, default = 'output.txt',
	metavar='<string>', help='naming output txt file')
parser.add_argument('--fxn', required=False, type=str, default = 'CA',
	metavar='<string>', help='choosing between bbcen and CA')
arg = parser.parse_args()

#reading mmcif as dictionary
#mmdict = MMCIF2Dict()

#print(mmdict['_pdbx_unobs_or_zero_occ_atoms.label_atom_id'])
#print(mmdict['_pdbx_unobs_or_zero_occ_residues.label_seq_id'])


f = open(arg.file, 'r')
lines = f.readlines()
header = lines[0]
lines = lines[1:]

#writing a new file with filtered pdbs
# newf = open(arg.name, 'w')
updated= []

#for i in range(len(lines)):
	#extracting the pdb ids and the chain id
for line in lines:
	pdbid = line[:4]
	chain = line[4]
	print(pdbid, chain)
	if len(updated) <= arg.n-1:
		#extracting file info online
		with urllib.request.urlopen(f'https://files.rcsb.org/header/{pdbid}.cif') as c:
			doc = cif.read_string(c.read())
			block = doc.sole_block()

			#checking chain for missing residues:
			#for making_bb+cen
			#if arg.fxn == 'bbcen':
			res_chain = [i for i in block.find_loop( '_pdbx_unobs_or_zero_occ_residues.auth_asym_id')]

			if len(res_chain) == 0:

				atom_chain = [i for i in block.find_loop( '_pdbx_unobs_or_zero_occ_atoms.auth_asym_id')]
				ac = [element for element in block.find_loop('_pdbx_unobs_or_zero_occ_atoms.label_atom_id')]
				print(atom_chain, ac)
				for a, r in zip(atom_chain, ac):
					if r == 'CA' and r == chain:
						continue
				updated.append(pdbid+chain)

			else:
				#for rc in res_chain:
				if res_chain[0] == chain:
					#print('skip')
					continue
					#for future more extensive referencing
					#res = [element for element in block.find_loop('_pdbx_unobs_or_zero_occ_residues.label_seq_id')]


			# #checking for CA
			#if arg.fxn == 'CA':
			# atom_chain = [i for i in block.find_loop( '_pdbx_unobs_or_zero_occ_atoms.auth_asym_id')]
			# ac = [element for element in block.find_loop('_pdbx_unobs_or_zero_occ_atoms.label_atom_id')]
			# print(atom_chain, ac)
			# for a, r in zip(atom_chain, ac):
			# 	if r == 'CA' and r == chain:
			# 		continue
			# updated.append(pdbid+chain)
				#sys.exit()
				# for ac in atom_chain:
				# 	ac = [element for element in block.find_loop('_pdbx_unobs_or_zero_occ_atoms.label_atom_id')]
				#print(i)
print(updated)
print(len(updated))
	#sys.exit()
