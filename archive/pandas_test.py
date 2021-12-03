#!/usr/bin/python3 

import argparse
import json
import pandas as pd
import sys

import Bio.Data.IUPACData as conv
from Bio.PDB import MMCIF2Dict
import Bio.PDB.MMCIFParser as mmcifparser
#from Bio.PDB.DSSP import DSSP as dssp
from molmass import Formula
import numpy as np

cifparser = mmcifparser()
aa_dict = (conv.protein_letters_3to1_extended)
def make_atom_frame(files):
	assert(type(files) == list)
	"""
	- pdb ids
	- 	model ids
	- 		chain ids
	- 			residue ids # integers?
	- 			res code
	- 			atom code
	- 			xyz
	- 			missing
	-          max 3,2,1 turn 3,2,1
	"""
	for cif in files:
		psplit = cif.split('/')
		xyz  = cifparser.get_structure(psplit[-1][:-4], cif)
		missing = dict()
		cifdic = MMCIF2Dict.MMCIF2Dict(cif)
		unobs_keys = [k for k in cifdic.keys() if 'unobs' in k]
		tot_skip = []
		if len(unobs_keys) > 0:
			if '_pdbx_unobs_or_zero_occ_residues.label_seq_id' in cifdic:
				tot_skip.append(cifdic['_pdbx_unobs_or_zero_occ_residues.auth_seq_id'])
			
			if '_pdbx_unobs_or_zero_occ_atoms.label_seq_id' in cifdic:
				tot_skip.append(cifdic['_pdbx_unobs_or_zero_occ_atoms.authq_seq_id'])
		
		for i, m in enumerate(xyz.get_models()):
			for j, c in enumerate(m.get_chains()):
				index = 0
				for r in c.get_residues():
					
					#sys.exit()
					index += 1
					res = r.get_resname()
					res = res.capitalize()
					#skipping non-residues
					if res not in aa_dict.keys(): continue
					print(r)
					#print(dir(r))
					print(r.get_full_id())
					print(r.id)
					for atom in r.get_atoms():
						print(atom.get_full_id())
						continue
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
		sys.exit()
	
	df = pd.DataFrame(seq_data, columns =['Index', 'Molecule_Name','Model_ID',
										  'Chain_ID', 'Residue_ID', 'Residue',
										  'Atom', 'Type', 'X', 'Y', 'Z'])

parser = argparse.ArgumentParser(description='testing pandas multi-index')
parser.add_argument(
	'--file', required=False, type=str, metavar='<path>', 
	help='path to cif file')
parser.add_argument(
	'--cifs', '-c', required=False, type=str, metavar='<path>', 
	help='path to file containing set of cifs to build')

arg = parser.parse_args()

if arg.file:
	atm_frame = make_atom_frame([arg.file])
	print(atm_frame.head(1000000))	
elif arg.cifs:
	with open(arg.cifs, 'r') as fp:
		cifs = json.load(fp)
	fp.close()
	
	atm_frame = make_atom_frame(cifs)
	print(atm_frame.head(10000000))
