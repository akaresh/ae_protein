#!/usr/bin/env python

import argparse
import json
import lzma
import os
import pickle
import sys

import Bio.PDB as bpdb
from fragment_generators import ca_fragment, bb_fragment, bbcen_fragment
from fragment_generators import bbsc_fragment
from fragment_generators import alphabet, aa_dict, aa_sc_atoms
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from turns_data_class import TurnFrame


def get_frag_xyz(struct, chain, ids, dssp, frag_type):
	if frag_type == 'ca':     return ca_fragment(struct, chain, ids, dssp)
	elif frag_type == 'bb':   return bb_fragment(struct, chain, ids, dssp)
	elif frag_type =='bbcen': return bbcen_fragment(struct, chain, ids, dssp)
	elif frag_type =='bbsc':  return bbsc_fragment(struct, chain, ids, dssp)
	else: return None


parser = argparse.ArgumentParser(description='make turn fragments')
parser.add_argument('--pdbs', required=True, type=str,
	metavar='<path>', help='path to chains directory')
parser.add_argument('--flank', required=True, type=int,
	metavar='<int>', help='number of flanking residues to take')
parser.add_argument('--type', required=False, type=str,
	default='ca', metavar='<str>', help='fragment type')
parser.add_argument('--save', required=True, type=str,
	default=None, metavar='<path>',
	help='folder and base name for turn fragment data set, dont add metadata')

arg = parser.parse_args()

assert(arg.type == 'ca' or arg.type == 'bb' or arg.type == 'bbcen'
	   or arg.type == 'bbsc')

frag_set = []
for file in os.listdir(arg.pdbs):
	if file.endswith('.turns.json'):
		# collect turns.json and dssp.json files for this structure
		info = file.split('.')
		pdbid = info[0]
		chainid = info[0][4:].capitalize()
		print(pdbid)
		
		pdb = pdbid+'.pdb'
		input_pdb = os.path.join(arg.pdbs, pdb)
		assert(os.path.isfile(input_pdb))
		
		dssp_json = os.path.join(arg.pdbs, pdbid+'.dssp.json')
		if not os.path.isfile(dssp_json): continue
		
		turns_json = os.path.join(arg.pdbs, file)
		turns = dict()
		with open(turns_json, 'r') as fp:
			turns = json.load(fp)
		
		ssdic = dict()
		with open(dssp_json, 'r') as fp:
			ssdic = json.load(fp)
		
		# read in pdb structures with PDBParser
		parser = bpdb.PDBParser()
		structure = parser.get_structure(pdbid, input_pdb)
			
		if chainid not in structure[0]:
			chainid = chainid.lower()
			if chainid not in structure[0]:
				continue
		
		for id, turn in turns.items():
			
			if turn['res1'] - arg.flank < 0: continue
			if str(turn['res4'] + arg.flank) not in ssdic: continue
			
			frag_ids = list(range(
				turn['res1'] - arg.flank,
				turn['res4'] + arg.flank + 1))
			
			seq = ''
			skip = False
			for fi in frag_ids:
				if str(fi) not in ssdic:
					skip = True
					break
			if skip: continue
			
			xyz = get_frag_xyz(structure, chainid, frag_ids, ssdic, arg.type)
			print(xyz)
			if xyz is None: continue
			print(xyz)
			seq    = ''
			secstr = ''
			for fi in frag_ids:
				residue = structure[0][chainid][fi]
				residue_name = aa_dict[residue.get_resname().lower()]
				seq += residue_name
				secstr += ssdic[str(fi)][1]
			
			mat = cdist(xyz, xyz, metric='euclidean')
			mat = np.reshape(mat, (1, mat.shape[0], mat.shape[1]))
			mat = np.divide(mat, 7.0).tolist()
			
			frag_dic = {
				'pdbid' : pdbid,
				'res_ids' : frag_ids,
				'seq' : seq,
				'sec_struct' : secstr,
				'coords' : xyz,
				'dmatrix' : mat,
				'turn_type' : turn['type1'],
				'turn_symbol' : turn['symbol']
			}
			frag_set.append(frag_dic)

df = pd.DataFrame(frag_set)	
print(df.head(5))
print(df.columns)
print(df.shape)

new_turn_frame = TurnFrame(df)
new_turn_frame.flank = arg.flank
new_turn_frame.turn_type = arg.type
new_turn_frame.split = None

save_dir = arg.save.split('/')
save_dir = '/'.join(save_dir[:-1])
assert(os.path.isdir(save_dir))

save_name = f"{arg.save}_full_flank{arg.flank}_type{arg.type}.pickle.xz"
print(save_name)
with lzma.open(save_name, 'wb') as fp:
	pickle.dump(new_turn_frame, fp)
