#!/usr/bin/python3 

import argparse
import json
import os
import sys

import Bio.PDB as bpdb
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

alphabet = [
	"G",
	"P",
	"D",
	"E",
	"K",
	"R",
	"H",
	"S",
	"T",
	"N",
	"Q",
	"A",
	"M",
	"Y",
	"W",
	"V",
	"I",
	"L",
	"F",
	"C"
]

aa_dict = {
	"ala" : "A", 
	"arg" : "R",
	"asn" : "N",
	"asp" : "D",
	"cys" : "C",
	"gln" : "Q",
	"glu" : "E",
	"gly" : "G",
	"his" : "H",
	"ile" : "I",
	"leu" : "L",
	"lys" : "K",
	"met" : "M",
	"phe" : "F",
	"pro" : "P",
	"ser" : "S",
	"thr" : "T",
	"trp" : "W",
	"tyr" : "Y",
	"val" : "V"
}

aa_sc_atoms = {
	"A" : ["CB"], 
	"R" : ["CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"],
	"N" : ["CB", "CG", "OD1", "ND2"],
	"D" : ["CB", "CG", "OD1", "OD2"],
	"C" : ["CB", "SG"],
	"Q" : ["CB", "CG", "CD", "OE1", "NE2"],
	"E" : ["CB", "CG", "CD", "OE1", "OE2"],
	"G" : ["CA"],
	"H" : ["CB", "CG", "ND1", "CD2", "CE1", "NE2"],
	"I" : ["CB", "CG1", "CG2", "CD1"],
	"L" : ["CB", "CG", "CD1", "CD2"],
	"K" : ["CB", "CG", "CD", "CE", "NZ"],
	"M" : ["CB", "CG", "SD", "CE"],
	"F" : ["CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
	"P" : ["CB", "CG", "CD", "N", "CA"],
	"S" : ["CB", "OG"],
	"T" : ["CB", "OG1", "CG2"],
	"W" : ["CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"],
	"Y" : ["CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"],
	"V" : ["CB", "CG1", "CG2"]
}

def ca_fragment(struct, chain, ids, dssp):
	pos = []
	
	for i in ids:
		if i not in struct[0][chain]: return None
		if "CA" not in struct[0][chain][i]: return None
		
		residue = struct[0][chain][i]
		if residue.get_resname().lower() not in aa_dict: return None
		
		grp = residue["CA"].get_coord().tolist()
		pos.append(grp)
	
	return pos


def bb_fragment(struct, chain, ids, dssp):
	bbatms = ['N', 'CA', 'C', 'O']
	pos = []
	
	for i in ids:
		if i not in struct[0][chain]: return None
		
		residue = struct[0][chain][i]
		if residue.get_resname().lower() not in aa_dict: return None
		
		for bb in bbatms:
			if bb not in struct[0][chain][i]: return None
			
			grp = residue[bb].get_coord().tolist()
			pos.append(grp)
	
	return pos


def bbcen_fragment(struct, chain, ids, dssp):
	pass


def bbsc_fragment(struct, chain, ids, dssp):
	pass


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

arg = parser.parse_args()

assert(arg.type == 'ca' or arg.type == 'bb' or arg.type == 'bbcen'
	   or arg.type == 'bbsc')

frag_set = []
for file in os.listdir(arg.pdbs):
	if len(frag_set) >= 1000: break
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
		#try:
		structure = parser.get_structure(pdbid, input_pdb)
		#except:
		#	continue
			
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
			if xyz is None: continue
			
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
			print(json.dumps(mat, indent=2))
			sys.exit()
			frag_set.append(frag_dic)

df = pd.DataFrame(frag_set)	
print(df.head(5))
print(df.columns)
print(df.shape)


"""

id
pdbid
res-ids
seq
ss
coords
dmatrix
type

PISCES
50 % sequence id
~28k
100 turn type distribution
80/20 70/30
80 -> training
20 -> testing
classes have to be balanced

MSE / gt relative square difference

70 / 30

"""



























