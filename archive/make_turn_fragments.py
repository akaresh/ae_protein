#!/usr/bin/python3 

import argparse
import json
import os
import sys

import Bio.PDB as bpdb
import pandas as pd

alphabet = ["G", "P", "D", "E", "K", "R", "H", "S", "T", "N", "Q", "A", "M", "Y", "W", "V", "I", "L", "F", "C"]

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


parser = argparse.ArgumentParser(description='make turn fragments')
parser.add_argument('--pdbs', required=True, type=str,
	metavar='<path>', help='path to chains directory')
parser.add_argument('--flank', required=True, type=int,
	metavar='<int>', help='number of flanking residues to take')

arg = parser.parse_args()

frag_set = []
for file in os.listdir(arg.pdbs):
	if file.endswith('.json'):
		input_json = os.path.join(arg.pdbs, file)
		info = file.split('.')
		pdbid = info[0]
		chainid = info[0][4:].capitalize()
		print(pdbid)
		#print(input_json)
		pdb = pdbid+'.pdb'
		input_pdb = os.path.join(arg.pdbs, pdb)
		if os.path.isfile(input_pdb):
			turns = dict()
			with open(input_json, 'r') as fp:
				turns = json.load(fp)
			fp.close()
			
			parser = bpdb.PDBParser()
			try:
				structure = parser.get_structure(pdbid, input_pdb)
			except:
				continue
			
			if chainid not in structure[0]:
				chainid = chainid.lower()
				if chainid not in structure[0]:
					continue
			
			for id, turn in turns.items():
				if len(turn) > 5: continue
				
				res = [r[1] for r in turn if r[1] not in alphabet]
				ids = [r[0] for r in turn]
				if len(res) > 0: continue
				
				turn_aa = dict()
				for id, r in zip(ids, res): 
					if id not in turn_aa: turn_aa[int(id)] = r
				
				try:
					ids = [int(ind) for ind in ids]
					ids = sorted(ids)
				except:
					continue
				
				frag_ids = list(range(
					ids[0] - arg.flank, ids[-1] + arg.flank + 1))
				
				frag_dict = {
					'pdbid' : pdbid,
					'seq' : None,
					'ids' : frag_ids,
					'xyz' : None
				}
				seq = ''	
				frag_xyz = []
				skip = False
				for fid in frag_ids:
					if fid not in structure[0][chainid]:
						skip = True
						break
					
					if "CA" not in structure[0][chainid][fid]:
						skip = True
						break
					
					residue = structure[0][chainid][fid]
					if residue.get_resname().lower() not in aa_dict:
						skip = True
						break
					residue_name = aa_dict[residue.get_resname().lower()]
					if fid in turn_aa: assert(turn_aa[fid] == residue_name)
					seq += residue_name
					
					grp = residue["CA"].get_coord().tolist()
					frag_xyz.append(grp)
				
				if not skip:
					frag_dict['seq'] = seq
					frag_dict['xyz'] = frag_xyz
					frag_set.append(frag_dict)
				else: continue

df = pd.DataFrame(frag_set)	
print(df.head(5))
print(df.columns)
print(df.shape)

outfrag = 'turnframe.pickle.xz'
df.to_pickle(outfrag, compression='xz')




				
				
			
			