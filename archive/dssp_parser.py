#!/usr/bin/python3 

import argparse
import json
import os
import sys

from Bio.PDB import MMCIF2Dict

parser = argparse.ArgumentParser(description='creating dssp dictionaries')
parser.add_argument('--pdbs', required=True, type=str,
	metavar='<path>', help='path to chains directory')

arg = parser.parse_args()

for pdb in os.listdir(arg.pdbs):
	if os.path.isfile(os.path.join(arg.pdbs, pdb)):
		if pdb.endswith('.dssp'):
			input_dssp = os.path.join(arg.pdbs, pdb)
			#print(input_dssp)
			print(pdb)
			ids = pdb.split('.')
			turn_defs = ids[0]+'.turns.json'
			turn_json = os.path.join(arg.pdbs, turn_defs)
			if os.path.isfile(turn_json): continue
			
			with open(input_dssp, 'r') as fp:
				lines = fp.readlines()
			
			turn_id = 0
			start = False
			begin = False
			turns = dict()
			for line in lines:
				if '#  RESIDUE AA' in line:
					begin = True
					continue
				if begin:
					line = line.rstrip()
					if start:
						if line[16] != 'T':
							start = False
							if len(turns[turn_id]) == 1:
								del turns[turn_id]
								turn_id -= 1
						elif line[16] == 'T':
							row = line.split()
							turns[turn_id].append((row[1], line[13]))
					else:
						if line[16] == 'T':
							start = True
							turn_id += 1
							row = line.split()
							turns[turn_id] = list()
							turns[turn_id].append((row[1], line[13]))
			
			turn_defs = ids[0]+'.turns.json'
			turn_json = os.path.join(arg.pdbs, turn_defs)
			with open(turn_json, 'w') as fp:
				json.dump(turns, fp)
			fp.close()
		else: continue