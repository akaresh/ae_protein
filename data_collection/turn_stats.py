#!/usr/bin/python3 

import argparse
import json
import os
import sys

alphabet = ["G", "P", "D", "E", "K", "R", "H", "S", "T", "N", "Q", "A", "M", "Y", "W", "V", "I", "L", "F", "C"]

parser = argparse.ArgumentParser(description='compute turn statistics')
parser.add_argument('--pdbs', required=True, type=str,
	metavar='<path>', help='path to chains directory')

arg = parser.parse_args()


aa_freq      = dict()
len_freq     = dict()
tot_turns    = 0
tot_chains = 0
for file in os.listdir(arg.pdbs):
	if file.endswith('.json'):
		input_json = os.path.join(arg.pdbs, file)
		
		with open(input_json, 'r') as fp:
			turns = json.load(fp)
		fp.close()
		#print(json.dumps(turns, indent=2))
		tot_chains += 1
		for id, turn in turns.items():
			#print(turn)
			skip = False
			for res in turn:
				ind, aa = res
				if aa not in alphabet:
					#print('bad aa')print(turn)
					skip = True
					break
			if skip: continue
			else:
				tot_turns += 1
				for res in turn:
					ind, aa = res
					if aa not in aa_freq: aa_freq[aa] = 0
					aa_freq[aa] += 1
				
				if len(turn) not in len_freq: len_freq[len(turn)] = 0
				len_freq[len(turn)] += 1

print(f'Chains with turns: {tot_chains}')
print(f'Total # of turns: {tot_turns}')
print(f'Turns per chain: {tot_turns / tot_chains:.2f}')
print()
print('Length distribution')
print('-------------------')
for l in sorted(len_freq):
	print(f'{l} : {len_freq[l]}')
print()
print('AA distribution')
print('---------------')
for aa in alphabet:
	print(f'{aa} : {aa_freq[aa]}')
"""
-1. turns per protein
0. total turns
1. length distribution
2. Amino acid distribution
"""