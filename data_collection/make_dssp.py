#!/usr/bin/python3 

import argparse
import os
import sys

parser = argparse.ArgumentParser(description='creating dssp files')
parser.add_argument('--pdbs', required=True, type=str,
	metavar='<path>', help='path to chains directory')

arg = parser.parse_args()

for pdb in os.listdir(arg.pdbs):
	if os.path.isfile(os.path.join(arg.pdbs, pdb)):
		if pdb.endswith('.pdb'):
			splits = pdb.split('.')
			pdbid = splits[0]
			#print(pdbid)
			input_pdb = os.path.join(arg.pdbs, pdb)
			output_dssp = os.path.join(arg.pdbs, f'{pdbid}.dssp')
			if not os.path.isfile(output_dssp):
				cmd = f'mkdssp -i {input_pdb} -o {output_dssp}'
				print(cmd)
				os.system(cmd)
				#sys.exit()
