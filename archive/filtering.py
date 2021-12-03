#!/usr/bin/python3 

import sys
import argparse
import os
import urllib.request
import json

from gemmi import cif #reading cif files

parser = argparse.ArgumentParser(description='choosing complete pdb')
parser.add_argument('--file', required=False, type=str,
	metavar='<path>', help='path to text file with all of the pdbs')
parser.add_argument('--dir', required=False, type=str,
	metavar='<path>', help='directory of cif files of chains')
parser.add_argument('--n', required=False, type=int, default = 10,
	metavar='<number>', help='number of pdbs for the dataset')

arg = parser.parse_args()

assert(arg.file or arg.dir)

updated = []

if arg.file:
	f = open(arg.file, 'r')
	lines = f.readlines()
	header = lines[0]
	lines = lines[1:]

	for line in lines:
		pdbid = line[:4]
		chain = line[4]
		if len(updated) <= arg.n-1:
			# extracting file info online
			with urllib.request.urlopen(f'https://files.rcsb.org/header/{pdbid}.cif') as c:
				doc = cif.read_string(c.read())
				block = doc.sole_block()
				updated.append(pdbid+chain)
		else: break
elif arg.dir:
	for fc in os.listdir(arg.dir):
		if fc.endswith(".cif"):
			updated.append(os.path.join(arg.dir, fc))
			if len(updated) >= arg.n: break

print(json.dumps(updated))
