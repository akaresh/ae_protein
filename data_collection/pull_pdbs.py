#!/usr/bin/python3 

import argparse
import csv
import os
import sys
import urllib.request

import Bio.PDB as bpdb
from Bio.PDB import MMCIF2Dict

class ChainSelector(bpdb.Select):
	def __init__(self, chainid):
		self.chainid = chainid
	def accept_model(self, model):
		if model.id == 0: return 1
		else:          return 0
	def accept_chain(self, chain):
		if chain.id == self.chainid: return 1
		else:                return 0

parser = argparse.ArgumentParser(description='creating chains/')
parser.add_argument('--pdbs', required=True, type=str,
	metavar='<path>', help='path to text file with pdb ids from PISCES')
parser.add_argument('--dir', required=True, type=str,
	metavar='<path>', help='directory to save chains/ and complete_structures')

arg = parser.parse_args()

path = os.path.join(arg.dir, 'complete_structures')
out  = os.path.join(arg.dir, 'chains')
with open(arg.pdbs, 'r') as fp:
	entries = fp.readlines()
fp.close()

entries = entries[1:]
for i, row in enumerate(entries):
	row = row.rstrip()
	r = row.split()
	assert(len(r) == 6) # make sure it is from PISCES
	pdbid = r[0][:4].lower()
	chainid = r[0][4:]
	print(pdbid, chainid)
	#print('downloading')
	
	if not os.path.isfile(f'{path}/{pdbid}.cif'):
		urllib.request.urlretrieve(
			f'http://files.rcsb.org/download/{pdbid}.cif',
			f'{path}/{pdbid}.cif')
	if not os.path.isfile(f'{path}/{pdbid}.pdb'):
		try:
			urllib.request.urlretrieve(
				f'http://files.rcsb.org/download/{pdbid}.pdb',
				f'{path}/{pdbid}.pdb')
		except:
			continue
	
# 	if not os.path.isfile(f'{out}/{pdbid}{chainid.lower()}.cif'):
# 		parser = bpdb.MMCIFParser()
# 		filep = path+'/'+pdbid+'.cif'
# 		try:
# 			structure = parser.get_structure(pdbid+chainid,filep)
# 		except:
# 			continue
# 		cifdic = MMCIF2Dict.MMCIF2Dict(filep)
# 		
# 		io = bpdb.MMCIFIO()
# 		io.set_dict(cifdic)
# 		io.set_structure(structure)
# 		
# 		outf = out+'/'+pdbid+chainid.lower()+'.cif'
# 		with open(outf, 'a') as fp:
# 			io._save_dict(fp)
# 			io._save_structure(
# 				fp,
# 				select=ChainSelector(chainid),
# 				preserve_atom_numbering=False)
# 		fp.close()
	
	if not os.path.isfile(f'{out}/{pdbid}{chainid.lower()}.pdb'):
		parser = bpdb.PDBParser()
		filep = path+'/'+pdbid+'.pdb'
		try:
			structure = parser.get_structure(pdbid+chainid,filep)
		except:
			continue
		
		io = bpdb.PDBIO()
		io.set_structure(structure)
		outf = out+'/'+pdbid+chainid.lower()+'.pdb'
		with open(outf, 'w') as fp:
			io.save(
				fp,
				select=ChainSelector(chainid),
				preserve_atom_numbering=False)
		fp.close()