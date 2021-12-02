#!/usr/bin/python3 

import argparse
import csv
import os
import sys
import urllib.request

import Bio.PDB as bpdb
from Bio.PDB import MMCIF2Dict

"""
Download structures culled from PISCES server.
Arguments
---------
	* pdbs: pdb_id list from PISCES server (required)
	* dir: directory to store the complete_structures/ and chains/ folders. 
	complete_structures/ is for saving the full pdb the culled PISCES chain
	comes from. chains/ just saves the chains. (required)
	* num: number of chains to download if a smaller dataset is desired. No
	shuffling. (optional)
Outputs
-------
	* downloaded protein chains culled from PISCES server into folder of choice.
"""

# needed to select chain from BioPython Structure object for saving chains
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
parser.add_argument('--num', '-n', required=False, type=int,
	metavar='<int>', help='number of chains to take', default=None)

arg = parser.parse_args()

full_pdbs   = os.path.join(arg.dir, 'complete_structures')
chain_pdbs  = os.path.join(arg.dir, 'chains')
with open(arg.pdbs, 'r') as fp:
	entries = fp.readlines()
fp.close()

entries = entries[1:]
for i, row in enumerate(entries):
	if arg.num is not None:
		if i > arg.num: break
	
	row = row.rstrip()
	r = row.split()
	assert(len(r) == 6) # make sure it is from PISCES
	pdbid = r[0][:4].lower()
	chainid = r[0][4:]
	print(pdbid, chainid)
	
	if not os.path.isfile(f'{full_pdbs}/{pdbid}.pdb'):
		try:
			urllib.request.urlretrieve(
				f'http://files.rcsb.org/download/{pdbid}.pdb',
				f'{full_pdbs}/{pdbid}.pdb')
		except:
			continue
	
	if not os.path.isfile(f'{chain_pdbs}/{pdbid}{chainid.lower()}.pdb'):
		parser = bpdb.PDBParser()
		filep = full_pdbs+'/'+pdbid+'.pdb'
		try:
			structure = parser.get_structure(pdbid+chainid,filep)
		except:
			continue
		
		io = bpdb.PDBIO()
		io.set_structure(structure)
		outf = chain_pdbs+'/'+pdbid+chainid.lower()+'.pdb'
		with open(outf, 'w') as fp:
			io.save(
				fp,
				select=ChainSelector(chainid),
				preserve_atom_numbering=False)
		fp.close()