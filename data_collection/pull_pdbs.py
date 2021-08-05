#!/usr/bin/python3 

"""
this takes in the pdbid list from pisces
this needs documenting
"""


import csv
import os
import sys
import urllib.request

import Bio.PDB as bpdb

class ChainSelector(bpdb.Select):
	def __init__(self, chainid):
		self.chainid = chainid
	def accept_model(self, model):
		if model.id == 0: return 1
		else:          return 0
	def accept_chain(self, chain):
		if chain.id == self.chainid: return 1
		else:                return 0

path = '/share/korflab/data/pdb40/complete_structures/'
out  = '/share/korflab/data/pdb40/chains/'

with open(sys.argv[1], 'r') as fp:
	entries = fp.readlines()
fp.close()

entries = entries[1:]
pdbl = bpdb.PDBList(pdb=path,obsolete_pdb=False)
for i, row in enumerate(entries):
	row = row.rstrip()
	r = row.split()
	pdbid = r[0][:4].lower()
	chainid = r[0][4:]
	
	if not os.path.isfile(f'{path}{pdbid}.cif'):
		print('downloading')
		urllib.request.urlretrieve(f'http://files.rcsb.org/download/{pdbid}.cif', 
								   f'{path}{pdbid}.cif')
	
	if os.path.isfile(f'{out}{pdbid}{chainid.lower()}.cif'):
		#print('skip')
		continue
	
	try:
		parser = bpdb.MMCIFParser()
		filep = path+pdbid+'.cif'
		print(row)
		print(filep)
		structure = parser.get_structure(pdbid+chainid,filep)
	
		io = bpdb.MMCIFIO()
		io.set_structure(structure)
		outf = out+pdbid+chainid.lower()+'.cif' 
		io.save(outf, select=ChainSelector(chainid))
	except KeyboardInterrupt:
		sys.exit()
	except:
		print('could not resolve parser')
		continue