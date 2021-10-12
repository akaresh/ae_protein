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

path = '/share/korflab/data/pdb25/complete_structures/'
out  = '/share/korflab/data/pdb25/chains/'

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
	print('downloading')
	urllib.request.urlretrieve(f'http://files.rcsb.org/download/{pdbid}.cif', 
								   f'{path}{pdbid}.cif')
	
	if os.path.isfile(f'{out}{pdbid}{chainid.lower()}.cif'):
		#print('skip')
		continue
	
	#try:
	parser = bpdb.MMCIFParser()
	filep = path+pdbid+'.cif'
	print(row)
	print(filep)
	structure = parser.get_structure(pdbid+chainid,filep)
	cifdic = MMCIF2Dict.MMCIF2Dict(filep)
	#print(cifdic)
	#print(list(cifdic.keys()))
	stop=False
	for k in cifdic.keys():
		if 'unobs' in k: 
			print(k)
			stop=True
			break
	io = bpdb.MMCIFIO()
	io.set_dict(cifdic)
	if stop: 
		print('we set it right???')
		print(list(io.dic.keys()))
		print()
	io.set_structure(structure)
	if stop: 
		print('dic after set_structure')
		print(list(io.dic.keys()))
		print()
	outf = out+pdbid+chainid.lower()+'.cif' 
	with open(outf, 'a') as fp:
		#io._save_structure(fp, select=ChainSelector(chainid), preserve_atom_numbering=True)
		io._save_dict(fp)
		io._save_structure(
			fp,
			select=ChainSelector(chainid),
			preserve_atom_numbering=True)
		if stop:
			print('io obj after saving')
			print(dir(io))
			print(list(io.dic.keys()))
			print(list(cifdic.keys()))
			sys.exit()
	fp.close()
#	except KeyboardInterrupt:
#		sys.exit()
#	except:
#		print('could not resolve parser')
#		continue
