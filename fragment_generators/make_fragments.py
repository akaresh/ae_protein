#!/usr/bin/python3 

"""
Script to generate dataset of protein structure fragments. 
"""

import argparse
import json
import pandas as pd
import sys

from fragment_lib import make_atom_frame, make_fragment_frame
import numpy as np

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)


parser = argparse.ArgumentParser(description='making fragments')
parser.add_argument(
	'--file', required=False, type=str, metavar='<path>', 
	help='path to cif file')
parser.add_argument(
	'--cifs', '-c', required=False, type=str, metavar='<path>', 
	help='path to file containing set of cifs to build')
parser.add_argument(
	'--atoms', '-a', required=False, type=str, default=None, metavar='<path>',
	help='path to prebuilt atom data frame')
parser.add_argument(
	'--size', required=False, type=int, default=3, metavar='<int>',
	help='length of fragments')
parser.add_argument(
	'--type', required=False, type=str, default='CA', metavar='<int>',
	help='type of fragments: CA, bb, bbcen, bbsc')
parser.add_argument(
	'--outatoms', '-u', required=False, type=str, default=None,
	metavar='<path>', help='save location for atom frame')
parser.add_argument(
	'--outfrag', '-o', required=False, type=str, default=None,
	metavar='<path>', help='save location for fragment frame')

arg = parser.parse_args()

assert(arg.type == 'CA' or arg.type == 'bb' or arg.type == 'bbcen'
	   or arg.type == 'bbsc')

if arg.file:
	if arg.atoms: atm_frame = pd.read_pickle(arg.atoms, compression='xz')
	else:
		atm_frame = make_atom_frame([arg.file])
		print(atm_frame.head(1000000))
		#sys.exit()
	
	frag_df = make_fragment_frame(atm_frame, arg.size, ftype=arg.type)
	
	print(frag_df.head(100000000))
	print(frag_df.columns)
	print(frag_df.shape)
	print(frag_df.xyz_set[0])
	print(np.array(frag_df.xyz_set[0]).shape)
	
	if arg.outfrag: frag_df.to_pickle(arg.outfrag, compression='xz')
	if arg.atoms == None:
		if arg.outatoms: atm_frame.to_pickle(arg.outatoms, compression='xz')
	
elif arg.cifs:
	with open(arg.cifs, 'r') as fp:
		cifs = json.load(fp)
	fp.close()
	
	if arg.atoms: atm_frame = pd.read_pickle(arg.atoms, compression='xz')
	else:
		with open(arg.cifs, 'r') as fp:
			cifs = json.load(fp)
		fp.close()
		atm_frame = make_atom_frame(cifs)
		print(atm_frame.head(10000000))
	
	print('finished atom frame')
	
	frag_df = make_fragment_frame(atm_frame, arg.size, ftype=arg.type)
	
	print(frag_df.head(100000))
	print(frag_df.columns)
	print(frag_df.shape)
	
	if arg.outfrag: frag_df.to_pickle(arg.outfrag, compression='xz')
	if arg.atoms == None:
		if arg.outatoms: atm_frame.to_pickle(arg.outatoms, compression='xz')
