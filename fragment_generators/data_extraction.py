#!/usr/bin/python3 

import argparse
import sys

from data_collection_library import make_CAframe, make_bbframe, make_bbcen, make_bbsc

parser = argparse.ArgumentParser(description='making fragments')
parser.add_argument('--file', required=False, type=str,
	metavar='<path>', help='path to cif files')
parser.add_argument('--cifs', '-c', required=False, type=str,
	metavar='<path>', help='path to file containing set of cifs to build')
parser.add_argument('--size', required=False, type=int, default = 3,
	metavar='<int>', help='length of fragments')
parser.add_argument('--t', required=False, type=str, default='CA',
	metavar='<int>', help='type of fragments: CA, bbframe, bbcen, bbsc')

arg = parser.parse_args()

# def make_fragment_frame(file, size, t):
# 	if t == 'CA':
# 		return make_CAframe(file, size)
# 	elif t == 'bb':
# 		return make_bbframe(file, size)
# 	elif t == 'bb+cen':
# 		return make_bbcen(file, size) # separate function to calculate a centroid
# 	elif t == 'bb+sc':
# 		return make_bbsc(file, size) # function to calculate cb-extended for glycine
# 	else:
# 		raise Exception(f'un-suppported fragment tyoe: {t}')
# 		return None
# 
# 
# make_fragment_frame(arg.file, arg.size, arg.t)
