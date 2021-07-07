#!/usr/bin/python3 

import argparse
import json
import sys
import os

parser = argparse.ArgumentParser(description='make fragment data set')
parser.add_argument('--cifs', '-c', required=False, type=str,
	metavar='<path>', help='path to file with all of cifs we want to use')

arg = parser.parse_args()

with open(arg.cifs) as fp:
	ids = json.load(fp)
fp.close()

print(ids[0])

for cif in ids:
	print(cif)
	sys.exit()

