#!/usr/bin/env python

import argparse
import os
import sys

import pandas as pd
from turns_data_class import TurnFrame

parser = argparse.ArgumentParser(description=''.join((
	'pipeline for creating multiple turn fragment save frames with different',
	' fragment types and number of flanking residues')))
parser.add_argument('--ids', required=True, type=str,
	metavar='<path>', help='path to text file with pdb ids from PISCES')
parser.add_argument('--dir', required=True, type=str,
	metavar='<path>', help='directory to save chains/ and complete_structures')
parser.add_argument('--num', required=False, type=int,
	metavar='<int>', help='number of chains to take', default=None)
parser.add_argument('--bt18', required=True, type=str,
	metavar='<path>', help='path to BetaTurns18.py2')
parser.add_argument('--dssp', required=True, type=str,
	metavar='<path>', help='path to binary for dssp')
parser.add_argument('--flanks', required=True, type=int, nargs='+',
	metavar='<int>', help='number of flanking residues to take')
parser.add_argument('--types', required=True, type=str, nargs='+',
	default='ca', metavar='<str>', help='fragment types')
parser.add_argument('--full', required=True, type=str, metavar='<path>',
	help='path to save full turn fragment dataset')
parser.add_argument('--train', required=True, type=str, metavar='<path>',
	help='path to save training turn fragment dataset')
parser.add_argument('--test', required=True, type=str, metavar='<path>',
	help='path to save testing turn fragment dataset')

arg = parser.parse_args()

# start the fqa conda environment
#os.system("conda init bash")
#os.system("source ~/.bash_profile")
#conda_cmd = "source activate fqa"
#os.system(conda_cmd)

# pull pdbs
assert(os.path.isfile(arg.ids))
assert(os.path.isdir(arg.dir))

if not os.path.isdir(os.path.join(arg.dir, 'chains')):
	os.mkdir(os.path.join(arg.dir, 'chains'))
if not os.path.isdir(os.path.join(arg.dir, 'complete_structures')):
	os.mkdir(os.path.join(arg.dir, 'complete_structures'))

if arg.num:
	pull_pdbs_cmd = ''.join((
		"conda run -n fqa",
		f" python pull_pdbs.py --pdbs {arg.ids} --dir {arg.dir} ",
		f"--num {arg.num}"))
else:
	pull_pdbs_cmd = ''.join((
		"conda run -n fqa",
		f" python pull_pdbs.py --pdbs {arg.ids} --dir {arg.dir}"))

os.system(pull_pdbs_cmd)

#os.system('which mkdssp')

# make the dssp for each collected pdb
chains = os.path.join(os.path.abspath(arg.dir), 'chains')
print(f"chains: {chains}")
dssp_cmd = f"conda run -n fqa python make_dssp.py --pdbs {chains}"
os.system(dssp_cmd)

# make the bt18 turn labels
# first reset the environment
#os.system("conda deactivate")
#os.system("source activate bt18")
#os.system("which python")
#chains = os.path.join(arg.dir, 'chains')
print(chains)
bt18_make_cmd = ''.join((
	"conda run -n bt18",
	" python make_bt18_labels.py",
	f" --pdbs {chains} --bt18 {arg.bt18} --dssp {arg.dssp}"))
os.system(bt18_make_cmd)

#os.system("conda deactivate")
#os.system("source activate fqa")

# parse the bt18 output

bt18_parse_cmd = f"conda run -n fqa python parse_bt18.py --pdbs {chains}"
os.system(bt18_parse_cmd)

# make turn frames

for flank in arg.flanks:
	for frag_type in arg.types:
		print(flank)
		full_set_cmd = ''.join((
			"conda run -n fqa",
			f" python turn_frames.py --pdbs {chains} --flank {flank}",
			f" --type {frag_type} --save {arg.full}"))
		os.system(full_set_cmd)
		full_result = f"{arg.full}_full_flank{flank}_type{frag_type}.pickle.xz"
		split_set_cmd = ''.join((
			"conda run -n fqa",
			f" python train_test_split.py --frame {full_result}",
			f" --train {arg.train} --test {arg.test}",
			f" --flank {flank} --type {frag_type}"))
		os.system(split_set_cmd)