#!/usr/bin/python3 

import argparse
import json
import pandas as pd
import sys

parser = argparse.ArgumentParser(description=''.join(
								('initial autencoder training with CA',
								 'fragments')))
parser.add_argument('--df', '-d', required=True, type=str,
	metavar='<path>', help='path to xz compressed pandas dataframe')

arg = parser.parse_args()

df = pd.read_pickle(arg.df, compression='xz')

print(df.head(3))
print(df.columns)
print(df.shape)

"""
- train.py --model --batchsize --epochs --splits --df
- models/ FragmentAutoencoder.py
	- models are dynamically constructed
	- fully connected multi-layer autoencoders
		- flatten the coordinate matrix
		- (7x3) initial size
		- convolutions on the 7x3 or
		- convolutions on the distance matrix
		- 7x7 matrix 
		- normalize distance values between 0-1

1. FC with flat coordinates
2. CNN with matrix coordinates
3. CNN with distance matrix

conceptual pt
learn lower dimensional representations of local structure
21->10->2->10->21
21->200->10->2->10->200->21

how are we going to test

"""
