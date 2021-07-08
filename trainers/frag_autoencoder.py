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


