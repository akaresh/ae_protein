#!/usr/bin/env python

import argparse
import json
import lzma
import os
import pickle
import sys

import pandas as pd
from turns_data_class import TurnFrame

parser = argparse.ArgumentParser(description='creating turn labels')
parser.add_argument('--frame', required=True, type=str,
	metavar='<path>', help='path to dataframe to split')
parser.add_argument('--train', required=True, type=str,
	metavar='<path>', help='path and base filename to save train frame')
parser.add_argument('--test', required=True, type=str,
	metavar='<path>', help='path and basr filename to save test frame')
parser.add_argument('--flank', required=True, type=int,
	metavar='<int>', help='number of flanking residues to take')
parser.add_argument('--type', required=False, type=str,
	default='ca', metavar='<str>', help='fragment type')

arg = parser.parse_args()

with lzma.open(arg.frame, 'rb') as fp:
	turn_data = pickle.load(fp)

assert(hasattr(turn_data, 'flank'))
assert(hasattr(turn_data, 'turn_type'))
assert(hasattr(turn_data, 'split'))
assert(hasattr(turn_data, 'df'))

assert(turn_data.df is not None)
assert(turn_data.flank == arg.flank)
assert(turn_data.turn_type == arg.type)
assert(turn_data.split is None)

df = turn_data.df

print(df.head(5))
print(df.columns)
print(df.shape)

types = df['turn_symbol'].unique().tolist()

print(json.dumps(types,indent=2))
print(len(types))

train_frames = list()
test_frames  = list()

for tt in types:
	subset = df[df['turn_symbol'] == tt].sample(frac=1.0).sample(frac=1.0)
	
	train_sub = subset.sample(frac=0.70).sort_index()
	
	test_sub = subset[~subset.index.isin(train_sub.index)].sort_index()
	
	train_frames.append(train_sub)
	test_frames.append(test_sub)

train_frame = pd.concat(train_frames).sample(frac=1.0).reset_index()
test_frame  = pd.concat(test_frames).sample(frac=1.0).reset_index()

print(train_frame.shape)
print(test_frame.shape)
print(df.shape)

train_data = TurnFrame(train_frame)
train_data.flank = arg.flank
train_data.turn_type = arg.type
train_data.split = 'train'

test_data = TurnFrame(test_frame)
test_data.flank = arg.flank
test_data.turn_type = arg.type
test_data.split = 'test'

save_dir = arg.train.split('/')
save_dir = '/'.join(save_dir[:-1])
assert(os.path.isdir(save_dir))

save_dir = arg.test.split('/')
save_dir = '/'.join(save_dir[:-1])
assert(os.path.isdir(save_dir))

train_save = f"{arg.train}_train_flank{arg.flank}_type{arg.type}.pickle.xz"
test_save  = f"{arg.test}_test_flank{arg.flank}_type{arg.type}.pickle.xz"
with lzma.open(train_save, 'wb') as fp:
	pickle.dump(train_data, fp)

with lzma.open(test_save, 'wb') as fp:
	pickle.dump(test_data, fp)
