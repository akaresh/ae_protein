#!/usr/bin/python3

"""
## Class definitions for AutoEncoders on protein structural fragments ##
"""
import sys

from torch import relu
from torch.nn import Dropout, Linear, Module, ModuleList
import torch.nn as nn


class SimpleAEcnn(Module):
	"""
	Class definition for a simple CNN autoencoder working on protein fragments.
	Following PyTorch AutoEncoder tutorial here:
	https://analyticsindiamag.com/how-to-implement-convolutional-autoencoder-in-pytorch-with-cuda/.
	All layers are convolutional.
	
	Parameters
	----------
	dropout: Dropout rate (optional)
		Optional dropout rate, between 0 and 1.
		float, 0 < dropout < 1
	
	Returns
	-------
	AutoEncoder model, PyTorch nn.Module object
	"""
	
	# done for CA and the fragment size of 7
	def __init__(self, dropout=None):

		if dropout is not None:
			assert(type(dropout) == float)
			assert(dropout < 1.0 and dropout > 0.0)
		
		super().__init__()
		
		# encoder
		self.conv1 = nn.Conv2d(
			in_channels=1,
			out_channels=10,
			kernel_size=2,
			padding=3)
		
		self.conv2 = nn.Conv2d(
			in_channels=10,
			out_channels=5,
			kernel_size=2,
			padding=3)
		
		self.pool = nn.MaxPool2d(2, 2)
		
		# decoder
		self.convt1 = nn.ConvTranspose2d(
			in_channels=5,
			out_channels=10,
			kernel_size=2,
			stride=1)
		
		self.convt2 = nn.ConvTranspose2d(
			in_channels=10,
			out_channels=1,
			kernel_size=2,
			stride=1)
		
		self.dropout = Dropout(dropout) if dropout is not None else None
	
	def forward(self, features):
		# conv1
		x = self.conv1(features)
		x = relu(x)
		x = self.pool(x)
		if self.dropout is not None: activate = self.dropout(x)
		
		# conv2
		x = self.conv2(x)
		x = relu(x)
		x = self.pool(x)
		if self.dropout is not None: activate = self.dropout(x)
		
		# convtranspose1
		x = self.convt1(x)
		x = relu(x)
		
		# convtranspose2
		x = self.convt2(x)
		reconstructed = relu(x)
		
		return reconstructed


class SimpleAEfc(Module):
	"""
	Class definition for a simple autoencoder working on protein fragments.
	Following PyTorch AutoEncoder tutorial here: https://gist.github.com/
	AFAgarap/4f8a8d8edf352271fa06d85ba0361f26.
	All layers are Fully Connected (FC) layers.
	
	Parameters
	----------
	inshape: Input shape of fragment data
		Model expects input fragment is a flattened array of coordinates.
		inshape > 0
	dropout: Dropout rate (optional)
		Optional dropout rate, between 0 and 1.
		float, 0 < dropout < 1
	
	Returns
	-------
	AutoEncoder model, PyTorch nn.Module object
	"""
	
	def __init__(self, inshape=None, dropout=None):
		
		assert(inshape is not None)
		assert(type(inshape) == int)
		if dropout is not None:
			assert(type(dropout) == float)
			assert(dropout < 1.0 and dropout > 0.0)
		
		super().__init__()
		
		self.encoder_hidden = Linear(in_features=inshape, out_features=128)
		self.encoder_out    = Linear(in_features=128, out_features=10)
		self.decoder_hidden = Linear(in_features=10, out_features=128)
		self.decoder_out    = Linear(in_features=128, out_features=inshape)
		
		self.dropout = Dropout(dropout) if dropout is not None else None
	
	def forward(self, features):
		"""
		Run forward pass for model.
		Forward pass goes as follows ->
			Input -> 128 FC hidden layer -> 10 FC latent code ->
				128 FC hidden layer -> Input FC layer
		Use torchinfo.summary to give summary of model
		"""
		activation = self.encoder_hidden(features)
		if self.dropout is not None: activation = self.dropout(activation)
		activation = relu(activation)
		
		latent = self.encoder_out(activation)
		if self.dropout is not None: latent = self.dropout(latent)
		latent = relu(latent)
		
		activation = self.decoder_hidden(latent)
		if self.dropout is not None: activation = self.dropout(activation)
		activation = relu(activation)
		
		reconstructed = self.decoder_out(activation)
		
		return reconstructed


class DynamicAEfc(Module):
	"""
	Class definition for AutoEncoder model with dynamic number of layers and
	units per layer. All layers are Fully Connected (FC) layers.
	
	Parameters
	----------
	inshape: Input shape of fragment data
		Model expects input fragment is a flattened array of coordinates.
		inshape > 0
	units: list for units per layer, requred
		list of units per layer, excluding size of input.
	function_list: list for non-linear functions applied at each layer, requred
		list of PyTorch non-linear functions to be applied at each layer.
		len(function_list) needs to each length of units list.
	dropouts: list of dropout probabilities per layer
		Not required.
		If specified, must be list of equal length to units and function_list.
	
	Returns
	-------
	AutoEncoder model, PyTorch nn.Module object
	"""
	
	def __init__(
		self,
		inshape=None,
		units=None,
		function_list=None,
		dropouts=None):
	
		assert(inshape is not None and type(inshape) == int)
		assert(units is not None and type(units) == list)
		assert(function_list is not None and type(function_list) == list)
		
		assert(len(units) == len(function_list))
		
		if dropouts is not None:
			assert(type(dropout) == list)
			assert(len(dropout) == len(units))
			self.dropout = []
			for dp in dropout:
				assert(dp < 1.0 and dp > 0.0)
				self.dropout.append(Dropout(dp))
		else: self.dropout = [None] * len(units)
		
		super().__init__()
		self.inshape = inshape
		self.units = units
		self.funs = function_list
		
		self.linears = []
		
		prev = self.inshape
		for uu in units:
			self.linears.append(Linear(prev, uu))
			prev = uu
		self.linears.append(Linear(prev, self.inshape))
		self.linears = ModuleList(self.linears)
	
	def forward(self, x):
		"""
		Perform forward pass in the dynamic AutoEncoder model.
		Use torchinfo.summary to find detailed summary of model.
		"""
		
		out = x
		for i in range(len(self.linears) - 1):
			out = self.linears[i](out)
			if self.dropout[i] is not None: out = self.dropout[i](out)
			out = self.funs[i](out)
		reconstructed = self.linears[-1](out)
		
		return reconstructed


if __name__ == '__main__':
	
	import argparse
	from math import floor
	import sys
	
	import numpy as np
	import pandas as pd
	from scipy.spatial.distance import cdist
	import torch
	import torch.nn as nn
	import torch.utils.data as data_utils
	import torch.optim as optim
	from torchinfo import summary
	import torchvision
	from training_tools import normalize_frag, distance_matrix, fit_model, pdb_writer
	
	parser = argparse.ArgumentParser(
		description='Test PyTorch Model definitions in library')
	parser.add_argument(
		'--df', '-d', required=True, type=str,
		metavar='<path>', help='path to xz compressed pandas dataframe')
	parser.add_argument(
		'--split', '-s', required=False, type=float,
		metavar='<float>', default=0.80, help='train/test split fraction')
	parser.add_argument(
		'--batchsize', '-b', required=False, type=int,
		metavar='<int>', default=1028, help='training batchsize')
	parser.add_argument(
		'--epochs', '-e', required=False, type=int,
		metavar='<int>', default=5, help='num of epochs to run')
	parser.add_argument(
		'--lrate', '-l', required=False, type=float,
		metavar='<float>', default=1e-3, help='learing rate')
	parser.add_argument(
		'--l2', '-w', required=False, type=float,
		metavar='<float>', default=0.0, help='l2 regularization weight')
	parser.add_argument(
		'--dropout', '-u', required=False, type=float,
		metavar='<float>', default=None, help='dropout rate')
	parser.add_argument(
		'--vis', '-v', required=False, type=bool,
		metavar='<string>', default=None, help='visualization')
	
	arg = parser.parse_args()
	
	df = pd.read_pickle(arg.df, compression='xz').sample(
		frac=1.0, random_state=42)
	
	df['norm_frag'] = df.xyz_set.apply(normalize_frag)
	df['dmatrix'] = df.xyz_set.apply(distance_matrix)
	# print(df.iloc[0])
	# sys.exit()
	
	fshape = df.norm_frag[0].shape[0]
	dshape = df.dmatrix[0].shape[0]
	
	print(df.head(3))
	print(df.tail(3))
	print(df.columns)
	print(df.shape)
	print()
	trn = floor(df.shape[0] * arg.split)
	
	# setting seed to just compare the results
	seed = 42
	# setting the random seed from pytorch random number generators
	torch.manual_seed(seed)
	# enabling benchmark mode in cudnn (GPU accelerated library of primitives
	# for deep neural net)
	torch.backends.cudnn.benchmark = False
	# making experiments reproducible
	torch.backends.cudnn.deterministic = True
	
	# use gpu if available
	device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
	
	# create a model from `AE` autoencoder class
	# load it to the specified device, either gpu or cpu
	
	model_aefc = SimpleAEfc(
		inshape=fshape, dropout=arg.dropout).to(device)
	s_aefc = summary(
		model_aefc, input_size=(arg.batchsize, 1, fshape), verbose=0)
	su_aefc = repr(s_aefc)
	print(su_aefc.encode('utf-8').decode('latin-1'))
	print()
	
	# Set the data loaders
	train_coords = np.array(df.norm_frag[:trn].to_list())
	test_coords  = np.array(df.norm_frag[trn:].to_list())


	# print(df.xyz_set[trn], len(df.xyz_set[trn]))
	# print(test_coords[0], len(test_coords[0]))

	# sys.exit()
	
	train = data_utils.TensorDataset(
		torch.Tensor(train_coords),
		torch.Tensor(train_coords))
	
	train_loader = data_utils.DataLoader(
		train,
		batch_size=arg.batchsize,
		shuffle=True)
	
	test = data_utils.TensorDataset(
		torch.Tensor(test_coords),
		torch.Tensor(test_coords))
	
	test_loader = data_utils.DataLoader(test, batch_size=1, shuffle=True)



	
	# Set loss and optimizer
	criterion = nn.L1Loss()
	optimizer = optim.Adam(
		model_aefc.parameters(),
		lr=arg.lrate,
		weight_decay=arg.l2)
	
	# Fit the model
	model_aefc = fit_model(
		model_aefc,
		train=train_loader,
		test=test_loader,
		optimizer=optimizer,
		criterion=criterion,
		device=device,
		epochs=arg.epochs)

		#if arg.vis:
	#print('Visualizer')
	# #initial
	# pdb_writer(coords=df.xyz_set[trn], seq=df.fragment_seq[trn],
	# 		   atoms=[df.fragment_type[trn][0]]*len(df.fragment_seq[trn]),
	# 		   chain=[df.chain_id[trn][0]]*len(df.fragment_seq[trn]))
	# #normalized
	# pdb_writer(coords=[df.norm_frag[trn][i:i+3] for i in range(0, len(df.norm_frag[trn]), 3)],
	# 		   seq=df.fragment_seq[trn],
	# 		   atoms=[df.fragment_type[trn][0]]*len(df.fragment_seq[trn]),
	# 		   chain=[df.chain_id[trn][0]]*len(df.fragment_seq[trn]))
	#normalized after training
	saved = model_aefc.forward(test[0][0])
	pdb_writer(coords=[saved[i:i+3] for i in range(0, len(saved),3)],
			   seq=df.fragment_seq[trn],
			   atoms=[df.fragment_type[trn][0]]*len(df.fragment_seq[trn]),
			   chain=[df.chain_id[trn][0]]*len(df.fragment_seq[trn]))
	#pic 1. how the frag look like, pic 2 (frag normalized before and after the model)

	
	model_dyn_aefc = DynamicAEfc(
		inshape=fshape,
		dropouts=arg.dropout,
		units=[128, 10, 128],
		function_list=[relu, relu, relu]).to(device)
	
	s_dyn_aefc = summary(
		model_dyn_aefc,
		input_size=(arg.batchsize, 1, fshape),
		verbose=0)
	su_dyn_aefc = repr(s_dyn_aefc)
	print(su_dyn_aefc.encode('utf-8').decode('latin-1'))
	print()
	
	# Set optimizer
	optimizer = optim.Adam(
		model_dyn_aefc.parameters(),
		lr=arg.lrate,
		weight_decay=arg.l2)
	
	# Fit the model
	model_dyn_aefc = fit_model(
		model_dyn_aefc,
		train=train_loader,
		test=test_loader,
		optimizer=optimizer,
		criterion=criterion,
		epochs=arg.epochs,
		device=device)
	
	model_aecnn = SimpleAEcnn(dropout=arg.dropout).to(device)
	s_aecnn = summary(
		model_aecnn,
		input_size=(arg.batchsize, 1, dshape, dshape),
		verbose=0)
	
	su_aecnn = repr(s_aecnn)
	print(su_aecnn.encode('utf-8').decode('latin-1'))
	print()
	
	# Set optimizer
	optimizer = optim.Adam(
		model_aecnn.parameters(),
		lr=arg.lrate,
		weight_decay=arg.l2)
	
	# Set the data loaders
	train_coords = np.array(df.dmatrix[:trn].to_list())
	test_coords  = np.array(df.dmatrix[trn:].to_list())
	
	train = data_utils.TensorDataset(
		torch.Tensor(train_coords),
		torch.Tensor(train_coords))
	
	train_loader = data_utils.DataLoader(
		train,
		batch_size=arg.batchsize,
		shuffle=True)
	
	test = data_utils.TensorDataset(
		torch.Tensor(test_coords),
		torch.Tensor(test_coords))
	
	test_loader = data_utils.DataLoader(test, batch_size=1, shuffle=True)
	
	# Fit the model
	model_aecnn = fit_model(
		model_aecnn,
		train=train_loader,
		test=test_loader,
		optimizer=optimizer,
		criterion=criterion,
		epochs=arg.epochs,
		device=device)

	#visualization here
