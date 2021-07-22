#!/usr/bin/python3 

"""
## Class definitions for AutoEncoders on protein structural fragments ##
"""

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
	
	#done for CA and the fragment size of 7
	def __init__ (self, 
				  dropout = None):

		if dropout != None:
			assert(type(dropout) == float)
			assert(dropout < 1.0 and dropout > 0.0)
		
		super().__init__()
		
		#structure
		
		#encoder
		self.conv1 = nn.Conv2d(in_channels = 1, out_channels = 10, kernel_size = 2, 
							   padding = 3)
		self.conv2 = nn.Conv2d(in_channels = 10, out_channels = 5, kernel_size = 2, 
		                       padding = 3)
		
		self.pool = nn.MaxPool2d(2, 2)
		
		#decoder
		self.convt1 = nn.ConvTranspose2d(in_channels = 5, out_channels = 10, kernel_size = 2,
										 stride = 1)
		self.convt2 = nn.ConvTranspose2d(in_channels = 10, out_channels = 1, kernel_size = 2,
										 stride = 1)
		
		self.dropout = Dropout(dropout) if dropout != None else None
	
	def forward(self, features):
		#conv1
		x = self.conv1(features)
		x = relu(x)
		x = self.pool(x)
		if self.dropout != None: activate = self.dropout(x)
		
		#conv2 
		x = self.conv2(x)
		x = relu(x)
		x = self.pool(x)
		if self.dropout != None: activate = self.dropout(x)
		
		#convtranspose1
		x = self.convt1(x)
		x = relu(x)
		
		#convtranspose2
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
	
	def __init__(self,
				 inshape=None,
				 dropout=None):
		
		assert(inshape != None)
		assert(type(inshape) == int)
		if dropout != None:
			assert(type(dropout) == float)
			assert(dropout < 1.0 and dropout > 0.0)
		
		super().__init__()
		
		self.encoder_hidden = Linear(in_features=inshape,out_features=128)
		self.encoder_out    = Linear(in_features=128, out_features=10)
		self.decoder_hidden = Linear(in_features=10, out_features=128)
		self.decoder_out    = Linear(in_features=128, out_features=inshape)
		
		self.dropout = Dropout(dropout) if dropout != None else None
	
	def forward(self, features):
		"""
		Run forward pass for model. 
		Forward pass goes as follows ->
			Input -> 128 FC hidden layer -> 10 FC latent code -> 
				128 FC hidden layer -> Input FC layer
		Use torchinfo.summary to give summary of model
		"""
		activation = self.encoder_hidden(features)
		if self.dropout != None: activation = self.dropout(activation)
		activation = relu(activation)
		
		latent = self.encoder_out(activation)
		if self.dropout != None: latent = self.dropout(latent)
		latent = relu(latent)
		
		activation = self.decoder_hidden(latent)
		if self.dropout != None: activation = self.dropout(activation)
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
	
	def __init__(self,
				inshape=None,
				units=None,
				function_list=None,
				dropouts=None):
		
		assert(inshape != None and type(inshape) == int)
		assert(units != None and type(units) == list)
		assert(function_list != None and type(functio_list) == list)
		
		assert(len(units) == len(function_list))
		
		if dropout != None:
			assert(type(dropout) == list)
			assert(len(dropout) == len(units))
			self.dropout = []
			for dp in dropout:
				assert(dp < 1.0 and dp > 0.0)
				self.dropout.append(Dropout(dp))
		else: self.dropout = [None]*len(units)
		
		super(DynamicAEff, self).__init__()
		self.inshape = inshape
		self.units = units
		self.funs = [x() for x in function_list]
		
		self.linears = []
		
		prev = self.inshape
		for uu in units:
			self.linears.append(Linear(prev, uu))
			prev = uu
		self.linears.append(prev, self.inshape)
		self.linears = ModuleList(self.linears)
	
	def forward(self, x):
		"""
		Perform forward pass in the dynamic AutoEncoder model.
		Use torchinfo.summary to find detailed summary of model. 
		"""
		
		out = x
		for i in range(len(self.units)-1):
			out = self.linears[i](out)
			if self.dropout[i] != None: out = self.dropout[i](out)
			out = self.funs[i](out)
		reconstructed = self.funs[-1](out)
		
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
	from training_tools import normalize_frag, distance_matrix
	
	parser = argparse.ArgumentParser(description=''.join(
									('Test training for PyTorch Model Class',
									'definitions in this library')))
	parser.add_argument('--df', '-d', required=True, type=str,
		metavar='<path>', help='path to xz compressed pandas dataframe')
	parser.add_argument('--split', '-s', required=False, type=float,
		metavar='<float>', default=0.80, help='train/test split fraction')
	parser.add_argument('--batchsize', '-b', required=False, type=int,
		metavar='<int>', default=1028, help='training batchsize')
	parser.add_argument('--epochs', '-e', required=False, type=int,
		metavar='<int>', default=20, help='num of epochs to run')
	parser.add_argument('--lrate', '-l', required=False, type=float,
		metavar='<float>', default=1e-3, help='learing rate')
	parser.add_argument('--l2', '-w', required=False, type=float,
		metavar='<float>', default=0.0, help='l2 regularization weight')
	parser.add_argument('--dropout', '-u', required=False, type=float,
		metavar='<float>', default=None, help='dropout rate')
	
	arg = parser.parse_args()
	
	df = pd.read_pickle(arg.df, compression='xz').sample(frac=1.0,
														 random_state=42)
	
	df['norm_frag'] = df.xyz_set.apply(normalize_frag)
	df['dmatrix'] = df.xyz_set.apply(distance_matrix)
	
	fshape = df.norm_frag[0].shape[0]
	
	print(df.head(3))
	print(df.columns)
	print(df.shape)
	
	trn = floor(df.shape[0]*arg.split)
	
	# 	setting seed to just compare the results
	seed = 42
	# 	setting the random seed from pytorch random number generators
	torch.manual_seed(seed)
	# 	enabling benchmark mode in cudnn (GPU accelerated library of primitives 
	# 	for deep neural net)
	torch.backends.cudnn.benchmark = False
	# 	making experiments reproducible
	torch.backends.cudnn.deterministic = True
	
	train_coords = np.array(df.dmatrix[:trn].to_list())
	test_coords  = np.array(df.dmatrix[trn:].to_list())
	
	train = data_utils.TensorDataset(torch.Tensor(train_coords), 
									 torch.Tensor(train_coords))
	train_loader = data_utils.DataLoader(train, 
										 batch_size=arg.batchsize,
										 shuffle=True)
	test = data_utils.TensorDataset(torch.Tensor(test_coords),
									torch.Tensor(test_coords))
	test_loader = data_utils.DataLoader(test,batch_size=1,shuffle=True)
	
	# 	use gpu if available
	device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
	
	# create a model from `AE` autoencoder class
	# load it to the specified device, either gpu or cpu
	model = SimpleAEcnn().to(device)
	#s = summary(model, input_size=(arg.batchsize,1,fshape), verbose=0)
	#print(s)
	
	#su = repr(s)
	#print(su.encode('utf-8').decode('latin-1'))
	
	# create an optimizer object
	# Adam optimizer with learning rate 1e-3
	optimizer = optim.Adam(model.parameters(),
						   lr=arg.lrate,
						   weight_decay=arg.l2)
	
	# mean-squared error loss
	criterion = nn.L1Loss()
	
	#training autoencoder for out specified number of epochs
	for epoch in range(arg.epochs):
		loss = 0
		for batch_features, _ in train_loader:
			# reshape mini-batch data to [N, 784] matrix
			# load it to the active device
			#batch_features = batch_features.view(-1, fshape).to(device)
			batch_features = batch_features.to(device)
		
			# reset the gradients back to zero
			# PyTorch accumulates gradients on subsequent backward passes
			optimizer.zero_grad()
			
			# compute reconstructions
			outputs = model(batch_features)
			
			# compute training reconstruction loss
			#print(batch_features.size())
			#print(outputs.size())
			train_loss = criterion(outputs, batch_features)

			
			# compute accumulated gradients
			train_loss.backward()
			
			# perform parameter update based on current gradients
			optimizer.step()
			
			# add the mini-batch training loss to epoch loss
			loss += train_loss.item()
			
		# compute the epoch training loss
		loss = loss / len(train_loader)
		
		vloss = 0
		for bv, _ in test_loader:
			#bv = bv.view(-1, fshape).to(device)
			bv = bv.to(device)
			
			outputs = model(bv)
			test_loss = criterion(outputs, bv)
			vloss += test_loss.item()
		vloss = vloss / len(test_loader)

		# display the epoch training loss
		print(''.join((f"epoch : {epoch+1}/{arg.epochs}, recon loss = {loss:.8f}",
					  f" test loss = {vloss:.8f}")))
	
	loss = 0
	for batch_features, _ in test_loader:
		#batch_features = batch_features.view(-1, fshape).to(device)
		batch_features = batch_features.to(device)
		
		outputs = model(batch_features)
		train_loss = criterion(outputs, batch_features)
		
		loss += train_loss.item()
	
	loss = loss / len(test_loader)
	
	print(f'testing loss: {loss}')