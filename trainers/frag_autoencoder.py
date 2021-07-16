#!/usr/bin/python3 

import argparse
import json
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

def normalize_frag(frag):
	frag = np.array(frag)
	start = frag[0]
	for i in range(1, frag.shape[0]):
		frag[i] -= start
	frag[0] -= start
	
	justx = frag[:,0]
	justy = frag[:,1]
	justz = frag[:,2]
	
	xf = np.amax(np.abs(justx))
	yf = np.amax(np.abs(justy))
	zf = np.amax(np.abs(justz))
	
	for i in range(frag.shape[0]):
		frag[i,0] /= xf
		frag[i,1] /= yf
		frag[i,2] /= zf
	
	return frag[1:].flatten()

class AE(nn.Module):
	def __init__(self, **kwargs):
		super().__init__()
		self.encoder_hidden_layer = nn.Linear(in_features=kwargs["input_shape"], out_features=128)
		self.encoder_output_layer = nn.Linear(in_features=128, out_features=10)
		self.decoder_hidden_layer = nn.Linear(in_features=10, out_features=128)
		self.decoder_output_layer = nn.Linear(in_features=128, out_features=kwargs["input_shape"])
		
		if kwargs["dropout"]: self.dropout = nn.Dropout(kwargs["dropout"])
		
	def forward(self, features):
		activation = self.encoder_hidden_layer(features)
		activation = self.dropout(activation)
		activation = torch.relu(activation)
		code = self.encoder_output_layer(activation)
		code = self.dropout(code)
		code = torch.relu(code)
		activation = self.decoder_hidden_layer(code)
		activatio = self.dropout(activation)
		activation = torch.relu(activation)
		activation = self.decoder_output_layer(activation)
		reconstructed = torch.relu(activation)
		return reconstructed

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description=''.join(
									('initial autencoder training with CA',
									'fragments')))
	parser.add_argument('--df', '-d', required=True, type=str,
		metavar='<path>', help='path to xz compressed pandas dataframe')
	parser.add_argument('--split', '-s', required=False, type=float,
		metavar='<float>', default=0.80, help='train/test split fraction')
	
	arg = parser.parse_args()
	
	df = pd.read_pickle(arg.df, compression='xz').sample(frac=1.0)
	
	df['norm_frag'] = df.xyz_set.apply(normalize_frag)
	
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
	
# 	setting constants (can convert to argparse later)
	batch_size = 1
	epochs = 100
	learning_rate = 1e-6
	
	train_coords = np.array(df.norm_frag[:trn].to_list())
	test_coords  = np.array(df.norm_frag[trn:].to_list())
	
	train = data_utils.TensorDataset(torch.Tensor(train_coords), 
									 torch.Tensor(train_coords))
	train_loader = data_utils.DataLoader(train, 
										 batch_size=batch_size,
										 shuffle=True)
	test = data_utils.TensorDataset(torch.Tensor(test_coords),
									torch.Tensor(test_coords))
	test_loader = data_utils.DataLoader(test,batch_size=1,shuffle=True)
	
# 	use gpu if available
	device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
	
	# create a model from `AE` autoencoder class
	# load it to the specified device, either gpu or cpu
	model = AE(input_shape=18,dropout=0.25).to(device)
	s = summary(model, input_size=(16,1,18), verbose=0)
	#print(s)
	
	su = repr(s)
	print(su.encode('utf-8').decode('latin-1'))
	
	# create an optimizer object
	# Adam optimizer with learning rate 1e-3
	optimizer = optim.Adam(model.parameters(),
						   lr=learning_rate,
						   weight_decay=1e-6)
	
	# mean-squared error loss
	criterion = nn.MSELoss()
	
	#training autoencoder for out specified number of epochs
	for epoch in range(epochs):
		loss = 0
		for batch_features, _ in train_loader:
			# reshape mini-batch data to [N, 784] matrix
			# load it to the active device
			batch_features = batch_features.view(-1, 18).to(device)
		
			# reset the gradients back to zero
			# PyTorch accumulates gradients on subsequent backward passes
			optimizer.zero_grad()
			
			# compute reconstructions
			outputs = model(batch_features)
			
			# compute training reconstruction loss
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
			bv = bv.view(-1, 18).to(device)
			
			outputs = model(bv)
			test_loss = criterion(outputs, bv)
			vloss += test_loss.item()
		vloss = vloss / len(test_loader)

		# display the epoch training loss
		print(''.join((f"epoch : {epoch+1}/{epochs}, recon loss = {loss:.8f}",
					  f" test loss = {vloss:.8f}")))
	
	for name, param in model.named_parameters():
		if param.requires_grad:
			print(name, param.data)
	
	loss = 0
	for batch_features, _ in test_loader:
		batch_features = batch_features.view(-1, 18).to(device)
		
		outputs = model(batch_features)
		train_loss = criterion(outputs, batch_features)
		
		loss += train_loss.item()
	
	loss = loss / len(test_loader)
	
	print(f'testing loss: {loss}') 
	
	
	
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
