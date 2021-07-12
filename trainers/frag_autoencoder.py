#!/usr/bin/python3 

import argparse
import json
import sys

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.utils.data as data_utils
import torch.optim as optim
import torchvision

class AE(nn.Module):
	def __init__(self, **kwargs):
		super().__init__()
		self.encoder_hidden_layer = nn.Linear(in_features=kwargs["input_shape"], out_features=10)
		self.encoder_output_layer = nn.Linear(in_features=10, out_features=2)
		self.decoder_hidden_layer = nn.Linear(in_features=2, out_features=10)
		self.decoder_output_layer = nn.Linear(in_features=10, out_features=kwargs["input_shape"])

	def forward(self, features):
		activation = self.encoder_hidden_layer(features)
		activation = torch.relu(activation)
		code = self.encoder_output_layer(activation)
		code = torch.relu(code)
		activation = self.decoder_hidden_layer(code)
		activation = torch.relu(activation)
		activation = self.decoder_output_layer(activation)
		reconstructed = torch.relu(activation)
		return reconstructed

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

# setting seed to just compare the results
seed = 42
# setting the random seed from pytorch random number generators
torch.manual_seed(seed)
# enabling benchmark moden in cudnn (GPU accelerated library of primitives for deep neural net)
torch.backends.cudnn.benchmark = False
# making experiments reproducible
torch.backends.cudnn.deterministic = True

#setting constants (can convert to argparse later)
batch_size = 1028
epochs = 25
learning_rate = 1e-3

coords = np.array(df.xyz_set[:-1].to_list())
print(coords.shape)

for c in coords:
	start = c[0]
	for i in range(1, c.shape[0]):
		c[i] -= start
	c[0] = c[0] - start

coords = [c.flatten() for c in coords]
coords = np.array(coords)

for i, c in enumerate(coords):
	mx = np.amax(np.abs(c))
	coords[i] = c/mx

print(coords.shape)

train = data_utils.TensorDataset(torch.Tensor(coords), torch.Tensor(coords))
train_loader = data_utils.DataLoader(train, batch_size=batch_size, shuffle=True)

#  use gpu if available
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# create a model from `AE` autoencoder class
# load it to the specified device, either gpu or cpu
model = AE(input_shape=21).to(device)

# create an optimizer object
# Adam optimizer with learning rate 1e-3
optimizer = optim.Adam(model.parameters(), lr=learning_rate)

# mean-squared error loss
criterion = nn.MSELoss()

#training autoencoder for out specified number of epochs
for epoch in range(epochs):
	loss = 0
	for batch_features, _ in train_loader:
		# reshape mini-batch data to [N, 784] matrix
		# load it to the active device
		batch_features = batch_features.view(-1, 21).to(device)

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

	# display the epoch training loss
	print("epoch : {}/{}, recon loss = {:.8f}".format(epoch + 1, epochs, loss))






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
