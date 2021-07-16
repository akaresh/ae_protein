#!/usr/bin/python3 

"""
## Class definitions for AutoEncoders on protein structural fragments ##
"""

from torch import relu
from torch.nn import Linear, Dropout

class SimpleAEff(nn.Module):
	"""
	Class definition for a simple autoencoder working on protein fragments.
	Following PyTorch AutoEncoder tutorial here: https://gist.github.com/
	AFAgarap/4f8a8d8edf352271fa06d85ba0361f26
	
	Parameters
	----------
	inshape: Input shape of fragment
		Model expects input fragment is a flattened array of coordinates. 
		inshape > 0
	dropout: Dropout rate (optional)
		Optional dropout rate, between 0 and 1. 
		float, 0 < dropout < 1
	
	Returns
	-------
	AutoEncoder Model PyTorch nn.Module object
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
		
		return reconstruted

class DynamicAEff(nn.Module):
	pass
