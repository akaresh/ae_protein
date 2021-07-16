#!/usr/bin/python3 

"""
## Class definitions for AutoEncoders on protein structural fragments ##
"""

from torch import relu
from torch.nn import Dropout, Linear, Module, ModuleList

class SimpleAEff(Module):
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
		
		return reconstruted

class DynamicAEff(Module):
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
	pass