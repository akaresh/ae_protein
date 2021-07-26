from torch import relu
from torch.nn import Dropout, Linear, Module, ModuleList

class SimpleAEcnn(Module):
	#done for CA and the fragment size of 7
	def __init__ (self, 
				 dropout = None):
		if droupout != None: 
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
		
		#self.dropout = Dropout(dropout) if dropout != None else None
	
	def forward(self, features):
		#conv1
		x = self.conv1(features)
		x = relu(x)
		x = self.pool(x)
		#if self.dropout != None: activate = self.dropout(x)
		
		#conv2 
		x = self.conv2(x)
		x = relu(x)
		x = self.pool(x)
		#if self.dropout != None: activate = self.dropout(x)
		
		#convtranspose1
		x = self.convt1(x)
		x = relu(x)
		
		#convtranspose2
		x = self.convt2(x)
		reconstructed = relu(x)		
		
		return reconstructed

class DynamicAEcnn(Module):
	"""
	Class definition for AutoEncoder model with dynamic number of layers and
	units per layer. All layers are Convolutional. 
	
	Parameters
	----------
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
				units=None,
				function_list=None,
				dropouts=None,
				maxpool=None):
		
		assert(units != None and type(units) == list)
		assert(function_list != None and type(function_list) == list)
		
		assert(len(units) == len(function_list))
		
		if dropout != None:
			assert(type(dropout) == list)
			assert(len(dropout) == len(units))
			self.dropout = []
			for dp in dropout:
				assert(dp < 1.0 and dp > 0.0)
				self.dropout.append(Dropout(dp))
		else: self.dropout = [None]*len(units)
		
		if maxpool != None:
			assert(type(maxpool) == list)
			assert(len(maxpool) == units.count(nn.Conv2d))
		else: self.maxpool = [None]*units.count(nn.Conv2d)
		
		super(DynamicAEcnn, self).__init__()
		self.units = units
		self.funs = [x() for x in function_list]
		
		self.layers = []
		
		prev = 1
		for uu in units:
			self.layers.append(Linear(prev, uu))
			prev = uu
		self.layers.append(prev, 1)
		self.layers = ModuleList(self.layers)
	
	def forward(self, x):
		"""
		Perform forward pass in the dynamic AutoEncoder model.
		Use torchinfo.summary to find detailed summary of model. 
		"""
		#CNN structure
		#conv
		#relu
		#pool
		#dropout
		
		out = x
		for i in range(len(self.units)-1):
			out = self.layers[i](out)
			out = self.funs[i](out)
			if self.maxpool[i] != None: out = self.maxpool[i](out)
			if self.dropout[i] != None: out = self.dropout[i](out)
		reconstructed = self.funs[-1](out)
		
		return reconstructed