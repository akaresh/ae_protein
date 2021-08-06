from torch import relu
from torch.nn import Dropout, Linear, Module, ModuleList


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
	# kernel size and stride
	
	def __init__(
		self,
		units=None,
		function_list=None,
		dropouts=None,
		maxpool=None):
		
		assert(units is not None and type(units) == list)
		assert(function_list is not None and type(function_list) == list)
		
		assert(len(units) == len(function_list))
		
		if dropout is not None:
			assert(type(dropout) == list)
			assert(len(dropout) == len(units))
			self.dropout = []
			for dp in dropout:
				assert(dp < 1.0 and dp > 0.0)
				self.dropout.append(Dropout(dp))
		else: self.dropout = [None] * len(units)
		
		if maxpool is not None:
			assert(type(maxpool) == list)
			assert(len(maxpool) == units.count(nn.Conv2d))
		else: self.maxpool = [None] * units.count(nn.Conv2d)
		
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
		# CNN structure
		# conv
		# relu
		# pool
		# dropout
		
		out = x
		for i in range(len(self.units) - 1):
			out = self.layers[i](out)  # switches from conv2d to contranspose2d
			out = self.funs[i](out)
			if self.maxpool[i] is not None: out = self.maxpool[i](out)
			if self.dropout[i] is not None: out = self.dropout[i](out)
		reconstructed = self.funs[-1](out)
		
		return reconstructed
