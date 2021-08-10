from torch import relu
from torch.nn import Dropout, Linear, Module, ModuleList, Conv2d, MaxPool2d
from torch.nn import ConvTranspose2d


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
		channels=None,
		function_list=None,
		dropout=None,
		maxpool_kernel=None,
		maxpool_stride=None,
		maxpool_padding=None,
		convd_kernel=None,
		convd_padding=None,
		convd_stride=None,
		convtd_kernel=None,
		convtd_padding=None,
		convtd_stride=None):
		
		assert(channels is not None and type(channels) == list)
		assert(function_list is not None and type(function_list) == list)
		assert(convd_kernel is not None and type(kernel) == list)
		assert(convtd_kernel is not None and type(kernel) == list)
		assert(maxpool_kernel is not None and type(maxpool_kernel) == list)
		assert(maxpool_padding is not None type(maxpool_padding) == list)
		assert(maxpool_stride is not None type(maxpool_stride) == list)
		
		super().__init__()
		self.funs = [x() for x in function_list]
		
		# checking length of input
		assert(len(channels) == len(function_list))
		
		# maxpool check
		self.maxpool = []
		for k, s, p in zip(maxpool_kernel, maxpool_stride, maxpool_padding):
			assert(k > 1 and s >= 1 and p >= 0)
			self.maxpool.append(MaxPool2d(k, stride=s, padding=p))
			
		# Conv2d
		self.conv2d = []
		prev = 1
		for c, k, s, p in zip(channels[:len(convd_kernel)+1], convd_kernel, convd_stride, convd_padding):
			assert(c > 1 and k > 1 and s >= 1 and p >=0)
			self.conv2d.append(
				Conv2d(
					in_channels=prev, out_channels=c, convd_kernel=k, convd_stride=s, 
					convd_padding=p))
			prev = c
			
		self.convt2d = []
		for c, k, s, p in zip(channels[len(self.conv2d):], convtd_kernel, convtd_stride, convtd_padding):
			assert(c > 1 and k > 1 and s >= 1 and p >=0)
			self.conv2d.append(
				ConvTranspose2d(
					in_channels=prev, out_channels=c, convd_kernel=k, convd_stride=s, 
					convd_padding=p))
			prev = c
		
		if dropout is not None:
			assert(type(dropout) == list)
			assert(len(dropout) == len(self.conv2d))
			self.dropout = []
			for dp in dropout:
				assert(dp < 1.0 and dp > 0.0)
				self.dropout.append(Dropout(dp))
		else: self.dropout = [None] * len(self.conv2d)
		
	
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
		
		#tranpose
		#relu
		
		#conv2d layers
		
		for c, f, p, d in zip(self.conv2d, self.funs, self.maxpool, self.dropout):
			x = c(x)
			x = f(x)
			x = p(x)
			if d is not None: x = d(x)
			
		#convt2d layers
		
		for c, f in zip(self.convt2d, self.funs[len(self.conv2d):]):
			x = c(x)
			x = f(x)
		
		return x
