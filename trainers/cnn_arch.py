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
	pass