#!/usr/bin/python3

"""
Shared library for functions for PyTorch NN training
"""

import sys

import datetime
import numpy as np
from scipy.spatial.distance import cdist
import torch
import Bio.Data.IUPACData as conv

# Dictionary for default values for convolutions in encoder networks
encoder_template = {
	'channels'      : 1,
	'conv_ks'       : (1, 1),
	'pool_ks'       : (1, 1),
	'conv_paddings' : (0, 0),
	'pool_paddings' : (0, 0),
	'conv_strides'  : (1, 1),
	'pool_strides'  : (1, 1)
}


# Dictionary for default values for convolution-transposes in decoder networks
decoder_template = {
	'channels'       : 1,
	'convt_ks'       : (1, 1),
	'convt_strides'  : (1, 1),
	'convt_paddings' : (0, 0)
}


class Dict2Obj(object):
	"""
	Simple class to turn a dictionary to an object.
	Useful for named attributes in convolution/convolution-transpose size
	calculations.
	
	Inputs
	------
	dic:      dictionary to object-ize.
	template: dictionary for default values.
		template dictionary is different dependening on convolution or
		convolution transpose. Only keys in template not in dic are added as
		object attributes.
	
	Returns
	------
	obj: obj with named attributes from input dictionary keys.
	"""
	def __init__(self, dic, template):
		for key in dic:
			setattr(self, key, dic[key])
		for key in template:
			if hasattr(self, key): continue
			setattr(self, key, template[key])


def normalize_frag(frag):
	"""
	Normalize structural fragment coordinates
	First atom of fragment is set to the origin.
	Each x,y,z axis normalized by the largest value in each respective
	dimension.
	
	Input
	-----
	frag: list of coordinates
	
	Returns
	-------
	frag: returns the normalized fragment.
		The first coordinate is removed since it is 0 and will not have effect
		the FC networks we want to build.
		The numpy array that is returned is flattened.
	"""
	assert(type(frag) == list)
	
	frag = np.array(frag)
	start = frag[0]
	for i in range(1, frag.shape[0]):
		frag[i] -= start
	frag[0] -= start
	return frag[1:].flatten()
	
	justx = frag[:, 0]
	justy = frag[:, 1]
	justz = frag[:, 2]
	
	xf = np.amax(np.abs(justx))
	yf = np.amax(np.abs(justy))
	zf = np.amax(np.abs(justz))
	
	for i in range(frag.shape[0]):
		frag[i, 0] /= xf
		frag[i, 1] /= yf
		frag[i, 2] /= zf
	
	return frag[1:].flatten()


def distance_matrix(frag):
	"""
	Create the matrix of distances between all pairs of atoms in fragment.
	Matrix that is built is symmetric.
	
	Example: input fragment is 8 coordinates, output matrix is 8x8
	Using cdist from scipy.spatial.distance and metric euclidean
	
	Input
	-----
	frag: list of coordinates for atoms in fragment
	
	Returns
	------
	mat: numpy symmetric matrix of distances.
		Matrix is reshaped to tensor for input into Torch CNN.
		Dimensions are inflated. Example: matrix 8x8, result is (1, 8, 8)
	"""
	assert(type(frag) == list)
	
	mat = cdist(frag, frag, metric='euclidean')
	mat = np.reshape(mat, (1, mat.shape[0], mat.shape[1]))
	return mat


def fit_model(
	model,
	train=None,
	test=None,
	optimizer=None,
	criterion=None,
	device=None,
	epochs=None):
	"""
	Fit function for all Torch Models
	
	Parameters
	----------
	train:     PyTorch dataloader
	test:      PyTorch dataloader
	device:    CPU/GPU device
	optimizer: PyTorch optimizer
	criterion: PyTorch loss
	epochs:    Number of epochs
	
	Returns
	-------
	Returns fitted model
	"""
	assert(model is not None)
	assert(hasattr(model, 'forward'))
	assert(train is not None)
	assert(hasattr(train, 'batch_size'))
	assert(type(epochs) == int and epochs > 0)
	assert(device is not None)
	assert(optimizer is not None)
	assert(hasattr(optimizer, 'zero_grad'))
	assert(criterion is not None)
	assert(hasattr(criterion, 'forward'))
	
	if test is not None: assert(hasattr(test, 'batch_size'))
	
	for epoch in range(epochs):
		loss = 0
		for batch_features, _ in train:
			batch_features = batch_features.to(device)
			
			# reset the gradients back to zero
			# PyTorch accumulates gradients on subsequent backward passes
			optimizer.zero_grad()
			
			# compute reconstructions
			outputs = model.forward(batch_features)
			
			# compute training reconstruction loss
			train_loss = criterion(outputs, batch_features)
			
			# compute accumulated gradients
			train_loss.backward()
				
			# perform parameter update based on current gradients
			optimizer.step()
			
			# add the mini-batch training loss to epoch loss
			loss += train_loss.item()
			
		# compute the epoch training loss
		loss = loss / len(train)
		
		vloss = 0
		for bv, _ in test:
			bv = bv.to(device)
			
			outputs = model.forward(bv)
			test_loss = torch.sqrt(criterion(outputs, bv))
			vloss += test_loss.item()
		vloss = vloss / len(test)
		
		# display the epoch training loss
		print(f"epoch : {epoch+1}/{epochs} recon loss = {loss:.8f} \
		test loss = {vloss:.8f}")
	return model


def pdb_writer(atoms=None, seq=None, chain=None, coords=None):
	"""
	Write fragment coordinates in pdb format
	
	Parameters
	----------
	atoms: List of atom-types. Ex: ['CA'] * 7
	seq: Amino acid identity foreach residue in the fragment.
		Needs to be same length as number of atoms in fragment.
	chain: Chain ids for each atom in the fragment. Can use dummy strings.
		Needs to be same length as number of atoms in fragment.
	coords: List of x,y,z coordinates for each atom in the fragment.
	
	Returns
	-------
	String formatted in PDB format.
	"""
	
	# asserts
	assert(atoms is not None and type(atoms) == list)
	assert(seq is not None and len(seq) >= 3)
	assert(chain is not None and type(chain) == list)
	assert(coords is not None and type(coords) == list)
	assert(len(atoms) == len(chain) == len(seq))
	
	# converting
	aa_dict = conv.protein_letters_1to3_extended
	
	# residue id
	res_id = 0
	prev   = None
	
	pdb_str     = ''
	connect_str = ''
	
	# printing the pdb
	for i, (a, s, ch, coo) in enumerate(zip(atoms, seq, chain, coords)):
		if s != prev:
			res_id += 1
			prev   = s
		
		x, y, z = coo
		
		begin = f'ATOM  {i+1:>5} {a:<4} {aa_dict[s].upper():>3} {ch}{res_id:>4}    '
		coordinates = f'{x:>8.3f}{y:>8.3f}{z:>8.3f}'
		end = '  1.00  0.00           C  '
		pdb_str += begin + coordinates + end + '\n'
		if i + 1 != len(atoms):
			connect = f'CONECT{(i+1):>5}{(i+2):>5}\n'
			connect_str += connect
	return pdb_str + connect_str


def layers_list(dic, template):
	"""
	Unroll dictionary that defines an encoder/decoder into a list of
	object containers for parameters in each layer.
	`layers_list` enforces that all values in dictionary are the same size.
	
	Input
	-----
	dic: Dictionary describing network layers.
	template: template dictionary with default values for parameters.
	
	Returns
	-------
	layers: list of objects for parameters per layer
	"""
	for v1 in dic.values():
		for v2 in dic.values(): assert(len(v2) == len(v1))
	
	layers = []
	kk = list(dic.keys())
	size = len(dic[kk[0]])
	for ii in range(size):
		new = {}
		for k, v in dic.items():
			if k not in new: new[k] = None
			new[k] = v[ii]
		newobj = Dict2Obj(new, template)
		layers.append(newobj)
	
	return layers


def conv_pool_out(hin, win, ksize=None, padding=None, stride=None):
	"""
	Compute resulting size of feature maps after convolution/pooling.
	
	Parameters
	----------
	* Positional args
		hin: height of input feature matrix
		win: width of input feature matrix
	* Keyword args
		ksize:   tuple for kernel dimensions
		padding: tuple of padding size for both dimensions
		stride:  tuple for stride step in each dimensions
	
	Returns
	------
	hout: resulting height of feature matrix
	wout: resulting width of feature matrix
	"""
	assert(hin is not None and win is not None)
	assert (ksize is not None and padding is not None and stride is not None)
	assert(type(hin) == int and type(win) == int)
	assert(
		type(ksize) == tuple and
		type(padding) == tuple and
		type(stride) == tuple)
	
	for k, p, s in zip(ksize, padding, stride):
		assert(type(k) == int and type(p) == int and type(s) == int)
	
	hout = hin + 2 * padding[0] - ksize[0]
	hout /= stride[0]
	hout += 1
	
	wout = win + 2 * padding[1] - ksize[1]
	wout /= stride[1]
	wout += 1
	
	if hout.is_integer() and wout.is_integer  : return int(hout), int(wout)
	else                                      : return None, None


def convt_out(hin, win, ksize=None, padding=None, stride=None):
	"""
	Compute resulting size of features after convolution transpose.
	
	Parameters
	----------
	* Positional args
		hin: height of input feature matrix
		win: width of input feature matrix
	* Keyword args
		ksize:   tuple for kernel dimensions
		padding: tuple of padding size for both dimensions
		stride:  tuple for stride step in each dimensions
	
	Returns
	------
	hout: resulting height of feature matrix
	wout: resulting width of feature matrix 
	"""
	assert(hin is not None and win is not None)
	assert (ksize is not None and padding is not None and stride is not None)
	assert(type(hin) == int and type(win) == int)
	assert(
		type(ksize) == tuple and
		type(padding) == tuple and
		type(stride) == tuple)
	
	for k, p, s in zip(ksize, padding, stride):
		assert(type(k) == int and type(p) == int and type(s) == int)
	
	hout = (hin - 1) * stride[0] - 2 * padding[0] + ksize[0]
	wout = (win - 1) * stride[1] - 2 * padding[1] + ksize[1]
	
	return hout, wout


def encoder_validator(hin, win, layers=None):
	"""
	Validate proposed encoder CNN layers produce valid dimensions for resulting
	feature maps. 
	
	Parameters
	----------
	* Positional args
		hin: height of input feature matrix
		win: width of input feature matrix
	* Keyword args
		layers: list of object containers holding parameters for each layer.
	
	Returns
	-------
	True/False if validation is successful or not
	(h, w) height and width of resulting feature map 
	"""
	assert(layers is not None)
	assert(type(hin) == int and type(win) == int)
	
	print(hin, win)
	h = hin
	w = win
	for l in layers:
		h, w = conv_pool_out(h, w, 
			ksize=l.conv_ks,
			padding=l.conv_paddings,
			stride=l.conv_strides)
		print(h, w)
		if h < 0 or w < 0: return False, None
		h, w = conv_pool_out(h, w,
			ksize=l.pool_ks,
			padding=l.pool_paddings,
			stride=l.pool_strides)
		print(h, w)
		if h < 0 or w < 0: return False, None
	
	return True, (h, w)


def decoder_validator(hin, win, layers=None):
	"""
	Validate proposed decoder CNN layers produce valide dimensions for resulting
	feature maps.
	
	Parameters
	----------
	* Positional args
		hin: height of input feature matrix
		win: width of input feature matrix
	* Keyword args
		layers: list of object containers holding parameters for each layer.
	
	Returns
	-------
	True/False if validation is successful or not
	(h, w) height and width of resulting feature map
	"""
	assert(type(hin) == int and type(win) == int)
	assert(layers is not None)
	
	print(hin, win)
	h = hin
	w = win
	for l in layers:
		h, w, = convt_out(h, w,
			ksize=l.convt_ks,
			padding=l.convt_paddings,
			stride=l.convt_strides)
		print('cont',h, w)
		if h < 0 or w < 0: return False, None
	
	return True, (h, w)


def cnn_ae_validator(inshape=None, encoder=None, decoder=None):
	"""
	Validate full CNN model for autoencoder architectures. 
	
	Parameters:
	inshape: tuple of input matrix dimensions
	encoder: dictionary of lists describing parameters for convolutions in
		encoder half.
	decoder: dictionary of lists describing parameters for convolution
		transposes in decoder half. 
	
	Returns:
	True if network specified is a valid CNN and autoencoder. 
	False otherwise
	"""
	hin, win = inshape
	encoder_layers = layers_list(encoder, encoder_template)
	decoder_layers = layers_list(decoder, decoder_template)
	status, latent = encoder_validator(hin, win, layers=encoder_layers)
	if status:
		status, out = decoder_validator(
			latent[0], latent[1], layers=decoder_layers)
		if status:
			print(out)
			out = (int(out[0]), int(out[1]))
			if (int(hin), int(win)) == out: return True
			else:                           return False
	else: return False


if __name__ == '__main__':
	# Test validator functions
	convd = {
		'conv_ks'       : [(4,4), (4,4)],
		'pool_ks'       : [(4,4), (4,4)],
		'conv_paddings' : [(1,1), (1,1)],
		'pool_paddings' : [(1,1), (1,1)],
		'pool_strides'  : [(1,1), (1,1)]
	}
	convtd = {
		'convt_ks'       : [(3,3), (3,3)],
		'convt_strides'  : [(1,1), (1,1)]
	}
	
	print(cnn_ae_validator(inshape=(7,7), encoder=convd, decoder=convtd))
