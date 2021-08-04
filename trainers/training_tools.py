#!/usr/bin/python3

"""
Shared library for functions for PyTorch NN training
"""

import numpy as np
from scipy.spatial.distance import cdist
import Bio.Data.IUPACData as conv


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
			test_loss = criterion(outputs, bv)
			vloss += test_loss.item()
		vloss = vloss / len(test)
		
		# display the epoch training loss
		print(f"epoch : {epoch+1}/{epochs} recon loss = {loss:.8f} \
		test loss = {vloss:.8f}")
	return model


def pdb_writer(atoms=None, seq=None, chain=None, coords=None):
	import sys
	print('atoms', atoms)
	print('seq', seq)
	print('chain', chain)
	print('coords', coords)

	#asserts
	assert(atoms!=None and type(atoms) == list)
	assert(seq!=None and len(seq) > 3)
	assert(chain!=None and type(chain) == list)
	assert(coords!=None and type(coords) == list)
	assert(len(atoms) == len(chain) == len(seq))


	###asserts and checks
	#check atoms, seq, chain, coords
	#if seq is type str, length < 4

	#converting
	aa_dict = conv.protein_letters_1to3_extended

	#residue id
	res_id = 0
	prev = None

	pdb_str = ''
	connect_str = ''

	#printing the pdb
	for i, (a, s, ch, coo) in enumerate(zip(atoms, seq, chain, coords)):

		if s != prev:
			res_id += 1
			prev = s

		x, y, z = coo

		begin = f'ATOM  {i+1:>5} {a:<4} {aa_dict[s].upper():>3} {ch}{res_id:>4}    '
		coordinates = f'{x:>8.3f}{y:>8.3f}{z:>8.3f}'
		end = '  1.00  0.00           C  '
		pdb_str += begin+coordinates+end+'\n'
		if i+1 != len(atoms):
			connect = f'CONECT{(i+1):>5}{(i+2):>5}\n'
			connect_str += connect
	print(pdb_str+connect_str)
	sys.exit()
	return True
