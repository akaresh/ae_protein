#!/usr/bin/python3 

"""
Shared library for functions for PyTorch NN training
"""

import numpy as np
from scipy.spatial.distance import cdist

def normalize_frag(frag):
	"""
	Normalize structural fragment coordinates
	"""
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
	
def distance_matrix(frag):
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
	Fit functions for all models
	
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
