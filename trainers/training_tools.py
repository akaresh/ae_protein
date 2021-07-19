#!/usr/bin/python3 

"""
Shared library for functions for PyTorch NN training

"""

import numpy as np

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