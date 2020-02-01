#!/usr/bin/env python

import numpy as np
from scipy.spatial.distance import cdist

def moran_i(data, w=None):
	
	n = len(data)
	values = data[:,2]
	
	# compute distances
	distances = cdist( data[:,:2], data[:,:2] )
	
	# create weights matrix
	weights = 1 / distances
	np.fill_diagonal(weights, 0)
	
	# or use weights submitted as parameter
	weights = w or weights
	
	# normalize weights
	rowsum = np.sum(weights, axis=0)
	rowsum[rowsum == 0] = 1
	weights = weights / rowsum

	nom   = 0.0
	denom = 0.0
	mean  = np.mean(values)
	
	# normalize weights
	rowsum = np.sum(weights, axis=0)
	rowsum[rowsum == 0] = 1
	weights = weights / rowsum
	
	# 
	y = values - mean
	z = np.outer( y , y )
	nom = np.sum( weights * z )

	denom = ( values - mean ) * ( values - mean )
	denom = np.sum( denom )
	
	S0 = np.sum(weights)
	
	nom   = nom * n
	denom = denom * S0
	
	I = nom / denom
	
	# expected E[I] value
	E = -1 / (n - 1)
	
	return I, E

def geary_c(data, w=None):
	
	n = len(data)
	values = data[:,2]
	
	# compute distances
	distances = cdist( data[:,:2], data[:,:2] )
	
	# create weights matrix
	weights = 1 / distances
	np.fill_diagonal(weights, 0)
	
	# or use weights submitted as parameter
	weights = w or weights
	
	# normalize weights
	rowsum = np.sum(weights, axis=0)
	rowsum[rowsum == 0] = 1
	weights = weights / rowsum

	nom   = 0.0
	denom = 0.0
	mean  = np.mean(values)
	
	# compute numerator's double sum
	for i in range(n):
		for j in range(n):
			nom = nom + weights[i,j] * ( values[i] - values[j] ) * ( values[i] - values[j] )

	denom = ( values - mean ) * ( values - mean )
	denom = np.sum( denom )
	
	S0 = np.sum(weights)
	
	nom   = nom * ( n - 1 ) 
	denom = denom * S0 * 2
	
	C = nom / denom
	
	return C
