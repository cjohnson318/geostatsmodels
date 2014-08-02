#!/usr/bin/env python

import numpy as np
import geostatsmodels.kriging as k
import random

def gridpath( xdim, ydim ):
    '''
    Input:  (xdim) iterable describing the start, stop, and no. steps
            (ydim) iterable describing the start, stop, and no. steps
    Output: (path) path through the 2D grid, a list of lists
                   each element has an index, an address in the gird,
                   and an address in space
    '''
    # dim = ( start, stop, steps )
    xrng = np.linspace( *xdim )
    yrng = np.linspace( *ydim )
    # total number of steps in the random path
    N = xdim[2] * ydim[2]
    # an array of indices
    idx = np.arange( N )
    # shuffle the indices
    random.shuffle( idx )
    # create a list for the path
    path = list()
    # create a counter
    t = 0
    # for each cell in the x dimension
    for i in range( xdim[2] ):
    	# for each cell in the y dimension
        for j in range( ydim[2] ):
            # record a shuffled index value, idx[t],
            # an integer cell address (i,j), and a
            # physical address in feet, ( xrng[i], yrng[j] )
            path.append( [ idx[t], (i,j), (xrng[i],yrng[j]) ] )
            # increment t for the shuffled indices, idx
	    t += 1
    # sort the shuffled indices
    # thereby shuffling the path
    path.sort()
    return path

def sgs( data, model, hs, bw, xs, ys=None, pad=0.0 ):
    '''
    Input:  (data)   <N,3> NumPy array of data
            (hs)     NumPy array of distances
            (bw)     bandwidth of the semivariogram
            (xs)     number of cells in the x dimension
            (ys)     number of cells in the y dimension
    Output: (M)      <xsteps,ysteps> NumPy array of data
                     representing the simulated distribution
                     of the variable of interest 
    '''
    # check for meshsize in second dimension
    if ys == None:
	      ys = xs
    # create path
    xdim = ( data[:,0].min()-pad, data[:,0].max()+pad, xs )
    ydim = ( data[:,1].min()-pad, data[:,1].max()+pad, ys )
    path = gridpath( xdim, ydim )
    # create array for the output
    M = np.zeros((xs,ys))
    # for each cell in the grid..
    for step in path :
        # grab the index, the cell address, and the physical location
        idx, cell, loc = step
        # perform the kriging
        kv = k.krige( data, model, hs, bw, loc, 4 )
        # add the kriging estimate to the output
        M[cell[0],cell[1]] = kv
        # add the kriging estimate to a spatial location
        newdata = [ loc[0], loc[1], kv ]
        # add this new point to the data used for kriging
        data = np.vstack(( data, newdata ))
    return M
