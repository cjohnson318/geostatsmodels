#!/usr/bin/env python
import numpy as np
from scipy.spatial.distance import cdist
from geostatsmodels.utilities import pairwise

def kmatrices( data, covfct, u, N=0 ):
    '''
    Input  (data)  ndarray, data
           (model) modeling function
                    - spherical
                    - exponential
                    - gaussian
           (u)     unsampled point
           (N)     number of neighboring points
                   to consider, if zero use all
    '''
    # u needs to be two dimensional for cdist()
    if np.ndim( u ) == 1:
        u = [u]
    # distance between u and each data point in P
    d = cdist( data[:,:2], u )
    # add these distances to P
    P = np.hstack(( data, d ))
    # if N>0, take the N closest points,
    if N > 0:
        P = P[d[:,0].argsort()[:N]]
    # otherwise, use all of the points
    else:
        N = len( P )

    # apply the covariance model to the distances
    k = covfct( P[:,3] )
    # check for nan values in k
    if np.any( np.isnan( k ) ):
        raise ValueError('The vector of covariances, k, contains NaN values')
    # cast as a matrix
    k = np.matrix( k ).T

    # form a matrix of distances between existing data points
    K = pairwise( P[:,:2] )
    # apply the covariance model to these distances
    K = covfct( K.ravel() )
    # check for nan values in K
    if np.any( np.isnan( K ) ):
        raise ValueError('The matrix of covariances, K, contains NaN values')
    # re-cast as a NumPy array -- thanks M.L.
    K = np.array( K )
    # reshape into an array
    K = K.reshape(N,N)
    # cast as a matrix
    K = np.matrix( K )

    return K, k, P

def simple( data, covfct, u, N=0, nugget=0 ):
    
    # calculate the matrices K, and k
    K, k, P = kmatrices( data, covfct, u, N )

    # calculate the kriging weights
    weights = np.linalg.inv( K ) * k
    weights = np.array( weights )

    # calculate k' * K * k for
    # the kriging variance
    kvar = k.T * weights

    # mean of the variable
    mu = np.mean( data[:,2] )
    
    # calculate the residuals
    residuals = P[:,2] - mu

    # calculate the estimation
    estimation = np.dot( weights.T, residuals ) + mu

    # calculate the sill and the 
    # kriging standard deviation
    sill = np.var( data[:,2] )
    kvar = float( sill + nugget - kvar )
    kstd = np.sqrt( kvar )

    return float( estimation ), kstd

def ordinary( data, covfct, u, N=0, nugget=0 ):

    # calculate the matrices K, and k
    Ks, ks, P = kmatrices( data, covfct, u, N )

    # if N is not set, determine from Ks
    if N == 0:
        N, N = Ks.shape

    # add a column and row of ones to Ks,
    # with a zero in the bottom, right hand corner
    K = np.matrix( np.ones(( N+1, N+1 )) )
    K[:N,:N] = Ks
    K[-1,-1] = 0.0

    # add a one to the end of ks
    k = np.matrix( np.ones(( N+1,1 )) )
    k[:N] = ks
    
    # calculate the kriging weights
    weights = np.linalg.inv( K ) * k
    weights = np.array( weights )

    # calculate k' * K * k for
    # the kriging variance
    kvar = k.T * weights

    # mean of the variable
    mu = np.mean( data[:,2] )
    
    # calculate the residuals
    residuals = P[:,2]

    # calculate the estimation
    estimation = np.dot( weights[:-1].T, residuals )

    # calculate the sill and the kriging standard deviation
    sill = np.var( data[:,2] )
    kvar = float( sill + nugget - kvar )
    kstd = np.sqrt( kvar )

    return float( estimation ), kstd

def krige(data, covfct, grid, method='simple', N=0, nugget=0):
    '''
    Krige an <Nx2> array of points representing a grid.
    
    Use either simple or ordinary kriging, some number N
    of neighboring points, and a nugget value.
    '''
    kriging_method = None
    if method == 'simple':
        kriging_method = simple 
    elif method == 'ordinary':
        kriging_method = ordinary
    M = len(grid)
    est = np.zeros((M, 1))
    kstd = np.zeros((M, 1))
    for i in range(M):
        est[i], kstd[i] = kriging_method(data, covfct, grid[i], N, nugget)
    return est, kstd
