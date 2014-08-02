#!/usr/bin/env python
import numpy as np
from scipy.spatial.distance import cdist
from geostatsmodels.utilities import pairwise

def krige( data, covfct, lags, tol, u, N=0 ):
    '''
    Input  (P)     ndarray, data
           (model) modeling function
                    - spherical
                    - exponential
                    - gaussian
           (lags)  kriging lag distances
           (tol)   kriging tolerance
           (u)     unsampled point
           (N)     number of neighboring points
                   to consider, if zero use all
    '''

    # mean of the variable
    mu = np.mean( data[:,2] )

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

    # calculate the kriging weights
    weights = np.linalg.inv( K ) * k
    weights = np.array( weights )

    # calculate the residuals
    residuals = P[:,2] - mu

    # calculate the estimation
    estimation = np.dot( weights.T, residuals ) + mu

    return float( estimation )
