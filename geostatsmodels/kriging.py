#!/usr/bin/env python

import numpy as np
from pandas import DataFrame, Series
from scipy.spatial.distance import pdist, cdist, squareform
#from numba import jit

def set_of_points_at_lag_h( data, lag, tol ):
    '''
    Input:  (data) NumPy array where the fris t two columns
                   are the spatial coordinates, x and y
            (lag)  the distance, h, between points
            (tol)  the tolerance we are comfortable with around (lag)
    Output: (ind)  list of tuples; the first element is the row of
                   (data) for one point, the second element is the row 
                   of a point (lag)+/-(tol) away from the first point,
                   e.g., (3,5) corresponds fo data[3,:], and data[5,:]
    '''
    # create a distance matrix
    p = squareform( pdist( data[:,:2] ) )
    # grab the coordinates in a given range: lag +/- tolerance
    i, j = np.where( ( p >= lag - tol )&( p <= lag + tol ) )
    # zip the coordinates into a list
    ind = zip( i, j )
    # take out the repeated elements,
    # since p is a *symmetric* distance matrix
    ind = [ i for i in ind if i[1] > i[0] ]
    return ind

def semivariogram_at_lag_h( data, lag, tol ):
    '''
    Input:  (data) NumPy array where the fris t two columns
                   are the spatial coordinates, x and y
            (lag)  the distance, h, between points
            (tol)  the tolerance we are comfortable with around (lag)
    Output:  (z)   semivariogram value at lag (h) +/- (tol)  
    '''
    # grab the indices of the points
    # that are lag +/- tolerance apart
    ind = set_of_points_at_lag_h( data, lag, tol )
    # take the squared difference between
    # the values of the variable of interest
    z = [ ( data[i,2] - data[j,2] )**2.0 for i,j in ind ]
    # half the mean squared difference
    z = np.mean( z ) / 2.0
    # return the semivariogram
    return z

def semivariogram( data, lags, tol ):
    '''
    Input:  (data) NumPy array where the fris t two columns
                   are the spatial coordinates, x and y
            (lag)  the distance, h, between points
            (tol)  the tolerance we are comfortable with around (lag)
    Output: (sv)   <2xN> NumPy array of lags and semivariogram values
    '''
    # calculate the semivarigram at different lags given some tolerance
    sv = [ semivariogram_at_lag_h( data, lag, tol ) for lag in lags ]
    # bundle the semivariogram values with their lags
    return np.array( zip( lags, sv ) ).T
 
def C( P, lag, tol ):
    '''
    Calculate the sill
    '''
    c0 = np.var( P[:,2] )
    if h == 0:
        return c0
    return c0 - semivariogram_at_lag_h( P, lag, tol )

def opt( fct, x, y, C0, parameterRange=None, meshSize=1000 ):
    '''
    Optimize parameters for a model of the semivariogram
    '''
    if parameterRange == None:
        parameterRange = [ x[1], x[-1] ]
    mse = np.zeros( meshSize )
    a = np.linspace( parameterRange[0], parameterRange[1], meshSize )
    for i in range( meshSize ):
        mse[i] = np.mean( ( y - fct( x, a[i], C0 ) )**2.0 )
    return a[ mse.argmin() ]
    
def spherical( h, a, C0 ):
    '''
    Spherical model of the semivariogram
    '''
    # if h is a single digit
    if type(h) == np.float64:
        # calculate the spherical function
        if h <= a:
            return C0*( 1.5*h/a - 0.5*(h/a)**3.0 )
        else:
            return C0
    # if h is an iterable
    else:
        # calcualte the spherical function for all elements
        a = np.ones( h.size ) * a
        C0 = np.ones( h.size ) * C0
        return map( spherical, h, a, C0 )
    
def exponential( h, a, C0 ):
    '''
    Exponential model of the semivariogram
    '''
    # if h is a single digit
    if type(h) == np.float64:
        # calculate the exponential function
        return C0*( 1.0 - np.exp( -3.0 * h / a ) )
    # if h is an iterable
    else:
        # calcualte the exponential function for all elements
        a = np.ones( h.size ) * a
        C0 = np.ones( h.size ) * C0
        return map( exponential, h, a, C0 )
        
def gaussian( h, a, C0 ):
    '''
    Gaussian model of the semivariogram
    '''
    # if h is a single digit
    if type(h) == np.float64:
        # calculate the Gaussian function
        return C0*( 1.0 - np.exp( -3.0 * h**2.0 / a**2.0 ) )
    # if h is an iterable
    else:
        # calcualte the Gaussian function for all elements
        a = np.ones( h.size ) * a
        C0 = np.ones( h.size ) * C0
        return map( gaussian, h, a, C0 )
         
def cvmodel( P, model, lags, tol ):
    '''
    Input:  (P)      ndarray, data
            (model)  modeling function
                      - spherical
                      - exponential
                      - gaussian
            (lags)   lag distances
            (tol)    tolerance
    Output: (covfct) function modeling the covariance
    '''
    # calculate the semivariogram
    sv = semivariogram( P, lags, tol )
    # calculate the sill
    C0 = C( P, lags[0], tol )
    # calculate the optimal parameters
    param = opt( model, sv[0], sv[1], C0 )
    # return a covariance function
    covfct = lambda h, a=param: C0 - model( h, a, C0 )
    return covfct

def krige( P, model, lags, tol, u, N ):
    '''
    Input  (P)     ndarray, data
           (model) modeling function
                    - spherical
                    - exponential
                    - gaussian
           (lags)  kriging lag distances
           (tol)   kriging tolerance
           (u)     unsampled point
           (N)     number of neighboring 
                   points to consider
    '''

    # covariance function
    covfct = cvmodel( P, model, lags, tol )
    # mean of the variable
    mu = np.mean( P[:,2] )

    # distance between u and each data point in P
	d = cdist( P[:,:2], u )
    # add these distances to P
    P = np.vstack(( P.T, d )).T
    # sort P by these distances
    # take the first N of them
    P = P[d.argsort()[:N]]

    # apply the covariance model to the distances
    k = covfct( P[:,3] )
    # cast as a matrix
    k = np.matrix( k ).T

    # form a matrix of distances between existing data points
    K = squareform( pdist( P[:,:2] ) )
    # apply the covariance model to these distances
    K = covfct( K.ravel() )
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
    
