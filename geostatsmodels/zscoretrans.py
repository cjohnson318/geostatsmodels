#!/usr/bin/env python
import scipy
import numpy as np
import scipy.stats

def cdf( d ):
    '''
    Input:  (d)    iterable, a data set
    Output: (f)    NumPy array with (bins) rows and two columns
                   the first column are values in the range of (d)
                   the second column are CDF values
                   alternatively, think of the columns as the
                   domain and range of the CDF function
            (finv) inverse of (f)
    ---------------------------------------------------------
    Calculate the cumulative distribution function
    and its inverse for a given data set
    '''
    # the number of data points
    N = float( len( d ) )
    # sorted array of data points
    xs = np.sort( d )
    # array of unique data points
    xu = np.unique( xs )
    # number of unique data points
    U = len( xu )
    # initialize an array of U zeros
    cdf = np.zeros((U))
    # for each unique data point..
    for i in range( U ):
        # count the number of points less than
        # this point, and then divide by the
        # total number of data points
    	cdf[i] = len( xs[ xs < xu[i] ] ) / N
    # f : input value --> output percentage describing
    # the number of points less than the input scalar
    # in the modeled distribution; if 5 is 20% into
    # the distribution, then f[5] = 0.20
    f = np.vstack((xu,cdf)).T
    # inverse of f
    # finv : input percentage --> output value that 
    # represents the input percentage point of the
    # distribution; if 5 is 20% into the distribution,
    # then finv[0.20] = 5
    finv = np.fliplr(f)
    return f, finv
    
   
def fit( d ):
    '''
    Input:  (d) NumPy array with two columns,
                a domain and range of a mapping
    Output: (f) function that interpolates the mapping d
    ----------------------------------------------------
    This takes a mapping and returns a function that
    interpolates values that are missing in the original
    mapping, and maps values outside the range* of the
    domain (d) to the maximum or minimum values of the
    range** of (d), respectively.
    ----------------------------------------------------
    *range - here, I mean "maximum minus minimum"
    **range - here I mean the output of a mapping
    '''
    x, y = d[:,0], d[:,1]
    def f(t):
        # return the minimum of the range
        if t <= x.min():
            return y[ np.argmin(x) ]
        # return the maximum of the range
        elif t >= x.max():
            return y[ np.argmax(x) ]
        # otherwise, interpolate linearly
        else:
            intr = scipy.interpolate.interp1d( x, y )
            return intr(t)
    return f
    
def to_norm( data ):
    '''
    Input  (data) 1D NumPy array of observational data
    Output (z)    1D NumPy array of z-score transformed data
           (inv)  inverse mapping to retrieve original distribution
    '''
    # look at the dimensions of the data
    dims = data.shape
    # if there is more than one dimension..
    if len( dims ) > 1:
        # take the third column of the second dimension
        z = data[:,2]
    # otherwise just use data as is
    else:
        z = data
    # grab the number of data points
    N = len( z )
    # grab the cumulative distribution function
    f, inv = cdf( z )
    # h will return the cdf of z
    # by interpolating the mapping f
    h = fit( f )
    # ppf will return the inverse cdf
    # of the standard normal distribution
    ppf = scipy.stats.norm(0,1).ppf
    # for each data point..
    for i in range( N ):
        # h takes z to a value in [0,1]
        p = h( z[i] )
        # ppf takes p (in [0,1]) to a z-score
        z[i] = ppf( p )
    # convert positive infinite values
    posinf = np.isposinf( z )
    z = np.where( posinf, np.nan, z )
    z = np.where( np.isnan( z ), np.nanmax( z ), z )
    # convert negative infinite values
    neginf = np.isneginf( z )
    z = np.where( neginf, np.nan, z )
    z = np.where( np.isnan( z ), np.nanmin( z ), z )
    # if the whole data set was passed, then add the
    # transformed variable and recombine with data
    if len( dims ) > 1:
        z = np.vstack(( data[:,:2].T, z )).T
    return z, inv

def from_norm( data, inv ):
    '''
    Input:  (data) NumPy array of zscore data
            (inv)  mapping that takes zscore data back
                   to the original distribution
    Output: (z)    Data that should conform to the 
                   distribution of the original data
    '''
    # convert to a column vector
    d = data.ravel()
    # create an interpolating function 
    # for the inverse cdf, mapping zscores
    # back to the original data distribution
    h = fit( inv )
    # convert z-score data to cdf values in [0,1]
    f = scipy.stats.norm(0,1).cdf( d )
    # use inverse cdf to map [0,1] values to the
    # original distribution, then add the mu and sd
    z = np.array( [ h(i) for i in f ] )
    # reshape the data
    z = np.reshape( z, data.shape )
    return z
