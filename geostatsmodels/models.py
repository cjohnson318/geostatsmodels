import numpy as np
import variograms

def opt( fct, x, y, c, parameterRange=None, meshSize=1000 ):
    '''
    Optimize parameters for a model of the semivariogram
    '''
    if parameterRange == None:
        parameterRange = [ x[1], x[-1] ]
    mse = np.zeros( meshSize )
    a = np.linspace( parameterRange[0], parameterRange[1], meshSize )
    for i in range( meshSize ):
        mse[i] = np.mean( ( y - fct( x, a[i], c ) )**2.0 )
    return a[ mse.argmin() ]

def spherical( h, a, c ):
    '''
    Spherical model of the semivariogram
    '''
    # if h is a single digit
    if type(h) == np.float64:
        # calculate the spherical function
        if h <= a:
            return c*( 1.5*h/a - 0.5*(h/a)**3.0 )
        else:
            return c
    # if h is an iterable
    else:
        # calcualte the spherical function for all elements
        a = np.ones( h.size ) * a
        c = np.ones( h.size ) * c
        return map( spherical, h, a, c )

def exponential( h, a, c ):
    '''
    Exponential model of the semivariogram
    '''
    # if h is a single digit
    if type(h) == np.float64:
        # calculate the exponential function
        return c*( 1.0 - np.exp( -3.0 * h / a ) )
    # if h is an iterable
    else:
        # calcualte the exponential function for all elements
        a = np.ones( h.size ) * a
        c = np.ones( h.size ) * c
        return map( exponential, h, a, c )

def gaussian( h, a, c ):
    '''
    Gaussian model of the semivariogram
    '''
    # if h is a single digit
    if type(h) == np.float64:
        # calculate the Gaussian function
        return c*( 1.0 - np.exp( -3.0 * h**2.0 / a**2.0 ) )
    # if h is an iterable
    else:
        # calcualte the Gaussian function for all elements
        a = np.ones( h.size ) * a
        c = np.ones( h.size ) * c
        return map( gaussian, h, a, c )

def covmodel( data, model, lags, tol ):
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
    sv = variograms.semivariogram( data, lags, tol )
    # calculate the sill
    c = np.var( data[:,2] )
    # calculate the optimal parameters
    param = opt( model, sv[0], sv[1], c )
    # return a covariance function
    covfct = lambda h, a=param: c - model( h, a, c )
    return covfct
