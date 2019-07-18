import numpy as np
import geostatsmodels.variograms as variograms

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
    
def typetest( h, a, lta, gta ):
    '''
    Input:  (h)   scalar or NumPy ndarray
            (a)   scalar representing the range parameter
            (lta) function to perfrom for values less than (a)
            (gta) function to perform for values greater than (a)
    Output:       scalar or array, depending on (h)
    '''
    # if (h) is a numpy ndarray, then..
    try:
        # apply lta() to elements less than a
        lt = lta( h[ np.where( h <= a ) ] )
        # apply gta() to elements greater than a
        gt = gta( h[ np.where( h > a ) ] )
        return np.hstack((lt,gt))
    # otherwise, if (h) is a scalar..
    except TypeError:
        if h <= a:
            return lta( h )
        else:
            return gta( h )
        
def nugget( h, a, c ):
    '''
    Nugget model of the semivariogram
    '''
    c = float(c)
    lta = lambda x: 0+x*0
    gta = lambda x: c+x*0
    return typetest( h, 0, lta, gta )

def linear( h, a, c ):
    '''
    Linear model of the semivariogram
    '''
    a, c = float(a), float(c)
    lta = lambda x: (c/a)*x
    gta = lambda x: c+x*0
    return typetest( h, a, lta, gta )

def spherical( h, a, c ):
    '''
    Spherical model of the semivariogram
    '''
    a, c = float(a), float(c)
    lta = lambda x: c*( 1.5*(x/a) - 0.5*(x/a)**3.0 )
    gta = lambda x: c+x*0
    return typetest( h, a, lta, gta )

def exponential( h, a, c ):
    '''
    Exponential model of the semivariogram
    '''
    a, c = float( a ), float( c )
    return c*( 1.0 - np.exp( -3.0*h/a ) )

def gaussian( h, a, c ):
    '''
    Gaussian model of the semivariogram
    '''
    a, c = float( a ), float( c )
    return c*( 1.0 - np.exp( -3.0*h**2.0/a**2.0 ) )
    
def power( h, w, c ):
    '''
    Power model of the semivariogram
    '''
    return c*h**w

def semivariance( fct, param ): 
    '''
    Input:  (fct)   function that takes data and parameters
            (param) list or tuple of parameters
    Output: (inner) function that only takes data as input
                    parameters are set internally
    '''
    def inner( h ):
        return fct(h,*param)
    return inner
    
def covariance( fct, param ): 
    '''
    Input:  (fct)   function that takes data and parameters
            (param) list or tuple of parameters
    Output: (inner) function that only takes data as input
                    parameters are set internally
    '''
    def inner( h ):
        return param[-1] - fct(h,*param)
    return inner

def fitmodel( data, fct, lags, tol ):
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
    a = opt( fct, sv[0], sv[1], c )
    # return a covariance function
    covfct = covariance( fct, ( a, c ) )
    return covfct
