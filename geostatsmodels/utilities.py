#!/usr/bin/env python

import scipy
import scipy.stats
import numpy as np
from scipy.spatial.distance import pdist, squareform

def readGeoEAS( fn ):
    '''
    Input:  (fn)   filename describing a GeoEAS file
    Output: (data) NumPy array
    --------------------------------------------------
    Read GeoEAS files as described by the GSLIB manual
    '''
    f = open( fn, "r" )
    # title of the data set
    title = f.readline()
    # number of variables
    nvar = int( f.readline() )
    # variable names
    columns = [ f.readline().strip() for i in range( nvar ) ]
    # make a list for the data
    data = list()
    # for each line of the data
    while True:
        # read a line
        line = f.readline()
        # if that line is empty
        if line == '':
            # the jump out of the while loop
            break
        # otherwise, append that line to the data list
        else:
            data.append( line )
    # strip off the newlines, and split by whitespace
    data = [ i.strip().split() for i in data ]
    # turn a list of list of strings into an array of floats
    data = np.array( data, dtype=np.float )
    # combine the data with the variable names into a DataFrame
    # df = pandas.DataFrame( data, columns=columns,  )
    return data
    
def pairwise( data ):
    '''
    Input:  (data) NumPy array where the first two columns
                   are the spatial coordinates, x and y
    '''
    # determine the size of the data
    npoints, cols = data.shape
    # give a warning for large data sets
    if npoints > 10000:
        print("You have more than 10,000 data points, this might take a minute.")
    # return the square distance matrix
    return squareform( pdist( data[:,:2] ) )

def degree_to_bearing( deg ):
    bearing = None
    if deg == None:
        raise ValueError('Input "deg" cannot be "None"')
    if( deg >= 0 )&( deg <= 90 ):
        bearing = 90 - deg
    elif( deg >= 90 )&( deg < 180 ):
        bearing = 90 + 360 - deg
    elif( deg >= 180 )&( deg < 270 ):
        bearing = 90 + 360 - deg
    elif( deg >= 270 )&( deg <= 360 ):
        bearing = 90 + 360 - deg
    return bearing

def bearing( p0, p1 ):
    '''
    Input:  (p0,p1) two iterables representing x and y coordinates
    Output:         the bearing, a degree in the range [0,360) from
                    due North clockwise aroung the compass face
    '''
    y = p1[1] - p0[1]
    x = p1[0] - p0[0]
    rad, deg = None, None
    if x != 0:
        rad = np.arctan( y / float( x ) )
        deg = np.rad2deg( rad )
    else:
        if y > 0:
            deg = 90
        elif y < 0:
            deg = 270
    if( x < 0 )&( y > 0 ):
        deg = 180 + deg
    elif( x < 0 )&( y < 0 ):
        deg = 270 - deg
    elif( x > 0 )&( y < 0 ):
        deg = 360 + deg
    return degree_to_bearing( deg )

def bearings( data, indices ):
    '''
    Input:  (data)    the original data set: x, y, variable
            (indices) the output of variograms.lagindices()
    Output:           a list of list of bearings corresponding
                      to variograms.lagindices()
    '''
    # indices is a list of lists
    # each list in indices represents a set of points at a certain lag
    # for a given lag, we have a list of index pairs
    # each pair represents the indices of two points in the data set having a given lag
    # --deep breath--
    # for each lag-set, we go through all of the index-pairs, and calculate the bearing
    # in the end, instead of a list of lag-sets containing index-pairs,
    # we have a list of lag-sets containing bearings
    return [ bearing( *data[ idx ] ) for idx in indices ]

def inangle( theta, angle, atol ):
    '''
    Input:  (theta) angle inquestion
            (angle) reference angle
            (atol)  tolerance about (angle)
    Output:         True or False, depending on whether
                    theta is in [angle-atol,angle+atol)
    '''
    lower = angle - atol
    upper = angle + atol
    if( lower >= 0 )&( upper <= 360 ):
        if( theta >= lower )&( theta < upper ):
            return True
    if( lower < 0 ):
        lower %= 360
        if( theta > upper )&( theta >= lower ):
            return True
        elif( theta >= 0 )&( theta < upper ):
            return True
    if( upper > 360 ):
        upper %= 360
        if( theta >= lower )&( theta < 0 ):
            return True
        if( theta >= 0 )&( theta < upper ):
            return True
    return False
