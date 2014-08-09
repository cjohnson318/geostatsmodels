#!/usr/bin/env python
import scipy
import numpy as np
import scipy.stats
import matplotlib
from pylab import *
from scipy.spatial.distance import pdist, squareform
import variograms

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
        print "You have more than 10,000 data points, this might take a minute."
    # return the square distance matrix
    return squareform( pdist( data[:,:2] ) )
    
def hscattergram( data, lag, tol, pwdist=None ):
	'''
	Input:  (data)    NumPy array with three columns, the first two 
	                  columns should be the x and y coordinates, and the 
	                  third should be the measurements of the variable
	                  of interest
	        (lag)     the lagged distance of interest
	        (tol)     the allowable tolerance about (lag)
	        (pwdist)  a square pairwise distance matrix
	Output:           h-scattergram figure showing the distribution of
	                  measurements taken at a certain lag and tolerance    
	'''
	# calculate the pairwise distances if they are not given
	if pwdist == None:
		pwdist = pairwise( data )
		indices = variograms.lagindices( pwdist, lag, tol )
	else:
		indices = variograms.lagindices( data, lag, tol )
	# collect the head and tail measurements
	head = data[ indices[:,0], 2 ]
	tail = data[ indices[:,1], 2 ]
	# create a scatterplot with equal axes
	fig, ax = subplots()
	ax.scatter( head, tail, marker="o", facecolor="none", edgecolor="b", alpha=0.5 );
	ax.set_aspect("equal");
	# set the labels and the title
	ax.set_ylabel("$z(u+h)$");
	ax.set_xlabel("$z(u)$");
	ax.set_title("Lags Between "+str(lag-tol)+" and "+str(lag+tol))
	# grab the limits of the axes
	xmin, xmax = ax.get_xlim();
	ymin, ymax = ax.get_ylim();
	# calculate the covariance and annotate
	cv = variograms.covariance( data, indices );
	ax.text( xmin*1.25, ymin*1.050, 'Covariance = {:3.2f}'.format(cv) );
	# calculate the semivariance and annotate
	sv = variograms.semivariance( data, indices );
	ax.text( xmin*1.25, ymin*1.025, 'Semivariance = {:3.2f}'.format(sv) );
	show();

# this is a colormap that ranges from yellow to purple to black
cdict = {'red':   ((0.0, 1.0, 1.0),
                   (0.5, 225/255., 225/255. ),
                   (0.75, 0.141, 0.141 ),
                   (1.0, 0.0, 0.0)),
         'green': ((0.0, 1.0, 1.0),
                   (0.5, 57/255., 57/255. ),
                   (0.75, 0.0, 0.0 ),
                   (1.0, 0.0, 0.0)),
         'blue':  ((0.0, 0.376, 0.376),
                   (0.5, 198/255., 198/255. ),
                   (0.75, 1.0, 1.0 ),
                   (1.0, 0.0, 0.0)) }
                   
YPcmap = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict,256)
