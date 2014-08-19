#!/usr/bin/env python

import numpy as np
from geostatsmodels import variograms

def hscattergram( data, pwdist, lag, tol ):
    '''
    Input:  (data)    NumPy array with three columns, the first two 
                      columns should be the x and y coordinates, and 
                      third should be the measurements of the variable
                      of interest
            (lag)     the lagged distance of interest
            (tol)     the allowable tolerance about (lag)
            (pwdist)  a square pairwise distance matrix
    Output:           h-scattergram figure showing the distribution of
                      measurements taken at a certain lag and tolerance
    '''
    # calculate the pairwise distances
    indices = variograms.lagindices( pwdist, lag, tol )
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

def laghistogram( data, lags, tol, pwdist=None ):
    '''
    Input:  (data)    NumPy array with three columns, the first two 
                      columns should be the x and y coordinates, and 
                      third should be the measurements of the variable
                      of interest
            (lags)    the lagged distance of interest
            (tol)     the allowable tolerance about (lag)
            (pwdist)  the pairwise distances, optional
    Output:           lag histogram figure showing the number of
                      distances at each lag
    '''
    # calculate the pairwise distances if they are not given
    if pwdist == None:
        pwdist = pairwise( data[:,:2] )
    # collect the distances at each lag
    indices = [ variograms.lagindices( pwdist, lag, tol ) for lag in lags ]
    # record the number of indices at each lag
    indices = [ len( i ) for i in indices ]
    # create a bar plot
    fig, ax = subplots()
    ax.bar( lags+tol, indices )
    ax.set_ylabel('Number of Lags')
    ax.set_xlabel('Lag Distance')
    ax.set_title('Lag Histogram')
    show();
