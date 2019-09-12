#!/usr/bin/env python
import numpy as np

import geostatsmodels.utilities as utilities


def lagindices(pwdist, lag, tol):
    '''
    Input:  (pwdist) square NumPy array of pairwise distances
            (lag)    the distance, h, between points
            (tol)    the tolerance we are comfortable with around (lag)
    Output: (ind)    list of tuples; the first element is the row of
                     (data) for one point, the second element is the row
                     of a point (lag)+/-(tol) away from the first point,
                     e.g., (3,5) corresponds fo data[3,:], and data[5,:]
    '''
    # grab the coordinates in a given range: lag +/- tolerance
    i, j = np.where((pwdist >= lag - tol) & (pwdist < lag + tol))
    # take out the repeated elements,
    # since p is a *symmetric* distance matrix
    indices=np.c_[i, j][np.where(j > i)]
    return indices


def anilagindices(data, pwdist, lag, tol, angle, atol):
    '''
    Input:  (data)   NumPy array where the frist two columns
                     are the spatial coordinates, x and y, and
                     the third column is the variable of interest
            (pwdist) square NumPy array of pairwise distances
            (lag)    the distance, h, between points
            (tol)    the tolerance we are comfortable with around (lag)
            (angle)  float, [0,360), North = 0 --> 360 clockwise
            (atol)   number of degrees about (angle) to consider
    '''
    index = lagindices(pwdist, lag, tol)
    brngs = utilities.bearings(data, index)
    bridx = list(zip(brngs, index))
    index = [idx.tolist() for br, idx in bridx if utilities.inangle(br, angle, atol)]
    # add 180 to the angle and take the modulus
    angle = (angle + 180) % 360
    index += [idx.tolist() for br, idx in bridx if utilities.inangle(br, angle, atol)]
    return np.array(index)


def semivariance(data, indices):
    '''
    Input:  (data)    NumPy array where the first two columns
                      are the spatial coordinates, x and y, and
                      the third column is the variable of interest
            (indices) indices of paired data points in (data)
    Output:  (z)      semivariance value at lag (h) +/- (tol)
    '''
    # take the squared difference between
    # the values of the variable of interest
    # the semivariance is half the mean squared difference
    i=indices[:, 0]
    j=indices[:, 1]
    z=(data[i, 2] - data[j, 2])**2.0
    return np.mean(z) / 2.0


def semivariogram(data, lags, tol):
    '''
    Input:  (data) NumPy array where the first two columns
                   are the spatial coordinates, x and y
            (lag)  the distance, h, between points
            (tol)  the tolerance we are comfortable with around (lag)
    Output: (sv)   <2xN> NumPy array of lags and semivariogram values
    '''
    return variogram(data, lags, tol, 'semivariogram')


def covariance(data, indices):
    '''
    Input:  (data) NumPy array where the first two columns
                   are the spatial coordinates, x and y
            (lag)  the distance, h, between points
            (tol)  the tolerance we are comfortable with around (lag)
    Output:  (z)   covariance value at lag (h) +/- (tol)
    '''
    # grab the indices of the points
    # that are lag +/- tolerance apart
    i=indices[:, 0]
    j=indices[:, 1]
    return np.cov(data[i, 2], data[j, 2])[0][1]


def covariogram(data, lags, tol):
    '''
    Input:  (data) NumPy array where the first two columns
                   are the spatial coordinates, x and y
            (lag)  the distance, h, between points
            (tol)  the tolerance we are comfortable with around (lag)
    Output: (cv)   <2xN> NumPy array of lags and covariogram values
    '''
    return variogram(data, lags, tol, 'covariogram')


def variogram(data, lags, tol, method):
    '''
    Input:  (data) NumPy array where the first two columns
                   are the spatial coordinates, x and y
            (lag)  the distance, h, between points
            (tol)  the tolerance we are comfortable with around (lag)
            (method) either 'semivariogram', or 'covariogram'
    Output: (cv)   <2xN> NumPy array of lags and variogram values
    '''
    # calculate the pairwise distances
    pwdist = utilities.pairwise(data)
    # create a list of lists of indices of points having the ~same lag
    index = [lagindices(pwdist, lag, tol) for lag in lags]
    # remove empty "lag" sets, prevents zero division error in [co|semi]variance()
    index = list(filter(lambda x: len(x) > 0, index))
    # calculate the variogram at different lags given some tolerance
    if method in ['semivariogram', 'semi', 'sv', 's']:
        v = [semivariance(data, indices) for indices in index]
    elif method in ['covariogram', 'cov', 'co', 'cv', 'c']:
        v = [covariance(data, indices) for indices in index]
    # bundle the semivariogram values with their lags
    return np.c_[lags, v].T
