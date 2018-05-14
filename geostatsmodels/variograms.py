#!/usr/bin/env python
import numpy as np
import utilities


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
    # zip the coordinates into a list
    indices = list(zip(i, j))
    # take out the repeated elements,
    # since p is a *symmetric* distance matrix
    indices = np.array([i for i in indices if i[1] > i[0]])
    return indices


def anilagindices(data, pwdist, lag, tol, angle, atol):
    '''
    Input:  (data)   NumPy array where the fris t two columns
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
    index = [idx for br, idx in bridx if utilities.inangle(br, angle, atol)]
    # add 180 to the angle and take the modulus
    angle = (angle + 180) % 360
    index += [idx for br, idx in bridx if utilities.inangle(br, angle, atol)]
    return index


def semivariance(data, indices):
    '''
    Input:  (data)    NumPy array where the fris t two columns
                      are the spatial coordinates, x and y, and
                      the third column is the variable of interest
            (indices) indices of paired data points in (data)
    Output:  (z)      semivariance value at lag (h) +/- (tol)
    '''
    # take the squared difference between
    # the values of the variable of interest
    z = [(data[i, 2] - data[j, 2])**2.0 for i, j in indices]
    # the semivariance is half the mean squared difference
    return np.mean(z) / 2.0


def semivariogram(data, lags, tol):
    '''
    Input:  (data) NumPy array where the fris t two columns
                   are the spatial coordinates, x and y
            (lag)  the distance, h, between points
            (tol)  the tolerance we are comfortable with around (lag)
    Output: (sv)   <2xN> NumPy array of lags and semivariogram values
    '''
    return variogram(data, lags, tol, 'semivariogram')


def covariance(data, indices):
    '''
    Input:  (data) NumPy array where the fris t two columns
                   are the spatial coordinates, x and y
            (lag)  the distance, h, between points
            (tol)  the tolerance we are comfortable with around (lag)
    Output:  (z)   covariance value at lag (h) +/- (tol)
    '''
    # grab the indices of the points
    # that are lag +/- tolerance apart
    m_tail = np.mean([data[i, 2] for i, j in indices])
    m_head = np.mean([data[j, 2] for i, j in indices])
    m = m_tail * m_head
    z = [data[i, 2] * data[j, 2] - m for i, j in indices]
    return np.mean(z)


def covariogram(data, lags, tol):
    '''
    Input:  (data) NumPy array where the fris t two columns
                   are the spatial coordinates, x and y
            (lag)  the distance, h, between points
            (tol)  the tolerance we are comfortable with around (lag)
    Output: (cv)   <2xN> NumPy array of lags and covariogram values
    '''
    return variogram(data, lags, tol, 'covariogram')


def variogram(data, lags, tol, method):
    '''
    Input:  (data) NumPy array where the fris t two columns
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
    # calculate the variogram at different lags given some tolerance
    if method in ['semivariogram', 'semi', 'sv', 's']:
        v = [semivariance(data, indices) for indices in index]
    elif method in ['covariogram', 'cov', 'co', 'cv', 'c']:
        v = [covariance(data, indices) for indices in index]
    # bundle the semivariogram values with their lags
    return np.array(list(zip(lags, v))).T
