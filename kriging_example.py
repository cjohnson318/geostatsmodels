#!/usr/bin/env python

import numpy as np
import geostatsmodels.utilities as u
import geostatsmodels.models as m
import geostatsmodels.kriging as k

# load data from the GeoEAS data format
data = u.readGeoEAS("data/cluster.dat")

# define a +/- tolerance for the lags
tol = 5

# create a NumPy array of lags
lags = np.arange( 10, 50, tol*2 )

# create a function to model the variance
covfct = m.covmodel( data, m.gaussian, lags, tol )

# specify a point to estimate
z = [ 25.0, 25.0 ]

# pass that function to the kriging unction
kv = k.krige( data, covfct, lags, tol, z )

print kv
