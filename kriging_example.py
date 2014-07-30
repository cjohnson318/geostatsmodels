#!/usr/bin/env python

import numpy as np

import geostatsmodels.utilities as u
import geostatsmodels.kriging as k

# grab the data from cluster.dat
d = u.readGeoEAS("cluster.dat")

# use only the first three columns
d = d[:,:3]

# lags is an array of distances to use in semivariogram modeling
lags = np.linspace( 10, 50, 10 )

# tolerance is a plus/minus distance around the lags
tol = 5

z = k.krige( d, k.spherical, lags, tol, [25.,25.] )
