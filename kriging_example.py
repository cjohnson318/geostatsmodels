#!/usr/bin/env python

import geostatsmodels.utilities as u
import geostatsmodels.kriging as k

# grab the data from cluster.dat
d = u.readGeoEAS("cluster.dat")

# use only the first three columns
d = d[:,:3]

# hs is an array of distances to use in semivariogram modeling
hs = np.linspace( 0, 50, 10 )

# bw is the bandwidth, the plus/minus distance around
# the distances listed in hs to use in semivariogram modeling
bw = 5

# for this example, z should be 3.9225013671953293
z = k.krige( d, k.spherical, hs, bw, [25.,25.], 4 )
