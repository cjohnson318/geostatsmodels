#!/usr/bin/env python

import utilities as u
import kriging as k
import simulation as s
from pylab import *

# read cluster.dat
d = u.readGeoEAS("cluster.dat")

# take the first three columns
d = d[:,:3]

# define the distances and bandwidth 
# for the semivariogram modeling
hs, bw = np.linspace( 0, 50, 10 ), 5

# if the data is not normally distributed,
# perfrom a z-score transformation
d, inv = u.to_norm( d )

# perform sequential Gaussian simulation using
# a spherical model, on a 5x5 grid
m = s.sgs( d, k.spherical, hs, bw, 5, 5 )

# use the [:,::-1].T to a) reverse the order of the columns
# and then b) transpose the data, this takes it from Python
# conventions, back to the way we normally think of spatial data
print m[:,::-1].T

'''
Printing m, you should see something like the following output.
Notice that these values look like they have a mean of zero, and a
standard deviation of one. This data has not been transformed back
to the original distribution.

array([[ 0.03499825, -0.56577758, -0.53807829,  0.45723896,  0.33959872],
       [-0.42029471, -0.27813954,  0.07044366,  0.39437271, -0.07625142],
       [-1.12918974, -0.36841119,  0.21116697, -1.11767303, -0.48140453],
       [ 0.39862306, -0.39414717, -0.95340368, -0.46121427, -0.08815428],
       [ 0.50258498, -0.91996145, -0.91149543, -0.61183492, -0.88997225]])
'''

# assuming you have matplotlib and  pylab installed, 
# you can visulize the untransformed data
matshow( m[:,::-1].T, cmap=u.YPcmap )

# perform the back-transformation
n = u.from_norm( m, inv )

print n[:,::-1].T

'''
The back-transformed data should look something like this:

array([[ 2.32771604,  0.89024448,  0.9233661 ,  3.63025147,  3.28306658],
       [ 1.01198912,  1.35326734,  2.33965581,  3.53342939,  1.94217347],
       [ 0.26234492,  1.10879666,  2.80242068,  0.26918933,  0.96348094],
       [ 3.54440162,  1.06344535,  0.33826991,  0.99125142,  1.89579485],
       [ 4.25739527,  0.36063006,  0.36684815,  0.82690577,  0.3828731 ]])
'''
