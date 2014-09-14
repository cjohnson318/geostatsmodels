# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

%matplotlib inline

# <markdowncell>

# This notebook is an effort to replicate the lessons found here:
# http://people.ku.edu/~gbohling/cpe940/Variograms.pdf
# 
# We'll do all of our imports here at the top.
# 
# *Note: This is a work in progress!*

# <codecell>

import numpy as np
from pandas import DataFrame, Series
from geostatsmodels import utilities, kriging, variograms, model, geoplot
import matplotlib.pyplot as plt
from scipy.stats import norm
import urllib2 
import os.path
import zipfile
import StringIO

# <markdowncell>

# We're going to fetch the data file we need for this exercise from the following URL:
# 
# http://people.ku.edu/~gbohling/geostats/WGTutorial.zip
# 
# Subsequent runs of this Notebook should use a local copy, saved in the current directory.

# <codecell>

clusterfile = 'ZoneA.dat'
if not os.path.isfile(clusterfile):
    fh = urllib2.urlopen('http://people.ku.edu/~gbohling/geostats/WGTutorial.zip')
    data = fh.read()
    fobj = StringIO.StringIO(data)
    myzip = zipfile.ZipFile(fobj,'r')
    myzip.extract(clusterfile)
    fobj.close()
    fh.close()
z = open(clusterfile,'r' ).readlines()
z = [ i.strip().split() for i in z[10:] ]
z = np.array( z, dtype=np.float )
z = DataFrame( z, columns=['x','y','thk','por','perm','lperm','lpermp','lpermr'] )
P = np.array( z[['x','y','por']] )

# <markdowncell>

# Let's make a plot of the data, so we know what we're dealing with.

# <codecell>

fig, ax = plt.subplots()
fig.set_size_inches(8,8)
cmap = geoplot.YPcmap
ax.scatter( z.x/1000, z.y/1000, c=z.por, s=64,cmap=cmap)
ax.set_aspect(1)
plt.xlim(-2,22)
plt.ylim(-2,17.5)

plt.xlabel('Easting [m]')
plt.ylabel('Northing [m]')
th=plt.title('Porosity %')

# <markdowncell>

# Let's verify that our data is distributed normally.

# <codecell>

hrange = (12,17.2)
mu, std = norm.fit(z.por)
ahist=plt.hist(z.por, bins=7, normed=True, alpha=0.6, color='c',range=hrange)

xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 100)
p = norm.pdf(x, mu, std)
plt.plot(x, p, 'k', linewidth=2)
title = "Fit results: mu = %.2f,  std = %.2f" % (mu, std)
th=plt.title(title)
xh=plt.xlabel('Porosity (%)')
yh=plt.ylabel('Density')
xl=plt.xlim(11.5,17.5)
yl=plt.ylim(-0.02,0.45)

# <codecell>

import scipy.stats as stats
qqdata = stats.probplot(z.por, dist="norm",plot=plt,fit=False)
xh=plt.xlabel('Standard Normal Quantiles')
yh=plt.ylabel('Sorted Porosity Values')
fig=plt.gcf()
fig.set_size_inches(8,8)
th=plt.title('')

# <markdowncell>

# What is the optimal "lag" distance between points?  Use the utilities scattergram() function to help determine that distance.

# <codecell>

pw = utilities.pairwise(P)
geoplot.hscattergram(P,pw,1000,500)
geoplot.hscattergram(P,pw,2000,500)
geoplot.hscattergram(P,pw,3000,500)

# <markdowncell>

# Here, we plot the semivariogram and overlay a horizontal line for the sill, $c$.

# <codecell>

tolerance = 250
lags = np.arange( tolerance, 10000, tolerance*2 )
sill = np.var(P[:,2])

geoplot.semivariogram( P, lags, tolerance )

# <markdowncell>

# Looking at the figure above, we can say that the semivariogram levels off around 4000, so we can set the range, $a$, to that value and model the covariance function.

# <codecell>

svm = model.semivariance( model.spherical, [ 4000, sill ] )
geoplot.semivariogram( P, lags, tolerance, model=svm )

# <markdowncell>

# We can visualize the distribution of the lagged distances with the `laghistogram()` function.

# <codecell>

geoplot.laghistogram( P, pw, lags, tolerance )

# <markdowncell>

# If we want to perform anisotropic kriging, we can visualize the distribution of the anisotropic lags using the `anisotropiclags()` function. Note that we use the bearing, which is measured in degrees, clockwise from North.

# <codecell>

geoplot.anisotropiclags( P, pw, lag=2000, tol=250, angle=45, atol=15 )

# <codecell>

geoplot.anisotropiclags( P, pw, lag=2000, tol=250, angle=135, atol=15 )

# <codecell>

geoplot.polaranisotropy( P, pw, lags, tolerance, nsectors=12 )

# <codecell>


