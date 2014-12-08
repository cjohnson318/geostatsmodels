# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# This notebook is the reproduction of an exercise found at http://people.ku.edu/~gbohling/cpe940/Kriging.pdf

# <codecell>

from geostatsmodels import utilities, variograms, model, kriging, geoplot
import pandas

# <markdowncell>

# We'll read the data from `ZoneA.dat`.

# <codecell>

z = utilities.readGeoEAS('../data/ZoneA.dat')

# <markdowncell>

# We want the first, second and fourth columns of the data set, representing the x and y spatial coordinates, and the porosity.

# <codecell>

P = z[:,[0,1,3]]

# <markdowncell>

# We'll be interested in determining the porosity at a point (2000,4700).

# <codecell>

pt = [ 2000, 4700 ]

# <markdowncell>

# We can plot our region of interest as follows:

# <codecell>

scatter( P[:,0], P[:,1], c=P[:,2], cmap=geoplot.YPcmap )
title('Zone A Subset % Porosity')
colorbar()
xmin, xmax = 0, 4250
ymin, ymax = 3200, 6250
xlim(xmin,xmax)
ylim(ymin,ymax)
for i in range( len( P[:,2] ) ):
    x, y, por = P[i]
    if( x < xmax )&( y > ymin )&( y < ymax ):
        text( x+100, y, '{:4.2f}'.format( por ) ) 
scatter( pt[0], pt[1], marker='x', c='k' )
text( pt[0]+100 , pt[1], '?')
xlabel('Easting (m)')
ylabel('Northing (m)') ;

# <markdowncell>

# We can determine the parameters for our model by looking at the semivariogram and trying to determine the appropriate range and sill.

# <codecell>

tolerance = 250
lags = np.arange( tolerance, 10000, tolerance*2 )
sill = np.var( P[:,2] )

# <markdowncell>

# The semivariogram plotting function, `svplot()`, plots sill as a dashed line, and the empirical semivariogram as determined from the data. It optionally plots a semivariance model.

# <codecell>

geoplot.semivariogram( P, lags, tolerance )

# <markdowncell>

# We can pass a model to this function using the optional `model` argument and see it plotted in red.

# <codecell>

svm = model.semivariance( model.spherical, ( 4000, sill ) )
geoplot.semivariogram( P, lags, tolerance, model=svm )

# <markdowncell>

# The covariance modeling function function will return a spherical covariance model that takes a distance as input, and returns an covariance estimate. We've used the global variance of the porosity in `ZoneA.dat` as the sill.

# <codecell>

covfct = model.covariance( model.spherical, ( 4000, sill ) )

# <markdowncell>

# We can then krige the data, using the covariance model, the point we are interested in, (2000,47000), and `N=6` signifying that we only want to use the six nearest points. The output of the simple and ordinary kriging functions below is the krigin estimate, and the standard deviation of the kriging estimate.

# <codecell>

kriging.simple( P, covfct, pt, N=6 )

# <codecell>

kriging.ordinary( P, covfct, pt, N=6 )

# <codecell>

est, kstd = kriging.krige( P, covfct, [[2000,4700],[2100,4700],[2000,4800],[2100,4800]], 'simple', N=6 )

# <codecell>

est

# <codecell>

kstd

