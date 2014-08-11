geostatsmodels
==============

I am implementing ideas from Clayton V. Deutsch's book Geostatistical Reservoir Modeling in order to better understand geostatistics. I am working through things first in Python in order to prototype things quickly and make sure I understand them. Later I will try to optimize the code for speed.

I used http://people.ku.edu/~gbohling/cpe940 as a starting point to understand kriging. There is also an R package for geostatistics named gstat, but I have not used it yet.

This software is licensed under the [MIT License](http://opensource.org/licenses/MIT).

Installation and Dependencies
-----------------------------

This package depends on:
 * numpy, the fundamental package for scientific computing with Python. <a href="http://www.numpy.org/">http://www.numpy.org/</a>  
 * matplotlib, a Python 2D plotting library which produces publication quality figures. <a href="<a href="http://matplotlib.org/index.html">http://matplotlib.org/index.html</a>
 * scipy, a Python library which provides many user-friendly and efficient numerical routines such as routines for numerical integration and optimization. <a href="<a href="http://www.scipy.org/scipylib/index.html">http://www.scipy.org/scipylib/index.html</a>
 * pandas, a Python library which implements excellent tools for data analysis and modeling. 

The best/easiest way to install these packages is to use one of the Python distributions described here:

<a href="http://www.scipy.org/install.html">http://www.scipy.org/install.html</a>

Anaconda has been successfully tested with geostatsmodels.

Most of those distributions should include <em>pip</em>, a command line tool for installing and 
managing Python packages.  You can use pip to install geostatsmodels itself.  
 
You may need to open a new terminal window to ensure that the newly installed versions of python and pip
are in your path.

To install geostatsmodels:

pip install git+git://github.com/cjohnson318/geostatsmodels.git

Uninstalling and Updating
-------------------------

To uninstall:

pip uninstall geostatsmodels

To update:

pip install -U git+git://github.com/cjohnson318/geostatsmodels.git

Usage
------
A notebook exploring some of the functionality of geostatsmodels is included in this repository.  

[Variogram Analysis](http://nbviewer.ipython.org/github/cjohnson318/geostatsmodels/blob/master/notebooks/VariogramAnalysis.ipynb)

[Kriging Example](http://nbviewer.ipython.org/github/cjohnson318/geostatsmodels/blob/master/notebooks/KrigingExample.ipynb)

More of these will follow.

