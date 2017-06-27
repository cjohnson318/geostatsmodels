geostatsmodels
==============

This is an extension of [cjohnson318](http://www.github.com/cjohnson318)'s implementation of concepts from Clayton V. Deutsch and Michael Pyrcz's book _Geostatistical Reservoir Modeling_, currently available via [Oxford University Press](https://global.oup.com/academic/product/geostatistical-reservoir-modeling-9780199731442?cc=us&lang=en&). 

This software is licensed under the [MIT License](http://opensource.org/licenses/MIT).

Installation and Dependencies
-----------------------------

This package depends on:
 * `numpy`, the fundamental package for scientific computing with Python. 
 * `matplotlib`, a Python 2D plotting library which produces publication quality figures. 
 * scipy, a Python library which provides many user-friendly and efficient numerical routines such as routines for numerical integration and optimization. 
 * `pandas`, a Python library which implements excellent tools for data analysis and modeling. 

The best/easiest way to install these packages is to use one of the Python distributions described [here](http://www.scipy.org/install.html). Anaconda has been successfully tested with geostatsmodels.

Most of those distributions should include `pip` or `conda`, which are command line tools for installing and 
managing Python packages.  You can use `pip` to install geostatsmodels itself.  
 
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
Some notebooks exploring the functionality of `geostatsmodels` are included below.  

[Variogram Analysis](http://nbviewer.ipython.org/github/cjohnson318/geostatsmodels/blob/master/notebooks/VariogramAnalysis.ipynb)

[Kriging Example](http://nbviewer.ipython.org/github/cjohnson318/geostatsmodels/blob/master/notebooks/KrigingExample.ipynb)

More to come!

