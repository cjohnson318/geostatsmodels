#!/usr/bin/env python

import unittest
from geostatsmodels import metrics as m
import numpy as np
import pandas as pd

meuse = pd.read_csv("data/meuse.csv", sep=',', header=0)
eps = 0.0001

# https://rstudio-pubs-static.s3.amazonaws.com/79757_3462f2bba8d745378545f0ea5bc38ee1.html

class metrics_cases( unittest.TestCase ):
	'''Tests for metrics.py'''
	
	def test_moran_i_meause_lead( self ):
		data = meuse[["x", "y", "lead"]].to_numpy()
		'''
		Does moran_i() return a result within
		eps=0.0001 of a benchmarked value?
		'''
		i,e = m.moran_i( data )
		'''
		R code used to retrieve the answer value:
			library(ape)
			library(sp)
			data(meuse)
			coordinates(meuse) <- ~x+y
			proj4string(meuse) <- CRS("+init=epsg:28992")
			w <- 1/as.matrix(dist(coordinates(meuse)))
			diag(w) <- 0
			Moran.I(meuse$lead, w)
		'''
		ans = 0.1007573
		print("I = ",i, ", answer = ", ans, ", expected = ", e)
		self.assertTrue( abs( ans - i ) < eps )
		
	def test_moran_i_meause_copper( self ):
		
		data = meuse[["x", "y", "copper"]].to_numpy()
		
		'''
		Does moran_i() return a result within
		eps=0.0001 of a benchmarked value?
		'''
		i,e = m.moran_i( data )
		'''
		R code used to retrieve the answer value:
			library(ape)
			library(sp)
			data(meuse)
			coordinates(meuse) <- ~x+y
			proj4string(meuse) <- CRS("+init=epsg:28992")
			w <- 1/as.matrix(dist(coordinates(meuse)))
			diag(w) <- 0
			Moran.I(meuse$copper, w)
		'''
		ans = 0.09301739
		print("I = ",i, ", answer = ", ans, ", expected = ", e)
		self.assertTrue( abs( ans - i ) < eps )

if __name__ == '__main__':
    unittest.main()
