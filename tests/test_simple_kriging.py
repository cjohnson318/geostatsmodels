#!/usr/bin/env python

import unittest
from geostatsmodels import kriging as k
import numpy as np

class KrigingTestCase( unittest.TestCase ):
	'''Tests for kriging.py'''
	
	def test_simple_kriging_spherical_model( self ):
		'''
		Does krige() return a result within
		eps=0.0001  of a benchmarked value?
		'''
		data = np.array([[0,0,1],[1,0,3],[0,1,3],[1,1,4]])
		hs = np.array([0.5,1.0,1.5,2.0])
		bw = 0.5
		u = [0.5,0.5]
		N = 2
		kv = k.krige( data, k.spherical, hs, bw, u, N )
		ans = 2.29891949337
		eps = 0.0001
		self.assertTrue( abs( ans - kv ) < eps ) 

if __name__ == '__main__':
    unittest.main()
