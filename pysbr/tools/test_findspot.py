#!/usr/bin/env python
# coding: utf-8
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

"""

test_findspot.py -

Copyright (c) 2013, code@oscaresteban.es (Oscar Esteban), and 
                    gw.fossdev@gmail.com (Gert Wollny)
                    with Biomedical Image Technology, UPM (BIT-UPM)
All rights reserved.
This file is part of PySBR.

"""

__author__ = "Gert Wollny"
__copyright__ = "Copyright 2013, Biomedical Image Technologies (BIT), \
                 Universidad Politécnica de Madrid"
__credits__ = ["Oscar Esteban", "Gert Wollny"]
__license__ = "FreeBSD"
__version__ = "0.1"
__maintainer__ = "Gert Wollny"
__email__ = "gw.fossdev@gmail.com"
__status__ = "Prototype"

import numpy as np
import unittest

# add the setup build path before loading the module 
from testtools import append_buildpath
append_buildpath()
import csbr

class TestSimple(unittest.TestCase):
    def runtest_findspot(self, o, expect, expect_volume, mean, bridge, volume=None):
        data =  np.array([
            [[   0,   0,   0,   0,   0,   0,   0], 
			 [   0,   0,   9,   0,   0,   0,   0],
			 [   0,   9,   7,   1,   0,   0,   0],
			 [   0,   8,   0,   1,   0,   0,   0],
			 [   0,   1,   1,   1,   1,   0,   0], 
			 [   0,   0,   0,   1,   0,   0,   0],
			 [   0,   0,   0,   1,   0,   0,   0],
			 [   0,   0,   0,   0,   0,   0,   0]],

			[[   0,   0,   0,   0,   0,   0,   0], 
			 [   0,  99,  99,  50,  90, 120,   0],
			 [   0, 100,  97,  50,  90,  89,   0],
			 [   0,  98,   0,   1,   0,  79,   0],
			 [   0,   1,   1,   1,   1,   1,   0], 
			 [   0,   0,   0,   1,   0,   0,   0],
			 [   0,   0,   0,   1,   0,   0,   0],
			 [   0,   0,   0,   0,   0,   0,   0]],
				  
			[[   0,   0,   0,   0,   0,   0,   0], 
			 [   0,  99,  10,   1,   0, 111,   0],
			 [   0,  99,  10,   1,   0,  89,   0],
			 [   0,  10,   1,   1,   1,  10,   0],
			 [   0,   1,   1,   1,   1,  21,   0], 
			 [   0,   0,   1,   1,   1,   0,   0],
			 [   0,   0,   0,   1,   0,   0,   0],
			 [   0,   0,   0,   0,   0,   0,   0]],

			[[   0,   0,   0,   0,   0,   0,   0], 
			 [   0,  89,   0,   1,   0, 100,   0],
			 [   0,  87,   0,   1,   0,  90,   0],
			 [   0,   0,   1,   0,   1,   0,   0],
			 [   0,   0,   1,   0,   1,   1,   0], 
			 [   0,   0,   1,   1,   0,   1,   0],
			 [   0,   0,   0,   0,   1,   1,   0],
			 [   0,   0,   0,   0,   0,   0,   0]]], dtype=np.int32, order=o)


    
        if volume is None:
            (result, nspots, hotpoints, volumes) = csbr.findspots3d(data, 0.7, mean, bridge)
        else:
            (result, nspots, hotpoints, volumes) = csbr.findspots3d(data, 0.7, mean, bridge, maxvolume=volume)

        
        self.assertEqual(nspots, 2)

        print (hotpoints)
        self.assertEqual(hotpoints[0][0],1)
        self.assertEqual(hotpoints[0][1],1)
        self.assertEqual(hotpoints[0][2],5)

        self.assertEqual(hotpoints[1][0],1)
        self.assertEqual(hotpoints[1][1],2)
        self.assertEqual(hotpoints[1][2],1)


        self.assertEqual(volumes[0], expect_volume[0])
        self.assertEqual(volumes[1], expect_volume[1])
        
        delta = np.array_equal(result, expect)

    
        self.assertEqual(delta, True)


    def runtest_findspot_mirrored(self, o, expect, expect_volume, mean, bridge, volume=None):
        data =  np.array([
            [[   0,   0,   0,   0,   0,   0,   0], 
			 [   0,   0,   9,   0,   0,   0,   0],
			 [   0,   9,   7,   1,   0,   0,   0],
			 [   0,   8,   0,   1,   0,   0,   0],
			 [   0,   1,   1,   1,   1,   0,   0], 
			 [   0,   0,   0,   1,   0,   0,   0],
			 [   0,   0,   0,   1,   0,   0,   0],
			 [   0,   0,   0,   0,   0,   0,   0]],

			[[   0,   0,   0,   0,   0,   0,   0], 
			 [   0,  99,  99,  50,  90, 100,   0],
			 [   0, 120,  97,  50,  90,  89,   0],
			 [   0,  98,   0,   1,   0,  79,   0],
			 [   0,   1,   1,   1,   1,   1,   0], 
			 [   0,   0,   0,   1,   0,   0,   0],
			 [   0,   0,   0,   1,   0,   0,   0],
			 [   0,   0,   0,   0,   0,   0,   0]],
				  
			[[   0,   0,   0,   0,   0,   0,   0], 
			 [   0,  99,  10,   1,   0,  99,   0],
			 [   0,  99,  10,   1,   0,  89,   0],
			 [   0,  10,   1,   1,   1,  10,   0],
			 [   0,   1,   1,   1,   1,  21,   0], 
			 [   0,   0,   1,   1,   1,   0,   0],
			 [   0,   0,   0,   1,   0,   0,   0],
			 [   0,   0,   0,   0,   0,   0,   0]],

			[[   0,   0,   0,   0,   0,   0,   0], 
			 [   0,  89,   0,   1,   0,  99,   0],
			 [   0,  87,   0,   1,   0,  90,   0],
			 [   0,   0,   1,   0,   1,   0,   0],
			 [   0,   0,   1,   0,   1,   1,   0], 
			 [   0,   0,   1,   1,   0,   1,   0],
			 [   0,   0,   0,   0,   1,   1,   0],
			 [   0,   0,   0,   0,   0,   0,   0]]], dtype=np.int32, order=o)


    
        if volume is None:
            (result, nspots, hotpoints, volumes) = csbr.findspots3d(data, 0.7, mean, bridge)
        else:
            (result, nspots, hotpoints, volumes) = csbr.findspots3d(data, 0.7, mean, bridge, maxvolume=volume)

        
        self.assertEqual(nspots, 2)


        self.assertEqual(hotpoints[0][0],1)
        self.assertEqual(hotpoints[0][1],1)
        self.assertEqual(hotpoints[0][2],5)

        self.assertEqual(hotpoints[1][0],1)
        self.assertEqual(hotpoints[1][1],2)
        self.assertEqual(hotpoints[1][2],1)


        self.assertEqual(volumes[0], expect_volume[0])
        self.assertEqual(volumes[1], expect_volume[1])
        
        delta = np.array_equal(result, expect)

        self.assertEqual(delta, True)

        
    def test_simple_limit_bridge(self):
        expect = np.array([
            [[   0,   0,   0,   0,   0,   0,   0], 
			 [   0,   0,   0,   0,   0,   0,   0],
			 [   0,   0,   0,   0,   0,   0,   0],
			 [   0,   0,   0,   0,   0,   0,   0],
			 [   0,   0,   0,   0,   0,   0,   0], 
			 [   0,   0,   0,   0,   0,   0,   0],
			 [   0,   0,   0,   0,   0,   0,   0],
			 [   0,   0,   0,   0,   0,   0,   0]],

            [[   0,   0,   0,   0,   0,   0,   0], 
			 [   0,   2,   2,   0,   1,   1,   0],
			 [   0,   2,   2,   0,   1,   0,   0],
			 [   0,   2,   0,   0,   0,   0,   0],
			 [   0,   0,   0,   0,   0,   0,   0], 
			 [   0,   0,   0,   0,   0,   0,   0],
			 [   0,   0,   0,   0,   0,   0,   0],
			 [   0,   0,   0,   0,   0,   0,   0]],
				  
			[[   0,   0,   0,   0,   0,   0,   0], 
			 [   0,   2,   0,   0,   0,   1,   0],
			 [   0,   2,   0,   0,   0,   0,   0],
			 [   0,   0,   0,   0,   0,   0,   0],
			 [   0,   0,   0,   0,   0,   0,   0], 
			 [   0,   0,   0,   0,   0,   0,   0],
			 [   0,   0,   0,   0,   0,   0,   0],
			 [   0,   0,   0,   0,   0,   0,   0]],

			[[   0,   0,   0,   0,   0,   0,   0], 
			 [   0,   0,   0,   0,   0,   1,   0],
			 [   0,   0,   0,   0,   0,   1,   0],
			 [   0,   0,   0,   0,   0,   0,   0],
			 [   0,   0,   0,   0,   0,   0,   0], 
			 [   0,   0,   0,   0,   0,   0,   0],
			 [   0,   0,   0,   0,   0,   0,   0],
			 [   0,   0,   0,   0,   0,   0,   0]]], dtype=np.short)

        for order in ['C', 'F']: 
            self.runtest_findspot(order, expect, (6,7), 0.7, 9.0/5.0 )
            self.runtest_findspot_mirrored(order, expect, (6,7), 0.7, 9.0/5.0)
		
    def test_simple_limit_volume(self):
        expect = np.array([
            [[   0,   0,   0,   0,   0,   0,   0], 
			 [   0,   0,   0,   0,   0,   0,   0],
			 [   0,   0,   0,   0,   0,   0,   0],
			 [   0,   0,   0,   0,   0,   0,   0],
			 [   0,   0,   0,   0,   0,   0,   0], 
			 [   0,   0,   0,   0,   0,   0,   0],
			 [   0,   0,   0,   0,   0,   0,   0],
			 [   0,   0,   0,   0,   0,   0,   0]],

            [[   0,   0,   0,   0,   0,   0,   0], 
			 [   0,   2,   2,   0,   1,   1,   0],
			 [   0,   2,   0,   0,   1,   0,   0],
			 [   0,   2,   0,   0,   0,   0,   0],
			 [   0,   0,   0,   0,   0,   0,   0], 
			 [   0,   0,   0,   0,   0,   0,   0],
			 [   0,   0,   0,   0,   0,   0,   0],
			 [   0,   0,   0,   0,   0,   0,   0]],
				  
			[[   0,   0,   0,   0,   0,   0,   0], 
			 [   0,   2,   0,   0,   0,   1,   0],
			 [   0,   2,   0,   0,   0,   0,   0],
			 [   0,   0,   0,   0,   0,   0,   0],
			 [   0,   0,   0,   0,   0,   0,   0], 
			 [   0,   0,   0,   0,   0,   0,   0],
			 [   0,   0,   0,   0,   0,   0,   0],
			 [   0,   0,   0,   0,   0,   0,   0]],

			[[   0,   0,   0,   0,   0,   0,   0], 
			 [   0,   0,   0,   0,   0,   1,   0],
			 [   0,   0,   0,   0,   0,   1,   0],
			 [   0,   0,   0,   0,   0,   0,   0],
			 [   0,   0,   0,   0,   0,   0,   0], 
			 [   0,   0,   0,   0,   0,   0,   0],
			 [   0,   0,   0,   0,   0,   0,   0],
			 [   0,   0,   0,   0,   0,   0,   0]]], dtype=np.short)

        for order in ['C', 'F']: 
            self.runtest_findspot(order, expect, (6,6), 0.7, 0.0, 6)
            self.runtest_findspot_mirrored(order, expect, (6,6), 0.7, 0.0, 6)            






            
    def test_volume(self):
        expect = np.array([
            [[   0,   0,   0,   0,   0,   0,   0], 
			 [   0,   0,   0,   0,   0,   0,   0],
			 [   0,   0,   0,   0,   0,   0,   0],
			 [   0,   0,   0,   0,   0,   0,   0],
			 [   0,   0,   0,   0,   0,   0,   0], 
			 [   0,   0,   0,   0,   0,   0,   0],
			 [   0,   0,   0,   0,   0,   0,   0],
			 [   0,   0,   0,   0,   0,   0,   0]],

            [[   0,   0,   0,   0,   0,   0,   0], 
			 [   0,   2,   2,   0,   1,   1,   0],
			 [   0,   2,   2,   0,   1,   1,   0],
			 [   0,   2,   0,   0,   0,   0,   0],
			 [   0,   0,   0,   0,   0,   0,   0], 
			 [   0,   0,   0,   0,   0,   0,   0],
			 [   0,   0,   0,   0,   0,   0,   0],
			 [   0,   0,   0,   0,   0,   0,   0]],
				  
			[[   0,   0,   0,   0,   0,   0,   0], 
			 [   0,   2,   0,   0,   0,   1,   0],
			 [   0,   2,   0,   0,   0,   1,   0],
			 [   0,   0,   0,   0,   0,   0,   0],
			 [   0,   0,   0,   0,   0,   0,   0], 
			 [   0,   0,   0,   0,   0,   0,   0],
			 [   0,   0,   0,   0,   0,   0,   0],
			 [   0,   0,   0,   0,   0,   0,   0]],

			[[   0,   0,   0,   0,   0,   0,   0], 
			 [   0,   2,   0,   0,   0,   1,   0],
			 [   0,   2,   0,   0,   0,   1,   0],
			 [   0,   0,   0,   0,   0,   0,   0],
			 [   0,   0,   0,   0,   0,   0,   0], 
			 [   0,   0,   0,   0,   0,   0,   0],
			 [   0,   0,   0,   0,   0,   0,   0],
			 [   0,   0,   0,   0,   0,   0,   0]]], dtype=np.short)

        for order in ['C', 'F']: 
            self.runtest_findspot(order, expect, (8,9), 2.2, 0.0)
            self.runtest_findspot_mirrored(order, expect, (8,9), 2.2, 0.0)

if __name__ == '__main__':
	unittest.main()
