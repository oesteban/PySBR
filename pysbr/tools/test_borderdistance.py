#!/usr/bin/env python
# coding: utf-8
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

"""

test_borderdistance.py - 

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
from math import sqrt

# add the setup build path before loading the module 
from testtools import append_buildpath
append_buildpath()
import csbr


class TestSimple(unittest.TestCase):
    def runtest_shapestats(self, o):
        (volume, gravcenter, idx) = self.create_test_shape('F')
        spot = (1,2,3)
        spotarray = np.array(spot, dtype=np.double)
        expect = (sqrt(48.0),7.3023054915, 3.6075617478)
        
        (gravcenterdist, distmean, distvar) = csbr.shapestats(volume, spot)
        self.assertAlmostEqual(gravcenterdist, expect[0])
        self.assertAlmostEqual(distmean, expect[1])
        self.assertAlmostEqual(distvar, expect[2])
        
        
    def create_test_shape(self, o):
        data =  np.zeros((12,13,14), dtype=np.bool, order=o)
        idx = [np.array([2,3,6]), np.array([8,9,6]), np.array([5,6,9])]

        for i in idx: 
            data[i[0],i[1],i[2]] = True
        gravcenter = np.array([5, 6, 7])
        return (data, gravcenter, idx)

    
    def test_order_fortran(self):
        self.runtest_shapestats('F')

            
    def test_order_C(self):
        self.runtest_shapestats('C')
    
    
        
if __name__ == '__main__':
	unittest.main()

    
