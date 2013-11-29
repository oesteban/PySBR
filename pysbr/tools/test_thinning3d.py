#!/usr/bin/env python
# coding: utf-8
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

"""

test_thinning3d.py -

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

from testtools import append_buildpath

append_buildpath()
import csbr

class TestSimple(unittest.TestCase):
	def runtest_thinning(self, o, expect):
		data =  np.array([
			[[0,0,0,0,0,0,0], 
			 [0,0,0,1,0,0,0],
			 [0,0,0,1,0,0,0],
			 [0,0,0,1,0,0,0],
			 [0,1,1,1,1,1,0], 
			 [0,0,0,1,0,0,0],
			 [0,0,0,1,0,0,0],
			 [0,0,0,0,0,0,0]],
				  
			[[0,0,0,0,0,0,0], 
			 [0,0,0,1,0,0,0],
			 [0,0,0,1,0,0,0],
			 [0,0,1,1,1,0,0],
			 [0,1,1,1,1,1,0], 
			 [0,0,1,1,1,0,0],
			 [0,0,0,1,0,0,0],
			 [0,0,0,0,0,0,0]],

			[[0,0,0,0,0,0,0], 
			 [0,1,1,1,0,0,0],
			 [0,1,0,1,1,0,0],
			 [0,1,1,1,1,1,0],
			 [0,1,1,1,1,1,0], 
			 [0,0,1,1,0,0,0],
			 [0,0,0,1,0,0,0],
			 [0,0,0,0,0,0,0]],

			[[0,0,0,0,0,0,0], 
			 [0,0,0,1,0,0,0],
			 [0,0,0,1,0,0,0],
			 [0,0,1,0,1,0,0],
			 [0,0,1,0,1,1,0], 
			 [0,0,1,1,0,1,0],
			 [0,0,0,0,1,1,0],
			 [0,0,0,0,0,0,0]]], dtype=np.bool, order=o)


		csbr.thinning3d(data)
		delta = np.logical_not(np.logical_xor(data, expect))
		
		self.assertEqual(np.all(delta), True)

		
		
	def test_simple(self):
		expect = np.array([
			[[0,0,0,0,0,0,0], 
			 [0,0,0,0,0,0,0],
			 [0,0,0,0,0,0,0],
			 [0,0,0,0,0,0,0],
			 [0,0,0,0,0,0,0], 
			 [0,0,0,0,0,0,0],
			 [0,0,0,0,0,0,0],
			 [0,0,0,0,0,0,0]],

			[[0,0,0,0,0,0,0], 
			 [0,0,0,0,0,0,0],
			 [0,0,0,1,0,0,0],
			 [0,0,0,1,0,0,0],
			 [0,0,0,1,0,0,0], 
			 [0,0,0,0,0,0,0],
			 [0,0,0,0,0,0,0],
			 [0,0,0,0,0,0,0]],

			[[0,0,0,0,0,0,0], 
			 [0,0,1,0,0,0,0],
			 [0,1,0,0,0,0,0],
			 [0,1,0,0,0,0,0],
			 [0,0,1,0,1,0,0], 
			 [0,0,0,1,0,0,0],
			 [0,0,0,0,0,0,0],
			 [0,0,0,0,0,0,0]],

			[[0,0,0,0,0,0,0], 
			 [0,0,0,0,0,0,0],
			 [0,0,0,0,0,0,0],
			 [0,0,0,0,0,0,0],
			 [0,0,0,0,0,0,0], 
			 [0,0,0,0,0,1,0],
			 [0,0,0,0,1,0,0],
			 [0,0,0,0,0,0,0]]], dtype=np.bool)

		for order in ['C', 'F']: 
			self.runtest_thinning(order, expect)

if __name__ == '__main__':
	unittest.main()
