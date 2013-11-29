#!/usr/bin/env python
# coding: utf-8
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

"""

test_linalg.py -

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


from  sbr_linalg import get_minimum_dist_point_from_plane

class TestSimple(unittest.TestCase):
    def test_closest_point_simple(self):
        hp = np.array([1,2,4])
        lm = np.array([1,2,-4])

        idx_list = []
        idx_list.append(np.array([1,3,1]))
        idx_list.append(np.array([2,3,2]))
        idx_list.append(np.array([2,4,3]))
        idx_list.append(np.array([1,4,0]))
        idx_list.append(np.array([2,3,-1]))

        md_point = get_minimum_dist_point_from_plane(hp, lm, idx_list)
        self.assertEqual(md_point[0], 1)
        self.assertEqual(md_point[1], 4)
        self.assertEqual(md_point[2], 0)

    def test_closest_point_two_in_plane(self):
        hp = np.array([1,2,4])
        lm = np.array([1,2,-4])

        idx_list = []
        idx_list.append(np.array([1,3,1]))
        idx_list.append(np.array([2,3,2]))
        idx_list.append(np.array([2,4,3]))
        idx_list.append(np.array([1,4,0]))
        idx_list.append(np.array([1,3,0]))
        idx_list.append(np.array([2,3,-1]))

        md_point = get_minimum_dist_point_from_plane(hp, lm, idx_list)
        self.assertEqual(md_point[0], 1)
        self.assertEqual(md_point[1], 3)
        self.assertEqual(md_point[2], 0)

    def test_closest_point_none_in_plane(self):
        hp = np.array([1,2,4])
        lm = np.array([1,2,-4])

        idx_list = []
        idx_list.append(np.array([1,3,1]))
        idx_list.append(np.array([2,3,2]))
        idx_list.append(np.array([2,4,3]))
        idx_list.append(np.array([2,3,-1]))

        md_point = get_minimum_dist_point_from_plane(hp, lm, idx_list)
        self.assertEqual(md_point[0], 1)
        self.assertEqual(md_point[1], 3)
        self.assertEqual(md_point[2], 1)

        
if __name__ == '__main__':
	unittest.main()


    
