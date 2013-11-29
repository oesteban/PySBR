#!/usr/bin/env python
# coding: utf-8
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

"""

sbr_linalg.py -

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
from math import sqrt
from math import fabs


def get_distance(p1, p2):
    delta = p1 - p2
    return sqrt(np.sum(delta * delta))
    

def get_minimum_dist_point_from_plane(hotpoint, landmark, idx_list):
    """Evaluate the plane that is orthogonal to the line through hotpoint 
    and landmark and select the point that is closest to this plane out of the idx_list.
    It two ploits have the same distance from the plane, select the one closer to the plane 
    foot point, which is the midlle point of the line from hotpoint to landmark. 
    hotpoint and landmark are given as numpy arrays and idx_list is a list of numpy arrays
    """
    direction = 0.5 * (hotpoint - landmark)
    foot = direction + landmark
    n = direction / sqrt(np.sum(direction * direction))
    d = -np.inner(n, foot)
    
    min_dist_point = None
    min_dist = 1e+32
    min_dist_to_foot = None
    
    for idx in idx_list:
        dist = fabs(np.inner(idx, n) + d); 
        if dist < min_dist:
            min_dist = dist
            min_dist_point = idx
            min_dist_to_foot = get_distance(idx, foot)
            
        elif dist == min_dist:
            dist_to_foot = get_distance(idx, foot)
            if dist_to_foot < min_dist_to_foot:
                min_dist_point = idx
                min_dist_to_foot = get_distance(idx, foot)
    
    return  min_dist_point


