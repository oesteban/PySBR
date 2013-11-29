#!/usr/bin/env python
# coding: utf-8
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

"""
Process.py: This file gathers the function calls and auxiliary methods
            for PySBR


Copyright (c) 2013, code@oscaresteban.es (Oscar Esteban), and 
                    gw.fossdev@gmail.com (Gert Wollny)
                    with Biomedical Image Technology, UPM (BIT-UPM)
All rights reserved.

This file is part of PySBR.

"""

__author__ = "Oscar Esteban"
__copyright__ = "Copyright 2013, Biomedical Image Technologies (BIT), \
                 Universidad Politécnica de Madrid"
__credits__ = ["Oscar Esteban", "Gert Wollny"]
__license__ = "FreeBSD"
__version__ = "0.1"
__maintainer__ = "Oscar Esteban"
__email__ = "code@oscaresteban.es"
__status__ = "Prototype"


import nibabel as nib
import numpy as np
import scipy.ndimage as im
import os.path as op
import os
import math

import sys
sys.path.insert(0, op.abspath( './thinning3d/build/lib.linux-x86_64-2.7/') )
import csbr

def find_secondary_spots( nuclei, hotspots ):
    coms = []
    furthest = [ (0,0,0), (0,0,0) ]
    maxdists_1 = [ 0, 0 ]
    second_furthest = [ (0,0,0), (0,0,0) ]
    maxdists_2 = [ 0, 0 ]

    for hotspotid in [ 1, 2 ]:
        nucl_msk = np.zeros( shape=nuclei.shape )
        nucl_msk[ nuclei==hotspotid ] = 1
        coms.append( tuple( [ int(val) for val in im.measurements.center_of_mass( nucl_msk )]) )
        
        msk_eroded = im.binary_erosion( nucl_msk ).astype( np.uint8 )
        boundary = nucl_msk - msk_eroded
        bids = np.swapaxes( np.array( np.where( boundary==1 ) ).astype( np.uint8 ), 0, 1)
           
        for bid in bids:
            dist = np.linalg.norm( bid - hotspots[hotspotid-1] )
            if dist > maxdists_1[hotspotid-1]:
                maxdists_1[hotspotid-1]=dist
                furthest[hotspotid-1] = tuple( bid )
                
        normal = furthest[hotspotid-1] - hotspots[hotspotid-1]
    
        for bid in bids:
            v = bid - hotspots[hotspotid-1]
            if math.fabs( np.dot( v, normal ) ) < 1e-2:
                dist = np.linalg.norm( bid - hotspots[hotspotid-1] )
                if dist > maxdists_2[hotspotid-1]:
                    maxdists_2[hotspotid-1]=dist
                    second_furthest[hotspotid-1] = tuple( bid )

    return (furthest, second_furthest)
