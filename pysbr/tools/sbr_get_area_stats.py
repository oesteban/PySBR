#!/usr/bin/env python
# coding: utf-8
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

"""

sbr_get_area_stats.py - Evaluate area stats

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

import os
import sys
import nibabel as nib
from math import sqrt
import numpy as np
import json 
from scipy.ndimage.morphology import binary_dilation
from scipy.ndimage.morphology import binary_closing 
import scipy.ndimage as ndimage

# add the setup build path before loading the module 
from testtools import append_buildpath
append_buildpath()
import csbr
from findsecspots import * 

def evaluate_stats(area_mask, hpoint):
        struct1 = ndimage.generate_binary_structure(3, 1)
        dmask = binary_dilation(area_mask, struct1)
        border = np.bitwise_xor(dmask, area_mask)
        return csbr.shapestats(border, hpoint)
        
def evaluate_labels_and_state_for_all(label_image, nspots, hotpoints):
        result = {}
        for n in range(nspots):
                hpoint = hotpoints[n]
                l = n+1
                area = (label_image == l)
                stats = evaluate_stats(area, hpoint)
                result[n] = stats
        return result

if __name__== '__main__':
    from argparse import ArgumentParser
    from argparse import RawTextHelpFormatter

    global_affine = None
    
    parser = ArgumentParser(description='Get the activation masks and their skeleton', 
                                             formatter_class=RawTextHelpFormatter)
    
    group = parser.add_argument_group('File I/O')
    
    group.add_argument("-i", "--input",  action="store", required=True,  
                       help="input image to be analyzed")
    
    group.add_argument("-I", "--inmask",  action="store", 
                       help="input image mask", dest='inmask')

    group.add_argument("-K", "--maskhelp",  action="store", 
                       help="mask helper image", dest='maskhelp')

    group = parser.add_argument_group('Pre-Processing')
    group.add_argument("-g", "--gauss-prefilter-sigma",  action="store", 
                       type=float, dest='gauss_sigma', 
                       help="Prefilter the image with this gaussian.")
    

    group = parser.add_argument_group('Processing')
    group.add_argument("-E", "--seed-intensity-ratio",  action="store", 
                       type=float, default=0.6, dest='rseedratio', 
                       help="relative intensity ratio with respect to the "
                       "maximum image intensity used to limit region grow seeds.")
    
    group.add_argument("-P", "--spot-mean-intensity-ratio",  action="store", 
                       type=float, default=1.4, dest='rspotmeanlimit', 
                       help="threshold to stop region growing with respect "
                       "to the ROI mean intensity ")
    
    group.add_argument("-B", "--spot-bridge-intensity-ratio",  action="store", 
                       type=float, default=0.7, dest='rspotbridgelimit', 
                       help="threshold to stop region growing with respect to"
                       " the minimum intensity found on the line connecting the two hotspots")

    group.add_argument("-U", "--uphill",  action="store", 
                       type=float, default=1.01, dest='uphill', 
                       help="factor to limit the region growing, the value describes by which "
                       "factor the intensity may increase to be still added in the region growing")
    
    volume_group = group.add_mutually_exclusive_group()
    
    volume_group.add_argument("-V", "--max-volume",  action="store", type=int, default=320000000, dest='maxvolume', 
                       help="Maximum volume (in pixel) to which a hot spot region can grow during during segmentation")
    
    volume_group.add_argument("-a", "--max-acceptable-volume",  action="store", type=int, default=6000, 
                              dest='max_acceptable_volume', help="Maximum activation volume (in pixel) "
                              "to accept for a segmentation, if larger, segmentation should be considered to have failed")

    landmark_group = parser.add_argument_group('Landmarks')
    landmark_group.add_argument("-l", "--landmark-add-volume",  action="store", type=int, default=30, 
                              dest='landmark_add_volume', help="volume threashold for adding additional landmarks")

    landmark_group.add_argument("-e", "--eigenvalue-ratio-add-volume",  action="store", 
                              type=float, default=3.0, dest='eigenvalue_ratio_add_volume', 
                              help="Eigenvalue ratio threshold for adding landmarks")

    debug_group = parser.add_argument_group('Debug')
    
    debug_group.add_argument("-d", "--debug",  action="store", 
                             help="store the output of the hotspot location and boundary combined with the image here")
                                
    options = parser.parse_args()
    
    
    try:
        study = nib.load(options.input)
        if options.inmask is not None: 
            inmask = nib.load(options.inmask)
            mask = inmask.get_data()
            pixcoord = np.where(mask > 0)
            bbox = (tuple([np.min(a) for a in pixcoord]), 
                    tuple([np.max(a) for a in pixcoord]))
        else:
            inmask = None 
            if options.maskhelp is not None:
                mask_creator = nib.load(options.maskhelp)
                mask = mask_creator.get_data()        
                pixcoord = np.where(mask > 0)
                bbox = (tuple([np.min(a) - 5 for a in pixcoord]), 
                        tuple([np.max(a) + 5 for a in pixcoord]))
                print(bbox) 
                
                
    
    except nib.spatialimages.ImageFileError as niberror:
        print ("Problem loading image: {}".format(niberror))
        
    else:
        orgimage = study.get_data()
        if study.get_header()['qform_code']>0:
            global_affine = study.get_qform()
        else:
            global_affine = study.get_affine()

        if inmask is not None: 
            m = inmask.get_data()
            image = np.zeros_like( orgimage )
            image[ m>0 ] = orgimage[ m>0]
        else:
            image = orgimage

        if options.gauss_sigma is not None:
            image = ndimage.gaussian_filter(image, sigma=options.gauss_sigma) 
        
        (labels, nspots, hotpoints, volumes) = csbr.findspots3d(image, rseedlimit=options.rseedratio, 
                                                      rspotmeanlimit=options.rspotmeanlimit, 
                                                      rspotbridgelimit=options.rspotbridgelimit, 
                                                      maxvolume=options.maxvolume, boundingbox=bbox)

        
        
        print("hotpoints=", hotpoints )
        (root, input_shortname) = os.path.split(options.input)
        for v in volumes:
            if v > options.max_acceptable_volume:
                print ("image({}): Segmentation of activation areas failed".format(input_shortname))
                sys.exit(1); 

        # improve hotspot position 
        hotpoints = adjust_spot_positions(image, labels, hotpoints, options.debug)
    
        print ("new hotpoints:",hotpoints)
    
        test_landmarks = find_secondary_spots3d_intensity_weighted( image, labels, hotpoints, options.landmark_add_volume, 
                                                                    options.eigenvalue_ratio_add_volume, True )

        print ("test_landmarks:",test_landmarks)
