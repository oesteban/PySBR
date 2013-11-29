#!/usr/bin/env python
# coding: utf-8
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

"""

sbr_get_landmarks_without_tests.py
 
     This program takes an activation image and evaluates the landmarks 
     without considering whether they actually make sense. 
     Output is a JSON file with the landmark locations for the two actuvation 
     areas.



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



import sys
import nibabel as nib
from math import sqrt
import numpy as np
import json 

# add the setup build path before loading the module 
from testtools import append_buildpath
append_buildpath()
import csbr
import sbr_linalg 


# this should go under test 
def get_landmarks_from_mask(hotspot, mask):
    thinned = csbr.thinning3d(mask)
    hp = np.array([hotpoints[l][0], hotpoints[l][1], hotpoints[l][2]], dtype=np.float)
    landmarks = [hp]

    idx = np.where(mask == True)

    # there must be a better way
    max_dist = 0
    max_dist_idx = None
    idx_list = []
    for i in range(len(idx[0])):
        ii = np.array([idx[0][i], idx[1][i], idx[2][i]])
        idx_list.append(ii)
        delta = ii - hp 
        dist = np.sum(delta * delta)
        if dist > max_dist:
            max_dist = dist
            max_dist_idx = ii
    
    lm_large = max_dist_idx
        
    lm_mid = sbr_linalg.get_minimum_dist_point_from_plane(hp, lm_large, idx_list)

    landmarks.append(lm_mid)
    landmarks.append(lm_large)
    return landmarks


# since JSON doesn't support numpy arrays the result must be translated to 
# a list of tuples 
def save_landmarks(filename, lmlist):

    outlist = []
    for lms in lmlist:
        sublist = []
        for lm in lms: 
            sublist.append((lm[0], lm[1], lm[2]))
        outlist.append (sublist)
    
    
    f = open(filename, 'w')
    f.write(json.dumps(outlist, indent=2))
    f.close()
    
    
from argparse import ArgumentParser
from argparse import RawTextHelpFormatter

parser = ArgumentParser(description='Get the activation masks and their skeleton', 
                                         formatter_class=RawTextHelpFormatter)

group = parser.add_argument_group('File I/O')

group.add_argument("-i", "--input",  action="store", required=True,  help="input image to be analyzed")
group.add_argument("-o", "--landmarks",   action="store", required=True,  dest='landmarkfile', 
                   help="output JSON file containing three landmarks for each spot")


group = parser.add_argument_group('Processing')
group.add_argument("-E", "--seed-intensity-ratio",  action="store", type=float, default=0.6, dest='rseedratio', 
                   help="relative intensity ratio with respect to the maximum image intensity used to limit region grow seeds ")

group.add_argument("-P", "--spot-mean-intensity-ratio",  action="store", type=float, default=1.5, dest='rmeanspotratio', 
                   help="with respect to the spot seed intensity to stop region growing")

group.add_argument("-B", "--spot-bridge-intensity-ratio",  action="store", type=float, default=1.1, dest='rbridgespotratio', 
                   help="with respect to the spot seed intensity to stop region growing")


group.add_argument("-U", "--uphill",  action="store", type=float, default=1.1, dest='uphill', 
                   help="factor to limit the region growing, the value describes by which factor the intensity may increase to be still added in the region growing")

volume_group = group.add_mutually_exclusive_group()

volume_group.add_argument("-V", "--max-volume",  action="store", type=int, default=320000000, dest='maxvolume', 
                   help="Maximum volume a hot spot region can assume during segmentation")

volume_group.add_argument("-m", "--max-acceptable-volume",  action="store", type=int, default=300, dest='max_sensible_volume', 
                   help="Maximum activation volume to accept for a segmentation")

options = parser.parse_args()


try:
    study = nib.load(options.input)

except nib.spatialimages.ImageFileError as niberror:
    print ("Problem loading image: {}".format(niberror))
    
else:
    image = study.get_data(); 

    (labels, nspots, hotpoints, volumes) = csbr.findspots3d(image, rseedlimit=options.rseedratio, 
                                                  rspotmeanlimit=options.rmeanspotratio, 
                                                  rspotbridgelimit=options.rbridgespotratio, 
                                                  maxvolume=options.maxvolume)

    for v in volumes: 
        if v > options.max_sensible_volume:
            print "Hotspot segmentation failed; decision based on hotspot volume"
            sys.exit(1)

    lm_sets = []
    for l in range(nspots):
        mask = np.equal(labels, l+1)
        landmarks  = get_landmarks_from_mask(hotpoints[l], mask)

        lm_sets.append(landmarks)

save_landmarks(options.landmarkfile, lm_sets)

