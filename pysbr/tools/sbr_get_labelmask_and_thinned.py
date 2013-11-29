#!/usr/bin/env python
# coding: utf-8
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

"""

sbr_get_labelmask_and_thinned.py - 

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
import numpy as np

# add the setup build path before loading the module 
from testtools import append_buildpath
append_buildpath()
import csbr


from argparse import ArgumentParser
from argparse import RawTextHelpFormatter

parser = ArgumentParser(description='Get the activation masks and their skeleton', 
                                         formatter_class=RawTextHelpFormatter)

group = parser.add_argument_group('File I/O')

group.add_argument("-i", "--input",  action="store", required=True,  help="input image to be analyzed")
group.add_argument("-I", "--inmask",  action="store", required=True,  help="input image mask", dest='inmask')

group.add_argument("-m", "--maskimage",   action="store", required=True,  dest='maskfile', 
                   help="output image containing the activation masks")
group.add_argument("-t", "--thin-maskimage",   action="store", required=False,  dest='thinfile', 
                   help="output image containing the thinned activation masks")

group.add_argument("-s", "--spots",   action="store", required=True,  dest='spotfile', 
                   help="output text file containing the hot spots used for mask initialization")


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
    
group.add_argument("-V", "--max-volume",  action="store", type=int, default=320000000, dest='maxvolume', 
                   help="Maximum volume (in pixel) to which a hot spot region can grow during during segmentation")

options = parser.parse_args()


try:
    study = nib.load(options.input)
    if options.inmask is not None: 
        inmask = nib.load(options.inmask)
    else:
        inmask = None

except nib.spatialimages.ImageFileError as niberror:
    print ("Problem loading image: {}".format(niberror))
    
else:
    orgimage = study.get_data(); 
    if inmask is not None: 
        m = inmask.get_data()
        image = np.zeros_like( orgimage )
        image[ m>0 ] = orgimage[ m>0]
    else:
        image = orgimage

    try:
        (labels, nspots, hotpoints, volumes) = csbr.findspots3d(image, rseedlimit=options.rseedratio, 
                                                     rspotmeanlimit=options.rspotmeanlimit, 
                                                     rspotbridgelimit=options.rspotbridgelimit, 
                                                     maxvolume=options.maxvolume)
    except csbr.error as e:
        print("Wrror running csbr.findspots3d: {} {}".format(e.errno, e.strerror))
        sys.exit(1)

    if options.maskfile is not None:
        try: 
                nib.save(nib.Nifti1Image(labels, study.get_affine()), options.maskfile)
        except IOError as e:
                print ("I/O error saving to '{0}'({1}): {2}".format(options.maskfile, e.errno, e.strerror))

    if options.spotfile is not None:
        f = open(options.spotfile, 'w')
        for h in hotpoints:
            f.write("{} {} {}\n".format(h[0], h[1], h[2]))
        f.write             
                
    if options.thinfile is not None:
            mask1 = np.equal(labels, 1)
            result1 = csbr.thinning3d(mask1)
            mask2 = np.equal(labels, 2)
            result2 = csbr.thinning3d(mask2)
            
            thin = np.logical_or(mask1, mask2)
            if result1 == 0 and result2 == 0:
                try:    
                    nib.save(nib.Nifti1Image(np.asarray(thin, dtype=np.byte), study.get_affine()), options.thinfile)
                except IOError as e:
                    print ("I/O error saving to '{0}'({1}): {2}".format(options.maskfile, e.errno, e.strerror))
            else:
                    print ("Thinning failed")
        

        

    
    






