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
import numpy as np


if __name__== '__main__':
    from argparse import ArgumentParser
    from argparse import RawTextHelpFormatter
    
    
    parser = ArgumentParser(description='Evaluate dice index for two given masks', 
                                             formatter_class=RawTextHelpFormatter)


    group = parser.add_argument_group('File I/O')
    
    group.add_argument("-s", "--source",  action="store", required=True,  help="input mask to be analyzed")
    group.add_argument("-r", "--reference",  action="store", required=True,  help="reference mask ")

    options = parser.parse_args()

    try:
        src_data = nib.load(options.source)
        ref_data = nib.load(options.reference)

        src = src_data.get_data().astype(np.bool); 
        ref = ref_data.get_data().astype(np.bool);

        union = np.logical_or(src, ref)
        schnitt = np.logical_and(src, ref)

        print ((1.0 * np.count_nonzero(schnitt))/ np.count_nonzero(union))
        
    except nib.spatialimages.ImageFileError as niberror:
        print ("Problem loading image: {}".format(niberror))


