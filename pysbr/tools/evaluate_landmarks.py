#!/usr/bin/env python
# coding: utf-8
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

"""

evaluate_landmarks.py - Shape-based regional statistics and metrics

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

def enum(**enums):
    return type('Enum', (), enums)


SyndromState = enum(ZERO=0, ONE=1, TWO=2, THREE=3)

import sys
import numpy as np
import ndimage
import scipy.morphology as ndmorph


def evaluate_stats(area_mask, hpoint):
        struct1 = ndimage.generate_binary_structure(2, 1)
        dmask = ndmorph.binare_dilation(area_mask, struct1)
        border = np.bitwise_xor(dmask, area_mask)
        return get_distance_stats(border, hpoint)
        
def evaluate_labels_and_state_for_all(label_image, nspots, hotpoints):
        result = []
        for n in range(nspots):
                hpoint = hotpoints[n]
                l = n+1
                area = (label_image == l)
                stats = evaluate_stats(area, hpoint))
                print ("gravdistance={}, dist-mean={}, dist-var={}".format(stats[0],stats[1],stats[2])


