#!/usr/bin/env python
# coding: utf-8
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

"""

findsecspots.py - Evaluate secondary landmarks from segmentation and intensity weights 

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
import nibabel as nib
from scipy.ndimage import binary_dilation,\
                          binary_closing,\
                          generate_binary_structure

# add the setup build path before loading the module 
from testtools import append_buildpath
append_buildpath()
import csbr

import os.path

global_affine = None 

def evaluate_main_direction_and_rate_it_weighted(area, spot):
    org_coords = np.argwhere( area != 0 )
    
    if np.shape( coords )[0] == 0:
        raise RuntimeError( 'Area is missing, volume=0' )

    # add a squared intensity weighting, the max intensity 
         
    intensities = area[np.where(area > 0)]
    max_val = np.max(intensities)
    delta_val = max_val - np.min(intensities)
    
    icurve = lambda x : 1.0 - ((max_val - x)/delta_val)**2

    coords = [ icurve(area[tuple(oc)]) * (oc-spot)  for oc in org_coords]
        
    maindir = np.sum(coords, 0)
    #print("main direction:", maindir)
    
    cov = np.cov(coords.T)
    w,u = np.linalg.eig(cov)
        # find largest and second largest value 
        
        
    max_evpos = 0
    sec_evpos = 1 
        
    if w[0] < w[1]:
        if w[0] < w[2]:
            if w[1] < w[2]:
                max_evpos = 2
                sec_evpos = 1
            else:
                max_evpos = 1
                sec_evpos = 2
        else:
            max_evpos = 1
            sec_evpos = 0
            
    else: # w[0] > w[1]

        if w[0] > w[2]:
            max_evpos = 0
            if w[1] < w[2]:
                sec_evpos = 2
            else:
                sec_evpos = 1
        else: # w[0] < w[2]
            max_evpos = 2
            sec_evpos = 0

    # correct the rirection of the eigenvector if it points 
    # the other way  
    if np.dot(u[max_evpos], maindir) < 0: 
        u = -u; 
    return (u[max_evpos], w[max_evpos], w[max_evpos] / w[sec_evpos])


def evaluate_main_directions_weighted_rated(image, label_image, hotpoints):
    result = []
    for n,spot in enumerate(hotpoints):
        area = np.zeros_like( image )
        label = n+1
        area[ label_image == label ] = image[ label_image == label]
        if np.sum(area.reshape(-1) ) > 0:
            result.append(evaluate_main_direction_and_rate_it(area, np.array(spot) ))
        else:
            print ( 'Warning: area associated to hotspot %d has zero volume' % n )

    return result

def adjust_spot_position(image, border, iscale, candidates, other_thin):
    """Evaluate the new spot as the point on the thinning that has the 
    maximum minimal distance to the border. 
    """
    if iscale == 0:
        raise RuntimeError( 'Division by 0.0 found in adjust_spot_position' )
 
    cand_coords   = np.argwhere( candidates == True )
    border_coords = np.argwhere( border == True )
    other_coords  = np.argwhere( other_thin == True )
    
    thin_location = {}

    for cc in cand_coords:
        sum_dist = 0
        for ci in np.array( cand_coords ):
            sum_dist = sum_dist + np.linalg.norm( ci-cc )
        
        omin = 10000
        for ci in other_coords:
            dist = np.linalg.norm( ci - cc ) 
            if dist < omin:
                omin = dist
        
        thin_location[tuple(cc)] = sum_dist / (omin + 1) / cand_coords.size
                    
    max_dist = 0
    max_dist_point = None

    dists = []

    for cc in cand_coords:
        min_dist = 100000
        
        for bc in border_coords:
            dist = np.linalg.norm( bc - cc )**2
            if dist < min_dist:
                min_dist = dist
        
        nmin_dist = min_dist * thin_location[tuple(cc)] * image[tuple(cc)] / iscale

        if nmin_dist > max_dist:
            dists+=  [nmin_dist ]
            max_dist = nmin_dist
            max_dist_point = tuple(cc)

    return max_dist_point



def adjust_spot_positions(image, label_image, hp, debug=None):
    """Re-evaluate the spot positions based on the segmentation. 
        Parameters: 
        image: The original image (can be masked) that was sent to findspot3d
        label_image: the label image containing two labels 
        hp: the original hotpoints
        debug: set to true to write out an image debugimg.nii.gz with the stuff
        """
        
    struct2 = generate_binary_structure(3, 2)
    struct1 = generate_binary_structure(3, 1)
    peak_points =[] 

    if debug is None:
        temp_path = os.getenv("PYSBR_TEMP")
        if temp_path is not None:
            debug = os.path.join(temp_path, "debug-labels.nii.gz")
    
    if debug is not None:
            debimg = image.copy()

    nlabels = label_image.max()

    if nlabels!=len(hp):
        raise RuntimeError( 'number of labels and hotspots should be the same' )

    tins = []
    for n in range(nlabels):
        label = n+1
        area = binary_closing(label_image == label, struct2)
        thiniter = np.sum(area.reshape(-1)) / 1500 + 1
        csbr.thinning3d(area, thiniter)
        tins.append(area)
   
    for n in range(nlabels):
        label = n+1
        
        #avoid that a single pixel breaks the evaluation by running a closing 
        area = label_image == label
        
        #evaluate the boundary 
        dmask = binary_dilation(area, struct1)
        border = np.bitwise_xor(dmask, area)
        
        p = adjust_spot_position(image, border, image[tuple(hp[n])], tins[n], tins[(n + 1) % 2])
        peak_points.append(p)

        if debug is not None:
            debimg[border>0] = 196
            debimg[p] = 0
            nib.save(nib.Nifti1Image(debimg, global_affine), debug)

    peak_points = np.array( peak_points )
    return peak_points


def find_secondary_spots3d_intensity_weighted( image, labels, hotspots, minvolume, evratio, lm_fixed_steps, adjust=False ):
    """ This function finds additional landmarks. If two hostspots are given then it adds one 
    in the middle of these two, then it addes landmarks to the hotspots along the direction
    of the main PCA eigenvector of the segmented area weighted by the pixel intensity according either 
    at fixed steps or according  to the size of the first eigenvalue. 
    Input values are
    image: original image 
    labels: the image containing the labeled areas
    hotspots: the coordinates corresponding to the hotspots in the cativation areas 
    minvolume: the minimum volume at which we start adding secondary landmarks 
    evratio: the minimum ration of the 1st ans 2nd eigenvalue at which we start adding landmarks 
    lm_fixed_steps: if True use a fixed step-width of 3.0 to add the landmarks, if False, use the 
    eigenvalue as scaling factor
    adjust: if True re-locates the 1st hotspot to a better position

    Returns: A map that lists all hotspot based landmarks like follows: 
       index 0: a list including the point in the middle between hotspots[0] and hotspots[1]
       index 1: a list starting with the 1st hotspot, then its secondary landmarks  
       index 2: a list starting with the 2st hotspot, then its secondary landmarks  
    """


    dirs_rated = evaluate_main_directions_weighted_rated( image, labels, hotspots )
    histo = np.bincount( labels.reshape(-1) )

    newspots = {}

    hotspots = np.array( np.atleast_2d( hotspots ), dtype=np.float32 )

    if adjust:
        # Improve hotvox locations
        hotspots = adjust_spot_positions( image, labels, hotspots )

    if len(hotspots) == 2:
        newspots['middle'] = np.average( hotspots, axis=0 )



    # check volume and spot ratio, add at most two additional landmarks  
    for i,d in enumerate( dirs_rated ):
        idx = i+1
        key = 'spots%d'%idx
        newspots[key] = [hotspots[i]]
        volume = histo[idx]
        n_additional_landmarks_volume_based = volume / minvolume
        if n_additional_landmarks_volume_based > 2:
            n_additional_landmarks_volume_based = 2

        n_additional_landmarks_ev_based = d[2] / evratio
        if n_additional_landmarks_ev_based > n_additional_landmarks_volume_based:
            n_additional_landmarks_ev_based = n_additional_landmarks_volume_based

        if n_additional_landmarks_ev_based >= 1:
            if (lm_fixed_steps):
                addon_spot = hotspots[i] + 3.0 * d[0]
                newspots[key].append( addon_spot )
            else:
                addon_spot = hotspots[i] + 0.5 * i * d[1] * d[0]
                newspots[key].append( addon_spot )


    return newspots


def evaluate_main_direction_and_rate_it(area, spot):
        coords = np.argwhere( area>0 )
        vol = np.shape( coords )[0]
        if vol == 0:
            raise RuntimeError( 'Area is missing, volume=0' )

        # single pixel area, bail out  
        if vol == 1:
            return (np.array([1.0,0.0,0.0]), 0, 0)

        coords -= spot
        
        maindir = np.sum(coords, 0)
        #print("main direction:", maindir)
        
        cov = np.cov(coords.T)
        w,u = np.linalg.eig(cov)
        # find largest and second largest value 
        
        
        max_evpos = 0
        sec_evpos = 1 

        if w[0] < w[1]:
                if w[0] < w[2]:
                        if w[1] < w[2]:
                                max_evpos = 2
                                sec_evpos = 1
                        else:
                                max_evpos = 1
                                sec_evpos = 2
                else:
                        max_evpos = 1
                        sec_evpos = 0
                        
        else: # w[0] > w[1]
                
                if w[0] > w[2]:
                        max_evpos = 0
                        if w[1] < w[2]:
                                sec_evpos = 2
                        else:
                                sec_evpos = 1
                else: # w[0] < w[2]
                        max_evpos = 2
                        sec_evpos = 0

        if np.dot(u[max_evpos], maindir) < 0: 
            u = -u; 
        return (u[max_evpos], w[max_evpos], w[max_evpos] / w[sec_evpos])


def evaluate_main_directions_rated(label_image, hotpoints):
    result = []
    for n,spot in enumerate(hotpoints):
        area = (label_image == n+1).astype( np.uint8 )
        if np.sum(area.reshape(-1) ) > 0:
            result.append(evaluate_main_direction_and_rate_it(area, np.array(spot) ))
        else:
            print ( 'Warning: area associated to hotspot %d has zero volume' % n )

    return result


def findsecondaryspots3d( labels, hotspots, minvolume, evratio, lm_fixed_steps ):
    """ This function finds additional landmarks. If two hostspots are given then it adds one 
    in the middle of these two, then it addes landmarks to the hotspots along the direction
    of the main PCA eigenvector of the segmented area according either at fixed steps or according 
    to the size of the first eigenvalue. 
    Input values are
    labels: the image containing the labeled areas
    hotspots: the coordinates corresponding to the hotspots in the cativation areas 
    minvolume: the minimum volume at which we start adding secondary landmarks 
    evratio: the minimum ration of the 1st ans 2nd eigenvalue at which we start adding landmarks 
    lm_fixed_steps: if True use a fixed step-width of 3.0 to add the landmarks, if False, use the 
    eigenvalue as scaling factor       
    Returns: A map that lists all hotspot based landmarks like follows: 
       index 0: a list including the point in the middle between hotspots[0] and hotspots[1]
       index 1: a list starting with the 1st hotspot, then its secondary landmarks  
       index 2: a list starting with the 2st hotspot, then its secondary landmarks  
    """


    dirs_rated = evaluate_main_directions_rated( labels, hotspots )
    histo = np.bincount( labels.reshape(-1) )

    newspots = {}

    hotspots = np.array( np.atleast_2d( hotspots ), dtype=np.float32 )

    if len(hotspots) == 2:
        newspots['middle'] = np.average( hotspots, axis=0 )


    # check volume and spot ratio, add at most two additional landmarks  
    for i,d in enumerate( dirs_rated ):
        idx = i+1
        key = 'spots%d'%idx
        newspots[key] = [hotspots[i]]
        volume = histo[idx]
        n_additional_landmarks_volume_based = volume / minvolume
        if n_additional_landmarks_volume_based > 2:
            n_additional_landmarks_volume_based = 2

        n_additional_landmarks_ev_based = d[2] / evratio
        if n_additional_landmarks_ev_based > n_additional_landmarks_volume_based:
            n_additional_landmarks_ev_based = n_additional_landmarks_volume_based

        if n_additional_landmarks_ev_based >= 1:
            if (lm_fixed_steps):
                addon_spot = hotspots[i] + 3.0 * d[0]
                newspots[key].append( addon_spot )
            else:
                addon_spot = hotspots[i] + 0.5 * i * d[1] * d[0]
                newspots[key].append( addon_spot )

    return newspots
