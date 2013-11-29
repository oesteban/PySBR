#!/usr/bin/env python
# coding: utf-8
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

"""
shape.py - Nipype compliant interface for landmark extraction

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


import os, os.path as op
import warnings
import numpy as np
import nibabel as nib
import scipy.ndimage as si
from nipype.interfaces.base import (TraitedSpec, File, InputMultiPath,
                                    OutputMultiPath, Undefined, traits,
                                    isdefined, OutputMultiPath, 
                                    CommandLineInputSpec, CommandLine,
                                    BaseInterface, BaseInterfaceInputSpec,
                                    traits, Directory )
from nipype.utils.filemanip import split_filename,fname_presuffix

import json

import csbr
import pysbr.tools.polar as polar
from pysbr.tools.sbr_get_area_stats import adjust_spot_positions
from pysbr.tools.findsecspots import find_secondary_spots3d_intensity_weighted

import itk
import matplotlib.pyplot as plt
import matplotlib.cm as cm



class FindSpotsInputSpec(BaseInterfaceInputSpec):
    in_file = File( exists=True, mandatory=True,
                    desc="Target study to be processed" )
    in_roi = File( exists=True, mandatory=True,
                   desc="ROI associated to the activity")
    rseedlimit = traits.Float( 0.6, usedefault=True,
                               desc='relative intensity (from mean image intensity) for possible hotspot seed pixel.' )
    rspotmeanlimit = traits.Float( 1.6, usedefault=True,
                                   desc='relative intensity for region growing (from mean masked area image intensity)' )
    rspotbridgelimit = traits.Float( 0.7, usedefault=True,
                                     desc='relative intensity for region growing (from minimum intensity of the pixels that form the bridge between the hotspots)' )
    maxvolume = traits.Int( 30000, usedefault=True, desc='Maxium volume of a hotspot area in pixels' )
    minvolume = traits.Int( 30, usedefault=True, desc='Minimum volume factor to consider secondary landmarks')
    evratio = traits.Float( 1.5, usedefault=True, desc='ratio multiplier to decide how many secondary landmarks are used based on the ratio of the first and second PCA eigenvalue') 
    
    uphill = traits.Float( 1.0, usedefault=True,
                           desc='factor that describes in region growing by which factor the intensity may increase to still add to the region, should be >= 1.0' )
    maxacceptablevol = traits.Int( 30000, usedefault=True, desc='Maximum acceptable volume' )
#    bbox = traits.Enum(
#                           desc='tuple ((x1,y1,z1),(x2,y2,z2)) that defines a bounding box to restrict the search for spots.' )
    out_labels = File( desc='output segmentation filename' )
    out_spotfile = File( desc='output hotspots filename' )

    gauss_sigma = traits.Float( 0.0, usedefault=True, desc="Smooth the input image before running the point segmentation")

class FindSpotsOutputSpec(TraitedSpec):
    out_labels   = File( desc='Output Segmentation of activity nuclei' )
    out_spotfile = File( desc='Spots coordinates' )

class FindSpots(BaseInterface):
    """
    find two hostspots in the  given image and return the label masks of the areas 

    Example
    -------

    >>> import pysbr_nipype as pysbr
    >>> fspots = pysbr.FindSpots()
    >>> fspots.inputs.in_file = tfilename
    >>> fspots.inputs.in_roi = tmplmsk
    >>> fspots.inputs.out_labels = maskfile
    >>> fspots.inputs.out_spotfile = spotfile
    >>> res1 = fspots.run()  # doctest: +SKIP

    Inputs
    ------

        image: the image to search the two spots in
        rseedlimit: relative intensity (from mean image intensity) for possible hotspot seed pixel.
        rspotmeanlimit: relative intensity for region growing (from mean masked area image intensity)
        rspotbridgelimit: relative intensity for region growing (from minimum intensity of the pixels that form the bridge between the hotspots)
        maxvolume: Maxium volume of a hotspot area in pixels
        maxacceptablevol: Maximum acceptable volume
        bbox: [NOT IMPLEMENTED] tuple ((x1,y1,z1),(x2,y2,z2)) that defines a bounding box to restrict the search for spots
        uphill: factor that describes in region growing by which factor the intensity may increase to sill add to the region, should be >= 1.0

    Outputs
    -------
        out_labels: the label mask image containing the hot areas. Lower label number corresponds to the hot spot with higher intensity.
        out_spotfile: text file with the pixel index of the location of the two spots

    """
    input_spec = FindSpotsInputSpec
    output_spec = FindSpotsOutputSpec

    def _run_interface(self, runtime):
        im = nib.load( self.inputs.in_file )
        roi = nib.load( self.inputs.in_roi )
        im_data = im.get_data()
        roi_data = roi.get_data()
        im_data[roi_data==0] = 0

        # Read images
        imread = itk.ImageFileReader[itk.Image[itk.F,3]].New()
        imread.SetFileName( self.inputs.in_file )
        imread.Update()
        itkim = imread.GetOutput()
        IndexType = itk.ContinuousIndex[ itk.D, 3 ]
        
        maxaccvolume = self.inputs.maxacceptablevol
        minvolume = self.inputs.minvolume
        evratio = self.inputs.evratio
        lm_fixed_steps = True # TODO expose this parameter

        if self.inputs.gauss_sigma > 0.0:
            im_data = si.gaussian_filter(im_data, sigma=self.inputs.gauss_sigma) 

        (labels, nspots, hotvox, volumes) = csbr.findspots3d(im_data, 
                                                             rseedlimit=self.inputs.rseedlimit,
                                                             rspotmeanlimit=self.inputs.rspotmeanlimit, 
                                                             rspotbridgelimit=self.inputs.rspotbridgelimit,
                                                             maxvolume=self.inputs.maxvolume)

        
        
        vol_check = [ v<maxaccvolume for v in volumes ]
        if not np.all( vol_check ):
            raise RuntimeError( 'Segmentation of activation areas failed' )

        hotvox = adjust_spot_positions(im_data, labels, hotvox)

        print ("Corrected hotspots are ", hotvox)
        #landmark_map = findsecondaryspots3d( labels, hotvox, minvolume, evratio, lm_fixed_steps ) 
        landmark_map = find_secondary_spots3d_intensity_weighted( im_data, labels, hotvox, minvolume, evratio, lm_fixed_steps, adjust=True ) 

        # translate the landmarks to physical coordinates 
        for k,v in landmark_map.iteritems():
            v = np.atleast_2d( v )
            landmark_map[k] = [ tuple( itkim.TransformContinuousIndexToPhysicalPoint( IndexType( [ float(val) for val in vox ] ) ) ) for vox in v ]

        temp_path = os.getenv("PYSBR_TEMP")
        # Find segmentations' pointcloud
        for l in range(1,nspots+1):
            seg = np.zeros_like( labels )
            seg[labels==l] = 1
            boundary = si.binary_dilation( seg ) - seg
            if temp_path is not None: 
                nib.save( nib.Nifti1Image( boundary, roi.get_affine(), roi.get_header() ), 
                          os.path.join(temp_path, ('label%d.nii.gz' % l)))
                          
            idxs = np.argwhere( boundary==1 )
            landmark_map['surf%d'%l] = [ tuple( itkim.TransformContinuousIndexToPhysicalPoint( IndexType( [ float(val) for val in vox ] ) ) ) for vox in idxs ]

        # Generate outputs
        if not isdefined( self.inputs.out_labels ):
            self.inputs.out_labels = self._gen_fname( self.inputs.in_file,
                                                      suffix='labels',
                                                      ext=None )

        if not isdefined( self.inputs.out_spotfile ):
            self.inputs.out_spotfile = self._gen_fname( self.inputs.in_file,
                                                        suffix='spots', 
                                                        ext='.txt' )
        with open( self.inputs.out_spotfile, 'w') as f:
            json.dump(landmark_map, f)

        nib.save(nib.Nifti1Image(labels, im.get_affine(), im.get_header() ), self.inputs.out_labels )
        return runtime
                        
    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs["out_labels"] = self.inputs.out_labels
        outputs["out_spotfile"] = self.inputs.out_spotfile
        return outputs

    def _gen_fname( self, refname, suffix=None, ext=None ):
        basename = op.basename( refname )
        
        refname, refext = op.splitext( basename )
        if refext =='.gz':
            refname, refext2 = op.splitext(refname)
            refext = refext2 + refext

        if ext is None:
            ext = refext

        if suffix is None:
            suffix = ''
        else:
            suffix = '_'+suffix

        fname = op.abspath( './' + refname + suffix + ext )
        return fname


class SBRInputSpec(BaseInterfaceInputSpec):
    src_spots = File( exists=True, mandatory=True,
                      desc="Source spots describing template shape" )
    trg_spots = File( exists=True, mandatory=True,
                      desc="Target spots describing subject's shape" )
    nlabels = traits.Int( 2, mandatory=True, usedefault=True,
                          desc='Number of labels to search for' )
    out_transform = File( dec='Output JSON file' )

class SBROutputSpec(TraitedSpec):
    out_transform = File( desc='Source spots location in target space' )

class SBR(BaseInterface):
    """

    Example
    -------

    >>> import pysbr_nipype as pysbr
    >>> sbreg = pysbr.SBR()
    >>> sbreg.inputs.src_spots = source.txt
    >>> sbreg.inputs.trg_spots = target.txt
    >>> sbreg.inputs.out_transform = tfmap.txt
    >>> res = sbreg.run()  # doctest: +SKIP

    Inputs
    ------


    Outputs
    -------

    """
    input_spec = SBRInputSpec
    output_spec = SBROutputSpec

    def _run_interface(self, runtime):
        from scipy.ndimage import affine_transform

        temp_path = os.getenv("PYSBR_TEMP")

        finaljson = { 'moving': [], 'static': [] }

        with open( self.inputs.trg_spots ) as f:
            json_t = json.load(f)

        with open( self.inputs.src_spots ) as f:
            json_s = json.load(f)

        for label in range(1,self.inputs.nlabels + 1):
            ht = np.array( json_t['spots%d'%label][0], dtype=np.float32 )
            hs = np.array( json_s['spots%d'%label][0], dtype=np.float32 )

            pointmapper = { 'source': [ json_s['spots%d'%label][0] ], 'target': [ json_t['spots%d'%label][0] ] }

            points_t = np.array( json_t['surf%d'%label] )
            points_s = np.array( json_s['surf%d'%label] )

            # Feature matching core
            X,Y,S = polar.polar_rep( points_s, hs, interp='nearest' ) # Polar map
            X,Y,T = polar.polar_rep( points_t, ht, interp='nearest' ) # Polar map

            for polarmap,name in zip( [ S, T ], ['source','target'] ):
                f = plt.figure()
                plt.imshow( polarmap, cmap=cm.gray )
                if temp_path is not None:
                    plt.savefig( os.path.join(temp_path, ('label%d_' % label + name + '.png')))

            inc = np.array( [ np.pi / np.shape(S)[0], 2.0*np.pi / np.shape(S)[1] ] )  # Polarmaps pixelsize
            reg = polar.Registration( S, T ) # SHAPE-BASED registration
            res = reg.optimize()             # Optimize :)
            angle = res.x * inc           # Convert polarmap pixels to radians


            im1_tf = affine_transform(T, np.identity(2) ,mode='wrap',offset=-res.x )
            f = plt.figure()
            plt.imshow( im1_tf, cmap=cm.gray )
            if temp_path is not None:
                plt.savefig( os.path.join(temp_path, ('label%d_' % label + 'target_tfmd.png')))

            if np.any( angle > 0.5 ) or (points_t.size < np.round(0.75 * points_s.size)):
                angle = [0.0,0.0]

            print 'Angular correction: ',angle,' translation=',(ht-hs)

            # newpoints = []
            # for p in points_s:
            #     pointmapper['source'].append( tuple(p) )
            #     p = p - hs
            #     r = np.linalg.norm( p )
            #     theta = np.arccos( p[2]/r ) + angle[0]
            #     phi = np.arctan2(p[1],p[0]) + angle[1]

            #     # Rotate and translate to NEW center (target hotspot)
            #     newp = (np.array( [ np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta) ] ) * r) + ht
            #     pointmapper['target'].append( tuple(newp) )

            finaljson['moving'].append( pointmapper )

        # Add fixedpoints
        finaljson['static'] = json_s['static']

        if not isdefined( self.inputs.out_transform ):
            self.inputs.out_transform = self._gen_fname( self.inputs.src_spots,
                                                       suffix='map',
                                                       ext='.txt' )        

        with open( self.inputs.out_transform, 'w') as f:
            json.dump(finaljson, f)

        return runtime

                        
    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs["out_transform"] = self.inputs.out_transform
        return outputs

    def _gen_fname( self, refname, suffix=None, ext=None ):
        basename = op.basename( refname )
        
        refname, refext = op.splitext( basename )
        if refext =='.gz':
            refname, refext2 = op.splitext(refname)
            refext = refext2 + refext

        if ext is None:
            ext = refext

        if suffix is None:
            suffix = ''
        else:
            suffix = '_'+suffix

        fname = op.abspath( './' + refname + suffix + ext )
        return fname



class ThinningInputSpec(BaseInterfaceInputSpec):
    in_file = File( exists=True, mandatory=True,
                    desc="segmentation file" )
    out_file = File( desc="output file with thinning for all labels" )

class ThinningOutputSpec(TraitedSpec):
    out_file = File( desc="output file with thinning for all lables" )

class Thinning(BaseInterface):
    input_spec = ThinningInputSpec
    output_spec = ThinningOutputSpec

    def _run_interface(self, runtime):
        im = nib.load( self.inputs.in_file )
        im_data = im.get_data()
        nlabels = np.max( im_data.reshape(-1) )
        result = np.zeros( shape=im_data.shape )

        for l in range(1,nlabels):
            mask = np.equal(im_data, l)
            exitcode = csbr.thinning3d(mask)
            if exitcode !=0:
                raise IOError( 'PySBR thinning failed for label %d' % l )  
            result = result + mask.astype( np.uint8 )

        if not isdefined( self.inputs.out_file ):
            self.inputs.out_file = self._gen_fname( self.inputs.in_file,
                                                    suffix='thin' )

        nib.save(nib.Nifti1Image( result, im.get_affine(), im.get_header() ), self.inputs.out_file )
        return runtime

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs["out_file"] = self.inputs.out_file
        return outputs

    def _gen_fname( self, refname, suffix=None, ext=None ):
        basename = op.basename( refname )
        
        refname, refext = op.splitext( basename )
        if refext =='.gz':
            refname, refext2 = op.splitext(refname)
            refext = refext2 + refext

        if ext is None:
            ext = refext

        if suffix is None:
            suffix = ''
        else:
            suffix = '_'+suffix

        fname = op.abspath( './' + refname + suffix + ext )
        return fname

