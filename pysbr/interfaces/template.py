#!/usr/bin/env python
# coding: utf-8
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

"""
template.py - Nipype compliant interface for templates

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
from nipype.interfaces.base import (TraitedSpec, File, InputMultiPath,
                                    OutputMultiPath, Undefined, traits,
                                    isdefined, OutputMultiPath, 
                                    CommandLineInputSpec, CommandLine,
                                    BaseInterface, BaseInterfaceInputSpec,
                                    traits, Directory )
from nipype.utils.filemanip import split_filename,fname_presuffix

import csbr
from pysbr.tools.findsecspots import find_secondary_spots3d_intensity_weighted
import pysbr.tools.polar as polar

import itk
import json

class TemplateSourceInputSpec(BaseInterfaceInputSpec):
    base_dir = Directory( '../../data/template', usedefault=True, mandatory=True,
                          desc='directory where template is stored' )
    name = traits.Str( 'normal', mandatory=True, usedefault=True,
                       desc='template name' )
    force_init = traits.Bool( False, usedefault=True,
                              desc='Force initialization of template' )
    rseedratio = traits.Float( 0.6, usedefault=True,
                               desc='relative intensity ratio with respect to the maximum image intensity to search second spot.' ) # (-E)
    rmeanspotratio = traits.Float( 1.6, usedefault=True,
                                   desc='with respect to the mean region intensity to stop region growing.' ) # (-P)
    rbridgespotratio = traits.Float( 0.7, usedefault=True,
                                     desc='with respect to the minimum intensity of thebridge between the two seed points to stop region growing.' ) # (-B)
    uphill = traits.Float( 1.0, usedefault=True,
                           desc='factor that describes in region growing by which factor the intensity may increase to still add to the region, should be >= 1.0' )
    maxvolume = traits.Int( 300000, usedefault=True, desc='Maximum volume a hot spot region can assume during segmentation.' ) # (-V)
    max_sensible_volume = traits.Int( 30000, usedefault=True, desc='Maximum activation volume to accept for a segmentation.')
    minvolume = traits.Int( 5, usedefault=True, desc='Minimum volume for proceeding with secondary spots detection.' )
    evratio = traits.Float( 1.75, usedefault=True,
                            desc='ratio multiplier to decide how many secondary landmarks are used based on the ratio of the first and second PCA eigenvalue') 

class TemplateSourceOutputSpec(TraitedSpec):
    out_file = File( desc='filename of template data file' )
    out_spots= File( desc='textfile containing landmarks' )
    out_bbox  = File( desc='filename of associated bounding box' )
    out_seg = File( desc='filename of associated segmentation' )
    out_rois = traits.List( File(), desc='filename of associated ground-truth ROIs' )
    out_spots_im = File( desc='filename of image with landmarks' )


class TemplateSource(BaseInterface):
    """ Template source interface gives easy access to template
        datasets. It also executes the template analysis to cache
        features in the first run
    """

    input_spec = TemplateSourceInputSpec
    output_spec = TemplateSourceOutputSpec

    def _run_interface(self, runtime):
        tmpl_dir = op.abspath( self.inputs.base_dir )
        tname = self.inputs.name

        self._outfnames = [ op.join( tmpl_dir, tname, i ) for i in ['template.nii.gz', 'landmarks.txt', 'bbox.nii.gz', 'segmentation.nii.gz' ] ]
        self._out_rois = [ op.join( tmpl_dir, 'anat', '%s.nii.gz' % i ) for i in [ 'caudate_l', 'putamen_l', 'caudate_r', 'putamen_r' ] ]
        check_files = [ op.exists( fname ) for fname in self._outfnames+self._out_rois ]

        # If template is not initialized, do it.
        if self.inputs.force_init or not np.all( check_files ):
            self._generate_template( tmpl_dir, tname )

        return runtime

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs['out_file']  = self._outfnames[0]
        outputs['out_spots'] = self._outfnames[1]
        outputs['out_bbox']  = self._outfnames[2]
        outputs['out_seg']   = self._outfnames[3]
        outputs['out_rois']  = self._out_rois
        return outputs

    def _generate_template( self, tmpl_dir, tname, lm_fixed_steps=True, smooth_mask=True ):
        import scipy.ndimage as im
        import nibabel as nib
        import numpy as np
        import os.path as op
        import sbr
        from math import sin, cos, asin, acos, atan
        from scipy.interpolate import griddata
        from scipy.stats import scoreatpercentile

        # check inputs
        tplfname = self._outfnames[0]
        labfname = op.join( tmpl_dir, tname, 'segmentation.nii.gz' )
        mskfname = op.join( tmpl_dir, 'anat', 'brainmask.nii.gz' )

        if not ( op.exists( tplfname ) and op.exists( mskfname ) ):
            raise RuntimeError( 'A template cannot be generated if either template data or corresponding brainmask are not present' )
       

        # Read template in ITK
        imread = itk.ImageFileReader[itk.Image[itk.F,3]].New()
        imread.SetFileName( tplfname )
        imread.Update()
        itktpl = imread.GetOutput()
        IndexType = itk.ContinuousIndex[itk.D, 3]

        im_tpl = nib.load( tplfname )
        hdr = im_tpl.get_header().copy()
        hdr['data_type'] = 2
        hdr.set_data_dtype( np.uint8 )

        tpl_data = im_tpl.get_data()
        if op.exists( labfname ):
            labim = nib.load( labfname )
            labdata = labim.get_data()
        else:
            # This code generates a preliminary segmentation to generate the bbox
            # It assumes that activity centers are rather clear, and separated at
            # the 99.5 percentile. Then it chooses the objects closer to the center
            # of the image.
            labdata = np.zeros_like( tpl_data )
            bgval = scoreatpercentile( tpl_data, 15 )
            val = scoreatpercentile( tpl_data[tpl_data>bgval], 99.5 )
            labdata[tpl_data>val] = 1
            label_im, nb_labels = im.label( labdata )
            if nb_labels<=0:
                raise RuntimeError( 'No activation areas detected' )
            elif nb_labels<=2:
                labdata[label_im>0] = 1
            elif nb_labels>2:
                centers = im.center_of_mass( labdata, label_im, range(1,nb_labels+1) )
                c = np.array( labdata.shape )*0.5
                dists = [ np.linalg.norm( v-c ) for v in centers ]
                min1 = np.min( dists )
                labdata = np.zeros_like( tpl_data )
                labdata[ label_im == (np.argwhere( dists==min1 ) + 1) ] = 1
                dists2 = np.ma.masked_equal( dists, min1 )
                min2 = dists2.min()
                labdata[ label_im == (np.argwhere( dists==min2 ) + 1) ] = 1

            labim = nib.Nifti1Image( labdata.astype(np.uint8), im_tpl.get_affine(), hdr )
            nib.save( labim, labfname )

        # Create roi mask (find bounding box)
        mask = np.zeros_like( labdata )

        indices = np.array( np.nonzero(labdata) )
        first = indices.min(1)
        last = indices.max(1) + 1

        mask[first[0]:last[0],first[1]:last[1],first[2]:last[2]] = 1

        roi = im.binary_dilation( mask, np.ones( (10,10,10) ) ).astype( np.uint8 )
        roinii = nib.Nifti1Image( roi, labim.get_affine(), hdr )
        nib.save( roinii, self._outfnames[2] )

        # Find spots by standard procedure
        im_data = im_tpl.get_data() * roi

        (labels, nspots, hotvox, volumes) = csbr.findspots3d(im_data, 
                                                 rseedlimit=self.inputs.rseedratio,
                                                 rspotmeanlimit=self.inputs.rmeanspotratio, 
                                                 rspotbridgelimit=self.inputs.rbridgespotratio, 
                                                 maxvolume=self.inputs.maxvolume)
        hotvox = adjust_spot_positions(image, labels, hotvox)
        
        vol_check = np.all( [v>0 for v in volumes] )

        if not vol_check:
            raise RuntimeError( 'Findspots did not succeed on detecting hotspots in template' )

        # Find secondary spots
        landmark_map = find_secondary_spots3d_intensity_weighted( im_data, labels, hotvox, self.inputs.minvolume, self.inputs.evratio, lm_fixed_steps, adjust=True )

        # translate the landmarks to physical coordinates 
        for k,v in landmark_map.iteritems():
            v = np.atleast_2d( v )
            landmark_map[k] = [ tuple( itktpl.TransformContinuousIndexToPhysicalPoint( IndexType( [ float(val) for val in vox ] ) ) ) for vox in v ]

        # Find segmentations' pointcloud
        for l in range(1,nspots+1):
            seg = np.zeros_like( labels )
            seg[labels==l] = 1
            seg_er = im.binary_erosion( seg )
            boundary = seg - seg_er
            idxs = np.argwhere( boundary==1 )
            landmark_map['surf%d'%l] = [ tuple( itktpl.TransformContinuousIndexToPhysicalPoint( IndexType( [ float(val) for val in vox ] ) ) ) for vox in idxs ]

        # Find fixed (pial) spots in brainmask
        bmaskimg = nib.load( mskfname )
        imshape = np.array( bmaskimg.get_shape() )
        bmask = bmaskimg.get_data()
        bmask[bmask>0] = 1.0

        # Smooth and preprocess mask
        if smooth_mask:
            size = (10,10,10)
            struc = np.zeros( size )
            c = ( np.array( size ) -1 )*0.5
            positions = np.argwhere( struc==0 )
            dist = np.array( [ np.linalg.norm(val-c)< 4.0 for val in positions ] )
            struc[ dist.reshape( size ) ] = 1
            bmask = im.binary_closing( bmask, structure=struc ).astype( np.uint8 )
            bmask = im.binary_fill_holes( bmask ).astype( np.uint8 )
            bmask = im.gaussian_filter( bmask.astype(np.float32), 8 ).astype( np.float32 )
            bmask[bmask>=0.3] = 1
            bmask[bmask<1] = 0

       
        bmask_eroded = im.binary_erosion( bmask ).astype( np.uint8 )
        boundary = bmask - bmask_eroded

        ntheta = 9
        nphi = 18

        # Read brainmask in ITK
        imread2 = itk.ImageFileReader[itk.Image[itk.F,3]].New()
        imread2.SetFileName( mskfname )
        imread2.Update()
        itkmsk = imread2.GetOutput()

        cidx = im.center_of_mass( bmask )
        cmass = tuple( itkmsk.TransformContinuousIndexToPhysicalPoint( IndexType( [ float(val) for val in cidx ] ) )) 
        Z = np.array( np.argwhere( boundary>0 ), dtype=np.uint8 )
        points = np.array( [ tuple( itkmsk.TransformContinuousIndexToPhysicalPoint( IndexType( [ float(val) for val in vox ] ) )) for vox in Z ] )

        X,Y,S = polar.polar_rep( points, cmass, feat=Z, res=[ ntheta, nphi ] )
        
        coords_list = [ tuple( itkmsk.TransformContinuousIndexToPhysicalPoint( IndexType( [ float(val) for val in vox ] ) )) for vox in S.reshape(-1, 3) ]
        landmark_map['static'] = coords_list

        # one should add the hotspot index and the landmark index to the output to be able to set the
        # proper correspondences 
        with open( self._outfnames[1], 'w') as f:
            json.dump(landmark_map, f)

        nib.save( nib.Nifti1Image( labels, labim.get_affine(), hdr ), op.join( tmpl_dir, tname, 'labels.nii.gz' ) )

        spotsimdata = np.zeros_like( bmask )
        lid = 1
        for k,v in landmark_map.iteritems():
            for p in v:
                idx = itkmsk.TransformPhysicalPointToIndex( [ float(val) for val in p ] ) 
                spotsimdata[tuple(idx)] = lid
            lid+=1

        hdr = bmaskimg.get_header().copy()
        hdr['data_type'] = 2
        hdr.set_data_dtype( np.uint8 )
        nib.save( nib.Nifti1Image( spotsimdata, bmaskimg.get_affine(), hdr ), op.join( tmpl_dir, tname, 'spotsim.nii.gz' ) )

