#!/usr/bin/env python
# coding: utf-8
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

"""
phantom.py - Nipype interface for generating phantom data in PySBR

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

from pysbr.tools.sbr_get_area_stats import findsecondaryspots3d

import json

class PhantomGenerateInputSpec(BaseInterfaceInputSpec):
    in_file = File( exists=True, mandatory=True,
                    desc='segmentation file, the aseg.mgz file from FreeSurfer' )
    grade = traits.List( traits.Int( 0, min=0, max=3),
                         desc='pathology grade to be mimicked, all grades by default' )
    out_prefix = traits.String( desc='set a prefix for outputs' )

class PhantomGenerateOutputSpec(TraitedSpec):
    out_grade = traits.List( traits.Int( 0, min=0, max=3),
                         desc='pathology grade to be mimicked, all grades by default' )

    out_files = traits.List( File( desc='list of generated phantoms, corresponding to\
                                        the grades list' ))

    out_rois = traits.List( File( desc='list of files containing activation areas' ))
    out_seg  = File( desc='filename of associated segmentation' )


class PhantomGenerate(BaseInterface):
    """ Template source interface gives easy access to template
        datasets. It also executes the template analysis to cache
        features in the first run
    """

    input_spec = PhantomGenerateInputSpec
    output_spec = PhantomGenerateOutputSpec

    def _run_interface(self, runtime):
        self._grades = range(0,4)

        if isdefined( self.inputs.grade ):
            self._grades = self.inputs.grade

        self._out_prefix = op.abspath( './dat_phantom' )
        if isdefined( self.inputs.out_prefix ):
            self._out_prefix = self.inputs.out_prefix

        self._out_files = [ '%s_grad%d.nii.gz' % (self._out_prefix,i) for i in self._grades ]
        self._out_rois = [ '%s_roi_%s.nii.gz' % (self._out_prefix,i) for i in [ 'c_l', 'p_l', 'c_r', 'p_r' ] ]
        self._out_seg = '%s_seg.nii.gz' % self._out_prefix


        self.generate_dat_phantom( self.inputs.in_file,
                                   self._grades )


        return runtime

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs['out_grade'] = self._grades
        outputs['out_files'] = self._out_files
        outputs['out_rois']  = self._out_rois
        outputs['out_seg']   = self._out_seg
        return outputs

    def generate_dat_phantom( self, in_seg, grades=[0,1,2,3] ):
        import os
        import os.path as op
        import nibabel as nib
        import numpy as np
        import scipy.ndimage as im
        
        def disk_structure(size=(2,2,2)):
            totalSize = np.array(size)*2+1
            struct = np.zeros(totalSize)
            coord = np.indices(totalSize).astype(np.float32)
            positions = coord.T - ((totalSize.astype(np.float32)-1)/2)
            terms = np.sum( ((positions**2)/(np.array(size)**2).astype(np.float32)).T, axis=0 )
            struct[terms<=1] = 1
            return struct.astype(np.bool)
            
        max_value = 255

        # Freesurfer's necessary labels
        left_caudate = 11
        left_putamen = 12
        right_caudate= 50
        right_putamen= 51
        left_cortex = 3
        right_cortex = 42
        left_wm = 2
        right_wm = 41
        brain_stem = 16
        c_callosum = 86
        csf = 24
        ventricle_3rd = 14
        ventricle_4th = 15
        ventricle_left = 4
        ventricle_right= 43
        left_cerebellum_cortex = 8
        right_cerebellum_cortex = 47

        wm_labels = [ left_wm, right_wm, brain_stem, c_callosum ]
        gm_labels = [ left_cortex, right_cortex, left_cerebellum_cortex, right_cerebellum_cortex ]
        no_wm = gm_labels + [ csf, ventricle_3rd, ventricle_4th, ventricle_left, ventricle_right, left_caudate, left_putamen, right_caudate, right_putamen ] 
        
        # Read images
        seg_im = nib.load( in_seg )
        seg_data = seg_im.get_data()

        msk = np.zeros_like( seg_data )
        msk[ seg_data>0.0] = 1

        hdr = seg_im.get_header()
        hdr['data_type'] = 2
        hdr.set_data_dtype( np.uint8 )

        # Generate wm, gm binary masks 
        seg_gm = np.zeros_like( seg_data )
        seg_gm[ np.any( [ seg_data==val for val in gm_labels ], axis=0) ] = 1
        seg_gm = seg_gm  * msk

        seg_wm = np.zeros_like( seg_data )
        seg_wm[ np.any( [ seg_data==val for val in no_wm ], axis=0 ) ] = 1
        seg_wm = msk - seg_wm
        seg_wm[ seg_wm<0 ] = 0
              
        # Generate putamen and caudate activity areas
        putamen_l = np.zeros_like( seg_data )
        putamen_r = np.zeros_like( seg_data )
        putamen_l[ seg_data==left_putamen ] = 1.0
        putamen_r[ seg_data==right_putamen ] = 1.0

        caudate_l = np.zeros_like( seg_data )
        caudate_r = np.zeros_like( seg_data )
        caudate_l[ seg_data==left_caudate ] = 1.0
        caudate_r[ seg_data==right_caudate ] = 1.0

        # Save original ROIs
        for data,fname in zip( [ caudate_l, putamen_l, caudate_r, putamen_r], self._out_rois ):
            nib.save( nib.Nifti1Image( data.astype( np.uint8 ), seg_im.get_affine(), hdr ), fname )

        putamen = putamen_l+putamen_r
        caudate = caudate_l + caudate_r
        
        dist_l = im.distance_transform_edt( (putamen_l*(-1.0)+1.0) )
        dist_r = im.distance_transform_edt( (putamen_r*(-1.0)+1.0) )
        distmap = (dist_l*putamen_r + dist_r*putamen_l)
       
        mindist = np.min( distmap[distmap>0].reshape(-1) )
        maxdist = np.max( distmap[distmap>0].reshape(-1) )

        caudatemap = np.zeros_like( distmap )
        caudatemap[ (distmap>0) & (distmap<= mindist+0.1*maxdist) ] = 1
        caudatemap = im.binary_dilation( caudatemap.astype(np.int) , structure=disk_structure((5,5,5)) ).astype(np.int)
        caudatemap_l = caudatemap * caudate_l
        caudatemap_r = caudatemap * caudate_r
        caudatemap_l = im.binary_dilation( caudatemap_l.astype(np.int) , structure=disk_structure((3,3,3)) ).astype(np.int)
        caudatemap_r = im.binary_dilation( caudatemap_r.astype(np.int) , structure=disk_structure((3,3,3)) ).astype(np.int)
        
        putamenmap_l = im.binary_dilation( putamen_l.astype(np.int) , structure=disk_structure((2,2,2)) ).astype(np.int)
        putamenmap_r = im.binary_dilation( putamen_r.astype(np.int) , structure=disk_structure((2,2,2)) ).astype(np.int)
        putamenmap_l[caudatemap_l!=0] = 0
        putamenmap_r[caudatemap_r!=0] = 0

        # Save segmentation
        rois = [ caudatemap_l, putamenmap_l, caudatemap_r, putamenmap_r ]
        seg_map = np.any( np.array(rois), axis=0 ).astype( np.uint8 )
        nib.save( nib.Nifti1Image( seg_map, seg_im.get_affine(), hdr ), self._out_seg )
        
        putamenmap_l = im.distance_transform_edt( putamenmap_l )
        putamenmap_r = im.distance_transform_edt( putamenmap_r )

        caudatemap_l = caudatemap_l/caudatemap_l.max()
        caudatemap_r = caudatemap_r/caudatemap_r.max()
        putamenmap_l = putamenmap_l/ putamenmap_l.reshape(-1).max()
        putamenmap_r = putamenmap_r/ putamenmap_r.reshape(-1).max()
        
        caudatemap = caudatemap_l+caudatemap_r
        putamenmap = putamenmap_l+putamenmap_r

        # Generate contrast
        basal_act = 10
        bg_data = np.random.rayleigh( scale=0.05*basal_act, size=seg_data.shape )
        gm_data = np.random.normal( loc=basal_act, scale=5.0, size=seg_data.shape )
        wm_data = np.random.normal( loc=0.5*basal_act, scale=10.0, size=seg_data.shape )

        seg_gm[ seg_map>0 ] = 0
        seg_wm[ seg_map>0 ] = 1
        
        # nib.save( nib.Nifti1Image( seg_wm.astype( np.float32 ), seg_im.get_affine(), hdr ), self._out_prefix + '_wm.nii.gz' )
        # nib.save( nib.Nifti1Image( seg_gm.astype( np.float32 ), seg_im.get_affine(), hdr ), self._out_prefix + '_gm.nii.gz' )

        brain_data = bg_data * (1.0-msk) + gm_data*seg_gm + wm_data*seg_wm
        brain_data[brain_data<0]=0.0

        #                 Grade 0        Grade I          Grade II       Grade III
        c_caudate =  [[ 1.00, 0.98 ], [ 1.00, 0.97 ], [ 0.87, 0.95 ], [ 0.97, 0.85 ] ]
        c_putamen =  [[ 1.00, 0.97 ], [ 0.90, 0.60 ], [ 0.35, 0.60 ], [ 0.15, 0.10 ] ]
        #c_caudate =  [[ 22, 22 ], [20, 15], [15,10], [3,3] ]
        #c_putamen =  [[ 20, 20 ], [15,  5], [ 5, 2], [1,0] ]
        
        
        hdr['data_type'] = 16
        hdr.set_data_dtype( np.float32 )
        norm = 20
        signal = np.random.normal( loc=2.0*basal_act, scale=4.0, size=seg_data.shape )
        nib.save( nib.Nifti1Image( im.gaussian_filter( brain_data, 4 ), seg_im.get_affine(), hdr ), self._out_prefix + '_brain.nii.gz' )

        for grade,fname in zip(grades,self._out_files):
            new_contrast = caudatemap_l * c_caudate[grade][0]* signal + \
                           putamenmap_l * c_putamen[grade][0]* signal + \
                           caudatemap_r * c_caudate[grade][1]* signal + \
                           putamenmap_r * c_putamen[grade][1]* signal + \
                           brain_data.astype(float)

            new_contrast[new_contrast<0] = 0
            new_contrast = im.gaussian_filter( new_contrast, 4 )
            new_contrast *= (max_value* 1.0)/ norm
            nib.save( nib.Nifti1Image( new_contrast.astype( np.float32 ), seg_im.get_affine(), hdr ),  fname )        
           
    

