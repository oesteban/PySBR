#!/usr/bin/env python
# coding: utf-8
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

"""
transform.py

Nipype compliant interface for PySBR transforms

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
import json
import itk



class TransformBaseInterfaceInputSpec( BaseInterfaceInputSpec ):
    in_initialization = traits.List( File( exists=True ), mandatory=True,
                           desc='list of files with initialization transforms' )

    in_init_flags = traits.List( traits.Bool(False),mandatory=False,
                                 desc='list of corresponding flags to invert transforms')

class TransformBaseInterface(BaseInterface):
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

    def _read_tforms( self, transforms, flags=None ):
        transforms = np.atleast_1d( transforms )
        if flags is None:
            flags = [ False for n in transforms ]
        flags = np.atleast_1d( flags )

        inittf = itk.CompositeTransform[ itk.D, 3].New()
        tfs = []
        for tfname,flag in zip(transforms, flags ):
            # Read transform file
            tf_read = itk.TransformFileReader.New()
            tf_read.SetFileName( tfname )
            tf_read.Update()

            nXforms = tf_read.GetTransformList().size()

            print 'Reading %s transform file, contains #%d transforms, invert=%s' % (tfname,nXforms,flag)

            for tf_id in range(nXforms):
                params = (tf_read.GetTransformList()[tf_id]).GetParameters()
                if params.Size() == 12:
                    tf = itk.MatrixOffsetTransformBase[itk.D,3,3].New()
                    tf.SetParameters( params )
                elif params.Size() == 6:
                    tf = itk.Euler3DTransform[itk.D].New()
                    tf.SetParameters( params )
                elif params.Size() == 3:
                    tf = itk.TranslationTransform[itk.D,3].New()
                    tf.SetParameters( params )
                else:
                    raise NotImplementedError( 'Trying to use a non-linear transform with %d parameters' % params.Size() )
                
                if flag:
                    inittf.AddTransform( tf.GetInverseTransform() )
                else:
                    inittf.AddTransform( tf )

        return inittf



    def _read_tform_file( self, tfname ):
        tfs = []
        # Read transform file
        tf_read = itk.TransformFileReader.New()
        tf_read.SetFileName( tfname )
        tf_read.Update()
        for tf_id in range(0,tf_read.GetTransformList().size() ):
            params = (tf_read.GetTransformList()[tf_id]).GetParameters()
            if params.Size() == 12:
                tf = itk.MatrixOffsetTransformBase[itk.D,3,3].New()
                tf.SetParameters( params )
            elif params.Size() == 6:
                tf = itk.Euler3DTransform[itk.D].New()
                tf.SetParameters( params )
            elif params.Size() == 3:
                tf = itk.TranslationTransform[itk.D,3].New()
                tf.SetParameters( params )
            else:
                raise NotImplementedError( 'Trying to use a non-linear transform with %d parameters' % params.Size() )
            tfs.append( tf )
        return tfs



class EBSTransformInputSpec(TransformBaseInterfaceInputSpec):
    in_target = File( exists=True, mandatory=True,
                      desc="reference/target image" )
    in_source = File( exists=True, mandatory=True,
                      desc="moving/source image" )
    in_transform = File( exists=True, mandatory=True,
                         desc="file containing displacements" )
    in_apply = traits.List( File( exists=True), mandatory=False,
                            desc='list of images in source space to apply the transform\
                                  too' )

    stiffness = traits.Float( 0.001, mandatory=True, usedefault=True,
                              desc='Stiffness of the spline. A stiffness of 0.0 results\
                                    in the standard interpolating spline. A non-zero stiffness\
                                    allows the spline to approximate rather than interpolate the\
                                    landmarks. Stiffness values are usually rather small,\
                                    typically in the range of 0.001 to 0.1. The approximating\
                                    spline formulation is based on the short paper by R. Sprengel,\
                                    K. Rohr, H. Stiehl. "Thin-Plate Spline Approximation for Image\
                                    Registration". In 18th International Conference of the IEEE\
                                    Engineering in Medicine and Biology Society. 1996.')

    alpha = traits.Float( 10.0, mandatory=True, usedefault=True,
                          desc='Set alpha. Alpha is related to Poisson\'s Ratio (\\nu) as\
                                \\alpha = 12 (1-\\nu)-1. \\nu_{gold} ~ 0.43; \\nu_{foam} ~ 0.30' ) 
    out_moving = File( desc='filename for the source interpolated into \
                            target\'s space' )
    out_init = File( desc='filename for the source interpolated into \
                           target\'s space (only initialization)' )
    out_transform = File( desc='filename containing the bspline transform' )

    out_apply = traits.List( File(), desc='output files with transforms applied' )


class EBSTransformOutputSpec(TraitedSpec):
    out_moving = File( desc='filename for the source interpolated into \
                            target\'s space' )
    out_transform = File( desc='filename containing the bspline transform' )
    out_init = File( desc='filename for the source interpolated into \
                           target\'s space (only initialization)' )
    out_apply = traits.List( File(), desc='output files with transforms applied' )

class EBSTransform(TransformBaseInterface):
    input_spec = EBSTransformInputSpec
    output_spec = EBSTransformOutputSpec

    def _run_interface(self, runtime):
        import itk

        with open(self.inputs.in_transform) as f:
            tfm_dict = json.load(f)

        if not isdefined( self.inputs.in_init_flags ):
            flags = None,
        else:
            flags = self.inputs.in_init_flags

        inittf = self._read_tforms( self.inputs.in_initialization,
                                    flags )

        # Read images
        image_t = itk.Image[itk.F,3]
        target_reader = itk.ImageFileReader[image_t].New()
        target_reader.SetFileName( self.inputs.in_target )
        target_reader.Update()
        target_im = target_reader.GetOutput()

        # Resample image (only initialization)
        if not isdefined( self.inputs.out_init ):
            self.inputs.out_init = self._gen_fname( self.inputs.in_source,
                                                      suffix='init' )
        source_reader = itk.ImageFileReader[image_t].New()
        source_reader.SetFileName( self.inputs.in_source )
        source_reader.Update()
        source_im = source_reader.GetOutput()
        resampler2 = itk.ResampleImageFilter[image_t,image_t].New()
        resampler2.SetInput( source_im )
        resampler2.SetReferenceImage( target_im )
        resampler2.UseReferenceImageOn()
        resampler2.SetTransform( inittf )
        resampler2.Update()
        writer2 = itk.ImageFileWriter[image_t].New()
        writer2.SetInput( resampler2.GetOutput() )
        writer2.SetFileName( self.inputs.out_init )
        writer2.Update()

        ebs = itk.ElasticBodySplineKernelTransform[itk.D,3].New()
        target_points = ebs.GetModifiableTargetLandmarks()
        source_points = ebs.GetModifiableSourceLandmarks()
        ebs.SetStiffness( self.inputs.stiffness )
        ebs.SetAlpha( self.inputs.alpha )

        # Debug code
        trgtim = nib.load( self.inputs.in_target )
        srcpoints = np.zeros_like( trgtim.get_data() )
        trgpoints = np.zeros_like( trgtim.get_data() )
   
        # Set points in transform
        idx = 0
        nStatic = len( tfm_dict['static'] )

        for pmap in tfm_dict['moving']:
            for sp,tp in zip( pmap['source'], pmap['target'] ):
                source_points.GetPoints().InsertElement( idx, [ float(val) for val in sp ] )
                target_points.GetPoints().InsertElement( idx, [ float(val) for val in tp ] )
                idx+=1
                
                trgpoints[ tuple( target_im.TransformPhysicalPointToIndex( [ float(v) for v in tp ] ) ) ] = nStatic + 10
                srcpoints[ tuple( target_im.TransformPhysicalPointToIndex( [ float(v) for v in sp ] ) ) ] = nStatic + 10

        # Set zero displacement for static landmarks
        for i,p in enumerate( tfm_dict['static'] ):
            target_points.GetPoints().InsertElement( idx, [ float(val) for val in p ] )
            source_points.GetPoints().InsertElement( idx, [ float(val) for val in p ] )
            idx+=1
            
            trgpoints[ tuple( target_im.TransformPhysicalPointToIndex( [ float(v) for v in p ] ) ) ] = idx
            srcpoints[ tuple( target_im.TransformPhysicalPointToIndex( [ float(v) for v in p ] ) ) ] = idx

        temp_path = os.getenv("PYSBR_TEMP")
        if temp_path is not None:
            nib.save( nib.Nifti1Image( trgpoints, trgtim.get_affine(), trgtim.get_header() ), 
                      os.path.join(temp_path, 'trgpoints.nii.gz' ))
            nib.save( nib.Nifti1Image( srcpoints, trgtim.get_affine(), trgtim.get_header() ), 
                      os.path.join(temp_path, 'srcpoints.nii.gz' ))


        # Apply EBS
        ebs.ComputeWMatrix() # Hold on, this is the magic

        # Write transform
        if not isdefined( self.inputs.out_transform ):
            self.inputs.out_transform = self._gen_fname( self.inputs.in_source,
                                                      suffix='ebs', ext='.tfm' )
        tf_writer = itk.TransformFileWriter.New()
        tf_writer.SetInput( ebs )
        tf_writer.SetFileName( self.inputs.out_transform )
        tf_writer.Update()

        # Composite transform applies the last added transform first (LIFO)
        compositetf = itk.CompositeTransform[ itk.D, 3].New()
        compositetf.AddTransform( ebs )
        compositetf.AddTransform(inittf)

        # Resample images
        if not isdefined( self.inputs.out_moving ):
            self.inputs.out_moving = self._gen_fname( self.inputs.in_source,
                                                      suffix='ebs' )

        in_files = [ self.inputs.in_source ]
        out_files = [ self.inputs.out_moving ]

        if isdefined( self.inputs.in_apply ) and (len( self.inputs.in_apply )>0):
            in_files+=self.inputs.in_apply

            if isdefined( self.inputs.out_apply ) and ( len( self.inputs.in_apply ) == len(self.inputs.out_apply )):
                out_files+=self.inputs.out_apply
            else:
                self.inputs.out_apply = []
                for i,f in enumerate(self.inputs.in_apply):
                    new_name = self._gen_fname( self.inputs.in_source, suffix=('apply%02d'%i) )
                    self.inputs.out_apply+= [new_name]
                    out_files+= [new_name]
                
        for fin,fout in zip(in_files,out_files):
            source_reader = itk.ImageFileReader[image_t].New()
            source_reader.SetFileName( fin )
            source_reader.Update()
            resampler = itk.ResampleImageFilter[image_t,image_t].New()
            resampler.SetInput( source_reader.GetOutput() )
            resampler.SetReferenceImage( target_im )
            resampler.UseReferenceImageOn()
            resampler.SetTransform( compositetf )
            resampler.Update()
            writer = itk.ImageFileWriter[image_t].New()
            writer.SetInput( resampler.GetOutput() )
            writer.SetFileName( fout )
            writer.Update()

        return runtime

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs['out_moving'] = self.inputs.out_moving
        outputs['out_init'] = self.inputs.out_init
        outputs['out_transform'] = self.inputs.out_transform
        if isdefined( self.inputs.in_apply):
            outputs['out_apply'] = self.inputs.out_apply
        return outputs


class TransformPointsInputSpec(TransformBaseInterfaceInputSpec):
    in_target = File( exists=True, mandatory=True,
                      desc="reference/target image" )
    in_source = File( exists=True, mandatory=True,
                      desc="moving/source image" )
    in_points = File( exists=True, mandatory=True,
                      desc="file containing source points" )
    out_points = File( desc='filename of written output that contains \
                             the points from in_points transformed to \
                             target coordinates' )
    #inverse = traits.Bool( False, mandatory=True, usedefault=True,
    #                       desc='perform inverse transform' )


class TransformPointsOutputSpec(TraitedSpec):
    out_points = File( desc='filename for the source interpolated into \
                            target\'s space' )

class TransformPoints(TransformBaseInterface):
    """ TransformpPoints is an interface that applies an affine transformation
        to points contained in an input file
    """
    input_spec = TransformPointsInputSpec
    output_spec = TransformPointsOutputSpec

    def _run_interface(self, runtime):
        with open(self.inputs.in_points) as f:
            points_dict = json.load(f)

        if not isdefined( self.inputs.in_init_flags ):
            flags = None,
        else:
            flags = self.inputs.in_init_flags

        inittf = self._read_tforms( self.inputs.in_initialization, flags )

        # Transform points
        out_points = {}
        for k,v in points_dict.iteritems():
            out_group = []
            for i,p in enumerate( v ):
                p_prime = inittf.TransformPoint( [ float(val) for val in p ] )
                out_group.append( [ val for val in p_prime ] )
            out_points[k] = out_group

        if not isdefined( self.inputs.out_points ):
            self.inputs.out_points = self._gen_fname( self.inputs.in_target,
                                                      suffix='points_ref', ext='.txt' )

        out_fname = self.inputs.out_points
        with open( out_fname, 'w') as f:
            json.dump(out_points, f)

        return runtime

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs["out_points"] = self.inputs.out_points
        return outputs


#class ApplyEBSTransformInputSpec(BaseInterfaceInputSpec):
#    in_target = File( exists=True, mandatory=True,
#                      desc="reference/target image" )
#    in_source = File( exists=True, mandatory=True,
#                      desc="moving/source image" )
#    in_initialization = traits.List( File(exists=True), mandatory=True,
#                           desc='filename list containing the input linear initialization transforms' )
#    in_transform = File(exists=True, mandatory=True,
#                        desc='filename list containing the input EBS transform' )
#    out_moving = File( desc='output filename' )
#
#
#class ApplyEBSTransformOutputSpec(TraitedSpec):
#    out_moving = File( desc='filename for the source interpolated into \
#                            target\'s space' )
#
#class ApplyEBSTransform(BaseInterface):
#    input_spec = ApplyEBSTransformInputSpec
#    output_spec = ApplyEBSTransformOutputSpec
#
#    def _run_interface(self, runtime):
#        import itk
#        # Read images
#        image_t = itk.Image[itk.F,3]
#        target_reader = itk.ImageFileReader[image_t].New()
#        target_reader.SetFileName( self.inputs.in_target )
#        target_reader.Update()
#        target_im = target_reader.GetOutput()
#
#        source_reader = itk.ImageFileReader[image_t].New()
#        source_reader.SetFileName( self.inputs.in_source )
#        source_reader.Update()
#        source_im = source_reader.GetOutput()
#       
#        # Read transforms 
#        tf_list = []
#        for tfname in np.atleast_1d( self.inputs.in_initialization ):
#            # Read transform file
#            tf_read = itk.TransformFileReader.New()
#            tf_read.SetFileName( tfname )
#            tf_read.Update()
#            for tf_id in range(0,tf_read.GetTransformList().size() ):
#                params = (tf_read.GetTransformList()[tf_id]).GetParameters()
#                if params.Size() == 12:
#                    tf = itk.MatrixOffsetTransformBase[itk.D,3,3].New()
#                    tf.SetParameters( params )
#                elif params.Size() == 6:
#                    tf = itk.Euler3DTransform[itk.D].New()
#                    tf.SetParameters( params )
#
#                tf_list.append( tf )
#
#        # Composite transform applies the last added transform first (LIFO)
#        compositetf = itk.CompositeTransform[ itk.D, 3].New()
#        for tf in tf_list.reverse():
#            compositetf.Add( tf )
#
#        # Resample image (affine)
#        resampler = itk.ResampleImageFilter[image_t,image_t].New()
#        resampler.SetInput( source_im )
#        resampler.SetReferenceImage( target_im )
#        resampler.UseReferenceImageOn()
#        resampler.SetTransform( compositetf )
#        resampler.Update()
#        affine_im = resampler.GetOutput()
#
#        # EBS transform
#        tf_read = itk.TransformFileReader.New()
#        tf_read.SetFileName( selt.inputs.in_transform )
#        tf_read.Update()
#        ebs = itk.ElasticBodySplineKernelTransform[itk.D,3].New()
#        ebs.SetParameters( tf_read.GetTransformList()[0].GetParameters() ) # read parameters from file
#
#        # Resample image (EBS)
#        if not isdefined( self.inputs.out_moving ):
#            self.inputs.out_moving = self._gen_fname( self.inputs.in_source,
#                                                      suffix='moving' )
#        resampler = itk.ResampleImageFilter[image_t,image_t].New()
#        resampler.SetInput( affine_im )
#        resampler.SetReferenceImage( target_im )
#        resampler.UseReferenceImageOn()
#        resampler.SetTransform( ebs )
#        resampler.Update()
#        output_im = resampler.GetOutput()
#        writer = itk.ImageFileWriter[image_t].New()
#        writer.SetInput( output_im )
#        writer.SetFileName( self.inputs.out_moving )
#        writer.Update()
#        return runtime
#
#    def _list_outputs(self):
#        outputs = self._outputs().get()
#        outputs['out_moving'] = self.inputs.out_moving
#        return outputs
#
#    def _gen_fname( self, refname, suffix=None, ext=None ):
#        basename = op.basename( refname )
#        
#        refname, refext = op.splitext( basename )
#        if refext =='.gz':
#            refname, refext2 = op.splitext(refname)
#            refext = refext2 + refext
#
#        if ext is None:
#            ext = refext
#
#        if suffix is None:
#            suffix = ''
#        else:
#            suffix = '_'+suffix
#
#        fname = op.abspath( './' + refname + suffix + ext )
#        return fname


