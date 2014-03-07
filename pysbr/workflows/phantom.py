#!/usr/bin/env python
# coding: utf-8
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

"""
phantom.py - The PySBR phantom generation workflow

Copyright (c) 2013, code@oscaresteban.es (Oscar Esteban), and 
                    gw.fossdev@gmail.com (Gert Wollny)
                    with Biomedical Image Technology, UPM (BIT-UPM)
All rights reserved.
This file is part of PySBR.

"""

__author__ = "Oscar Esteban, Gert Wollny"
__copyright__ = "Copyright 2013, Biomedical Image Technologies (BIT), \
                 Universidad Politécnica de Madrid"
__credits__ = ["Oscar Esteban", "Gert Wollny"]
__license__ = "FreeBSD"
__version__ = "0.1"
__maintainer__ = "Oscar Esteban"
__email__ = "code@oscaresteban.es"
__status__ = "Prototype"

import sys
import os.path as op
import os
import glob

root_dir = op.abspath( '../../' )
sys.path.append( op.join( root_dir, 'pysbr') )

import interfaces as pysbr



import nipype.pipeline.engine as pe
import nipype.interfaces.freesurfer as fs
import nipype.interfaces.fsl.preprocess as fsl
import nipype.interfaces.utility as niu
import nipype.interfaces.io as nio

def to_orig( name='to_orig' ):
    to_orig = pe.Workflow( name=name )
    inputnode = pe.Node( niu.IdentityInterface( fields=['in_file', 'in_orig' ] ), name='inputnode' )
    outputnode = pe.Node( niu.IdentityInterface( fields=['out_file', 'out_orig' ] ), name='outputnode' )
    reorient = pe.Node( fs.MRIConvert( out_orientation='LPS'), name='toLPS' )
    applyreg = pe.Node( fs.ApplyVolTransform(reg_header=True, interp='nearest'), name='applyreg' )
    convert = pe.Node( fs.MRIConvert( out_type='niigz' ), name='convert' )
    to_orig.connect([
                      ( inputnode, reorient, [ ( 'in_orig', 'in_file' ) ])
                     ,( reorient,  applyreg, [ ( 'out_file', 'target_file' ) ] )
                     ,( inputnode, applyreg, [ ( 'in_file', 'source_file' ) ] )
                     ,( applyreg,   convert, [ ( 'transformed_file','in_file' ) ] )
                     ,( convert,   outputnode, [ ( 'out_file', 'out_file' ) ])
                     ,( reorient,   outputnode, [ ( 'out_file', 'out_orig' ) ])
                     ])
    
    return to_orig

def gen_phantom_workflow( name="GeneratePhantom", pixsize=(2.56,2.56,2.56) ):
    gen = pe.Workflow( name=name )
    inputnode = pe.Node( niu.IdentityInterface( fields=['in_file', 'in_mask', 'in_aseg' ] ), name='inputnode' )
    outputnode = pe.Node( niu.IdentityInterface( fields=['out_files','out_rois','out_seg', 'out_mri' ] ), name='outputnode' )
    seg_to_orig = to_orig(name='seg_to_orig')
    gen_data = pe.Node( pysbr.PhantomGenerate(), name='Gen' )
    resample0 = pe.MapNode( fs.MRIConvert(vox_size=pixsize, out_datatype='float'), name="resample0", iterfield=['in_file'] )
    resample1 = pe.MapNode( fs.MRIConvert(vox_size=pixsize, out_datatype='float'), name="resample1", iterfield=['in_file'] )
    resample2 = pe.Node   ( fs.MRIConvert(vox_size=pixsize, out_datatype='float'), name="resample2" )
    
    gen.connect([
                  (  inputnode, seg_to_orig, [ ('in_aseg', 'inputnode.in_file' ), ('in_file', 'inputnode.in_orig' )] )
                 ,(seg_to_orig,    gen_data, [ ('outputnode.out_file', 'in_file' )])
                 ,(   gen_data,   resample0, [ ('out_files','in_file' )])
                 ,(   gen_data,   resample1, [ ('out_rois','in_file' )])
                 ,(   gen_data,   resample2, [ ('out_seg','in_file' )])
                 ,(  resample0,  outputnode, [ ('out_file', 'out_files' ) ])
                 ,(  resample1,  outputnode, [ ('out_file', 'out_rois' ) ])
                 ,(  resample2,  outputnode, [ ('out_file', 'out_seg' ) ])
                 ,(seg_to_orig,  outputnode, [ ('outputnode.out_orig', 'out_mri' ) ])
                ])
    return gen

def generate_all( subjects_list, results_dir, ixi_dir, name='PhantomGeneration' ):
    processall = pe.Workflow( name=name )

    def _basename( fnames ):
        import numpy as np
        import os.path as op
        fnames = np.atleast_1d( fnames )
        
        res = [ op.basename( n ) for n in fnames ]
        print res
        return res

    def copyfile( in_files, subject_id, out_dir ):
        import os.path as op
        from os import makedirs
        from shutil import copyfile as copy
        from numpy import atleast_1d

        path = op.join( out_dir, subject_id )

        if not op.exists( path ):
            makedirs( path )

        out_files = []
        for in_file in atleast_1d( in_files ):
            out_file = op.join( path, op.basename( in_file ) )
            out_files+= [ out_file ]
            print '%s -> %s' % (in_file, out_file )
            copy( in_file, out_file )
        return out_files
            
    
    #initiate the infosource node
    infosource = pe.Node(niu.IdentityInterface(fields=['subject_id']),
                         name="inputnode" )
    
    #define the list of subjects your pipeline should be executed on
    infosource.iterables = ('subject_id', subjects_list )
    
    datasource = pe.Node( nio.DataGrabber(sort_filelist=False,infields=['subject_id'], outfields=['t1', 'aseg', 'brainmask']), name="dsource" )
    datasource.inputs.template = ixi_dir + '%s/mri/%s.mgz'

    datasink0 = pe.Node(niu.Function( input_names=['in_files','subject_id','out_dir' ], output_names=['out_files'], function=copyfile ), name='dsink0' )
    datasink0.inputs.out_dir=results_dir 
    datasink1 = pe.Node(niu.Function( input_names=['in_files','subject_id','out_dir' ], output_names=['out_files'], function=copyfile ), name='dsink1' )
    datasink1.inputs.out_dir=results_dir 
    datasink2 = pe.Node(niu.Function( input_names=['in_files','subject_id','out_dir' ], output_names=['out_files'], function=copyfile ), name='dsink2' )
    datasink2.inputs.out_dir=results_dir

    datasink3 = pe.Node(niu.Function( input_names=['in_files','subject_id','out_dir' ], output_names=['out_files'], function=copyfile ), name='dsink3' )
    datasink3.inputs.out_dir=results_dir 

    
    info = dict(t1=[['subject_id', 'T1' ]],
                aseg=[['subject_id','aseg'  ]],
                brainmask=[['subject_id', 'brainmask' ]] )
    
    datasource.inputs.template_args = info
    wf = gen_phantom_workflow()
    outputnode = pe.Node( niu.IdentityInterface( fields=['out_files','out_rois','out_seg','out_mri'] ), name='outputnode' )
    
    processall.connect([
                         ( infosource, datasource, [ ( 'subject_id', 'subject_id' ) ])
                        ,( datasource, wf, [ ('t1', 'inputnode.in_file'),
                                             ('aseg','inputnode.in_aseg') ])
                        ,( infosource, datasink0, [ ('subject_id', 'subject_id' ) ] )
                        ,( infosource, datasink1, [ ('subject_id', 'subject_id' ) ] )
                        ,( infosource, datasink2, [ ('subject_id', 'subject_id' ) ] )
                        ,( infosource, datasink3, [ ('subject_id', 'subject_id' ) ] )
                        ,( wf, datasink0, [ ('outputnode.out_files', 'in_files' )])
                        ,( wf, datasink1, [ ('outputnode.out_rois', 'in_files' )])
                        ,( wf, datasink2, [ ('outputnode.out_seg', 'in_files' )])
                        ,( wf, datasink3, [ ('outputnode.out_mri', 'in_files' )])
                        ,( datasink0, outputnode, [ ('out_files','out_files' )])
                        ,( datasink1, outputnode, [ ('out_files','out_rois' )])
                        ,( datasink2, outputnode, [ ('out_files','out_seg' )])
                        ,( datasink3, outputnode, [ ('out_files','out_mri' )])
                        ])
    return processall

if __name__== '__main__':
    ixi_dir = '/media/mnemea/IXI_dataset/FREESURFER/'
    subjects_list = [ op.basename( sub ) for sub in glob.glob( op.join( ixi_dir, 'S*' ) ) ]
    data_dir = op.abspath( '/media/mnemea/MINDt-Quantidopa/PySBR-PhantomImages/' )
    tmpl_dir = op.join( pysbr.get_datapath(), 'template' )
    work_dir = op.join( root_dir, 'temp')
    results_dir = op.join( data_dir, 'Phantoms' )
    
    if not op.exists( work_dir ):
        os.makedirs( work_dir )
        
    if not op.exists( results_dir ):
        os.makedirs( results_dir )

    wf = generate_all( subjects_list[0:60], results_dir, ixi_dir )
    wf.base_dir=work_dir
    wf.write_graph( format='pdf' )
#    wf.run()

