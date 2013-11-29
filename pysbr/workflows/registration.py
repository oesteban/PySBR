#!/usr/bin/env python
# coding: utf-8
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

"""
registration.py - The PySBR registration workflow

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

import nipype.pipeline.engine as pe
import nipype.interfaces.utility as niu
import os
import sysconfig
import sys
import os.path as op
import nipype.interfaces.ants as ants           # ANTS
import nipype.interfaces.freesurfer as fs
import nipype.interfaces.fsl as fsl
import pysbr.interfaces as sbrifs


def normalization( name='PySBR_Normalization', seg_params={}, reg_params={}, evaluate=False, conform_size=(2.0,2.0,2.0) ):
    """
    Creates the main PySBR registration workflow

    """

    if evaluate:
        name+= 'Ev'

    pipeline = pe.Workflow( name=name )

    def _aslist( tname ):
        import numpy as np
        return np.atleast_1d( tname ).tolist()

    inputnode = pe.Node( niu.IdentityInterface(
                fields=[ 'in_file', 'template' ] ),
                name='inputnode' )

    outputnode = pe.Node( niu.IdentityInterface(
                 fields=[ 'out_registered',
                          'out_transform',
                          'out_init_bbox',
                          'out_rois',
                          'out_init',
                          'out_init_rois',
                          'out_init_transform',
                          'out_init_composite' ] ),
                 name='outputnode' )
    tdir = op.join( sbrifs.get_datapath(), 'template' )
    tmpl_src = pe.Node( sbrifs.TemplateSource(base_dir=tdir), name='Template' )

    resample = pe.Node( fs.MRIConvert(vox_size=conform_size, resample_type='cubic', out_datatype='float'), name='Conform' )


    # ANTS Registration and transform application
    init = initialization()
    init_bbox = pe.Node( ants.ApplyTransforms(dimension=3,interpolation='NearestNeighbor' ), name='InitBBox' )

    # Transform template's hotspots and skullspots
    tfm_spots = pe.Node( sbrifs.TransformPoints(), name='InitLandmarks' )

    # Findspots
    # Inputs:  in_file, in_roi, rseedratio, rspotratio, maxvolume, uphill
    # Outputs: out_labels, out_spotfile
    fspots = pe.Node( sbrifs.FindSpots(), name='FindSubjectLandmarks' )

    for k in seg_params.keys():
        fspots.set_input(k, seg_params[k])
    
    sbreg = pe.Node( sbrifs.SBR(), name='ShapeMatching' )
    
    # Apply registration. 
    # Inputs:  in_target, in_source, in_target_points, in_source_points, in_fixed_points, stiffness, alpha
    # Outputs: out_moving
    ebsreg = pe.Node( sbrifs.EBSTransform(), name='EBSTransform' )
    
    for k in reg_params.keys():
        ebsreg.set_input(k, reg_params[k])

    pipeline.connect([
                         ( inputnode,     tmpl_src, [ ('template','name') ])
                        ,( inputnode,     resample, [ ('in_file', 'in_file' ) ])
                        ,( resample,          init, [(('out_file',_aslist),'moving_image') ])
                        ,( tmpl_src,          init, [(('out_file',_aslist),'fixed_image') ])
                        ,( resample,     init_bbox, [ ('out_file','reference_image') ])
                        ,( tmpl_src,     init_bbox, [ ('out_bbox', 'input_image') ])
                        ,( init,         init_bbox, [ ('forward_transforms','transforms'),('reverse_invert_flags','invert_transform_flags') ])
                        ,( init,         tfm_spots, [ ('forward_transforms','in_initialization'), ('forward_invert_flags', 'in_init_flags') ])
                        ,( inputnode,    tfm_spots, [ ('in_file','in_target') ])
                        ,( tmpl_src,     tfm_spots, [ ('out_spots','in_points'),('out_file','in_source') ])
                        ,( init_bbox,       fspots, [ ('output_image', 'in_roi' ) ])
                        ,( resample,        fspots, [ ('out_file', 'in_file' ) ])
                        ,( fspots,           sbreg, [ ('out_spotfile','trg_spots') ])
                        ,( tfm_spots,        sbreg, [ ('out_points','src_spots') ])
                        # Connect EBS registration
                        ,( resample,        ebsreg, [ ('out_file','in_target') ])
                        ,( tmpl_src,        ebsreg, [ ('out_file','in_source') ])
                        ,( init,            ebsreg, [ ('reverse_transforms', 'in_initialization' ), ('reverse_invert_flags','in_init_flags')   ])
                        ,( sbreg,           ebsreg, [ ('out_transform','in_transform') ])
                        # Connect outputs
                        ,( ebsreg,      outputnode, [ ('out_transform','out_transform'),('out_moving','out_registered'),('out_init','out_init') ])
                        ,( init_bbox,   outputnode, [ ('output_image','out_init_bbox' ) ] )
                        ,( init,        outputnode, [ ('forward_transforms','out_init_transform' ),('composite_transform', 'out_init_composite') ])
                    ])

    if evaluate:
        # ROIs mapping for evaluation
        init_rois = pe.MapNode( ants.ApplyTransforms(dimension=3), name='InitROIs', iterfield=['input_image'] )
        pipeline.connect([
                            # Initialize ROIs
                             ( resample,     init_rois, [ ('out_file','reference_image') ])
                            ,( tmpl_src,     init_rois, [ ('out_rois','input_image') ])
                            ,( init,          init_rois, [ ('forward_transforms','transforms'),('reverse_invert_flags','invert_transform_flags') ])
                            ,( init_rois,   outputnode, [ ('output_image','out_init_rois' ) ] )
                            # Connect EBS registration
                            ,( tmpl_src,        ebsreg, [ ('out_rois','in_apply') ])
                            # Connect registered ROIs to outputnode
                            ,( ebsreg,      outputnode, [ ('out_apply','out_rois') ])
                        ])

    return pipeline

def initialization( name='Initialization' ):
    """
    Creates a proxy node for the ANTs initialization

    """
    reg = pe.Node( ants.Registration() , name=name )
    reg.inputs.transforms = ['Rigid','Affine']
    reg.inputs.transform_parameters = [(1.0,),(2.0,), (1.0,)]
    reg.inputs.number_of_iterations = [[200],[100], [50] ]
    reg.inputs.dimension = 3
    reg.inputs.metric = ['Mattes']*2
    reg.inputs.metric_weight = [1.0]*2
    reg.inputs.radius_or_number_of_bins = [64]*3
    reg.inputs.sampling_strategy = ['Regular','Regular','Random']
    reg.inputs.sampling_percentage = [None,None,0.2]
    reg.inputs.convergence_threshold = [1.e-5,1.e-7,1.e-8]
    reg.inputs.convergence_window_size = [20,10,5]
    reg.inputs.smoothing_sigmas = [[16.0],[8.0],[4.0]]
    reg.inputs.sigma_units = ['mm']*3
    reg.inputs.shrink_factors = [[4],[2],[1]]
    reg.inputs.use_estimate_learning_rate_once = [True]*3
    reg.inputs.use_histogram_matching = [True]*3
    reg.inputs.winsorize_lower_quantile = 0.25
    reg.inputs.winsorize_upper_quantile = 0.95
    reg.inputs.initial_moving_transform_com = False
    reg.inputs.collapse_output_transforms = True
    return reg


def ants_normalization( name='ANTs_Normalization' ):
    """
    Creates the main ANTs non-linear, intensity-based registration workflow

    """

    pipeline = pe.Workflow( name=name )

    def _aslist( tname ):
        import numpy as np
        return np.atleast_1d( tname ).tolist()

    inputnode = pe.Node( niu.IdentityInterface(
                fields=[ 'in_file', 'template', 'initial_transforms' ] ),
                name='inputnode' )

    outputnode = pe.Node( niu.IdentityInterface(
                 fields=[ 'out_registered', 
                          'out_transforms',
                          'out_rois' ] ),
                 name='outputnode' )

    tdir = op.join( sbrifs.get_datapath(), 'template' )
    tmpl_src = pe.Node( sbrifs.TemplateSource(base_dir=tdir), name='Template' )   


    reg = pe.Node( ants.Registration() , name="NonlinearRefinement" )
    reg.inputs.transforms= [ 'SyN', 'SyN', 'SyN' ]
    reg.inputs.transform_parameters = [(1.0,4.0,6.0),(1.0,2.5,3.0),(0.5,0.8,1.0) ]
    reg.inputs.number_of_iterations = [ [20],[15],[15] ]
    reg.inputs.metric= [ ['Mattes'], ['CC'], ['CC'] ]
    reg.inputs.metric_weight= [ [1.0] , [1.0], [1.0] ]
    reg.inputs.radius_or_number_of_bins= [[64],[6],[1]]
    reg.inputs.sampling_strategy= [ ['Regular'], ['Regular'], ['Regular'] ]
    reg.inputs.sampling_percentage= [ [1.0], [1.0], [1.0] ]
    reg.inputs.convergence_threshold= [1.e-5, 1.e-7, 1.e-8]
    reg.inputs.convergence_window_size= [ 20, 10 ]
    reg.inputs.smoothing_sigmas= [ [8.0], [4.0], [1.0]  ]
    reg.inputs.sigma_units= ['mm','mm', 'mm' ]
    reg.inputs.shrink_factors= [[2],[2],[1]]
    reg.inputs.use_estimate_learning_rate_once= [True]*3
    reg.inputs.use_histogram_matching= [True]*3
    reg.inputs.winsorize_lower_quantile = 0.15

    applytfm = pe.Node( ants.ApplyTransforms(dimension=3), name='Apply' )
    tfm_bbox = pe.MapNode( ants.ApplyTransforms(dimension=3), iterfield=['input_image'], name='ApplyROIs' )

    pipeline.connect([
                         ( inputnode,     tmpl_src, [ ('template','name') ])
                        ,( inputnode,          reg, [(('in_file',_aslist),'moving_image'), ('initial_transforms','initial_moving_transform') ])
                        ,( tmpl_src,           reg, [(('out_file',_aslist),'fixed_image') ])
                        ,( reg,         outputnode, [ ('forward_transforms','out_transforms' ) ])
                        ,( inputnode,     applytfm, [ ('in_file','reference_image') ])
                        ,( tmpl_src,      applytfm, [ ('out_file', 'input_image') ])
                        ,( reg,           applytfm, [ ('reverse_transforms','transforms'),('reverse_invert_flags','invert_transform_flags') ])
                        ,( applytfm,    outputnode, [ ('output_image', 'out_registered' ) ])
                        ,( inputnode,     tfm_bbox, [ ('in_file','reference_image') ])
                        ,( tmpl_src,      tfm_bbox, [ ('out_rois', 'input_image') ])
                        ,( reg,           tfm_bbox, [ ('reverse_transforms','transforms'),('reverse_invert_flags','invert_transform_flags') ])
                        ,( tfm_bbox,    outputnode, [ ('output_image', 'out_rois' ) ])

                    ])
    return pipeline



# In case of command line run, execute one case
if __name__== '__main__':
    from argparse import ArgumentParser
    from argparse import RawTextHelpFormatter
    from os import getcwd
    
    parser = ArgumentParser(description='Run registration for one case', 
                            formatter_class=RawTextHelpFormatter)

    group = parser.add_argument_group('Input')
    
    group.add_argument("-i", "--input",  action="store",  required=True, 
                       help="input file to be analyzed")

    group.add_argument("-T", "--template",  action="store",  choices=['simulated', 'normal'], 
                       help="template to be used", default='simulated')

    group.add_argument("-w", "--workdir", action="store", default="/tmp", 
                       help="directory to store intermediate results")


    group_seg = parser.add_argument_group('Segmentation')

    group.add_argument("-g", "--gauss-prefilter-sigma",  action="store", default=0.0, 
                       type=float, dest='gauss_sigma', 
                       help="Prefilter the image with this gaussian.")
    
    group_seg.add_argument("-P", "--spot-mean-intensity-ratio",  action="store", 
                           type=float, default=1.6, dest='rspotmeanlimit', 
                           help="threshold to stop region growing with respect "
                           "to the ROI mean intensity ")

    group_seg.add_argument("-B", "--spot-bridge-intensity-ratio",  action="store", 
                           type=float, default=0.7, dest='rspotbridgelimit', 
                           help="threshold to stop region growing with respect to"
                           " the minimum intensity found on the line connecting the two hotspots")
    
    volume_group = group_seg.add_mutually_exclusive_group()
    
    volume_group.add_argument("-V", "--max-volume",  action="store", type=int, dest='maxvolume', 
                       help="Maximum volume (in pixel) to which a hot spot region can grow during during segmentation")
    
    volume_group.add_argument("-c", "--max-acceptable-volume",  action="store", type=int, 
                              dest='max_acceptable_volume', help="Maximum activation volume (in pixel) "
                              "to accept for a segmentation, if larger, segmentation should be considered to have failed")

    
    group_reg = parser.add_argument_group('Registration')
    
    group_reg.add_argument("-f", "--stiffness", action="store", type=float, default="0.001", 
                           help='Stiffness of the spline. A stiffness of 0.0 results '
                                    'in the standard interpolating spline. A non-zero stiffness '
                                    'allows the spline to approximate rather than interpolate the '
                                    'landmarks. Stiffness values are usually rather small, '
                                    'typically in the range of 0.001 to 0.1. The approximating '
                                    'spline formulation is based on the short paper by R. Sprengel, '
                                    'K. Rohr, H. Stiehl. "Thin-Plate Spline Approximation for Image '
                                    'Registration". In 18th International Conference of the IEEE '
                                    'Engineering in Medicine and Biology Society. 1996.')
    group_reg.add_argument("-a", "--alpha", action="store", type=float, default="10.0", 
                           help="Alpha is related to Poisson's Ratio (\\nu) as "
                           "\\alpha = 12 (1-\\nu)-1. \\nu_{gold} ~ 0.43; \\nu_{foam} ~ 0.30")
                                
    group_other = parser.add_argument_group('Other')

    group_other.add_argument( '-E', '--evaluate', action='store', type=bool, default=False,
                              dest='evaluate', help='apply registration results to ROIs' )

    group_other.add_argument( '-G', '--output-graphs', action='store', type=bool, default=False,
                               dest='out_graphs', help='write out nipype flowcharts' )

    options = parser.parse_args()

    # copy segmentation params 
    segmentation_params = {}
    if options.maxvolume is not None:
        segmentation_params['maxvolume'] = options.maxvolume
    if options.max_acceptable_volume is not None:
        segmentation_params['maxacceptablevol'] = options.max_acceptable_volume
    segmentation_params['gauss_sigma'] = options.gauss_sigma
    segmentation_params['rspotmeanlimit'] = options.rspotmeanlimit
    segmentation_params['rspotbridgelimit'] = options.rspotbridgelimit

    # copy registration params 
    registration_params = {}
    registration_params['stiffness'] = options.stiffness
    registration_params['alpha'] = options.alpha
    
    work_dir = options.workdir
    
    if not op.exists( work_dir ):
        os.makedirs( work_dir )
    
    wf = normalization( evaluate=options.evaluate, 
                        seg_params=segmentation_params, 
                        reg_params=registration_params )
    
    wf.base_dir = work_dir

    wf.inputs.inputnode.template = options.template
    wf.inputs.inputnode.in_file = op.abspath( options.input )

    if options.out_graphs:
        wf.write_graph( format='svg' )
        wf.write_graph( format='pdf' )

    wf.run()
