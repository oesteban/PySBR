#!/usr/bin/env python
# coding: utf-8
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

"""
evaluation.py - The PySBR evaluation workflow

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
import nipype.interfaces.io as nio
import os
import os.path as op
import glob
import nipype.interfaces.ants as ants           # ANTS
import nipype.interfaces.freesurfer as fs
import nipype.algorithms.misc as misc

import pysbr.interfaces as sbrifs
from pysbr.workflows.registration import normalization,ants_normalization_2


def evaluation_workflow( name='DaTSCAN_Evaluation' ):
    """
    Creates the main DaTSCAN evaluation workflow

    """

    pipeline = pe.Workflow( name=name )

    def _append_csv( in_row11, in_row12, in_row21, in_row22, in_row31, in_row32, subject_id, out_file=None ):
        import csv
        import os.path as op
        import numpy as np

        if out_file is None:
            out_file = op.abspath( 'stats.csv' )

        row1 = [ 'PySBR', subject_id ] + np.atleast_1d( in_row11 ).tolist() + np.atleast_1d( in_row12 ).tolist()
        row2 = [ 'ANTS' , subject_id ] + np.atleast_1d( in_row21 ).tolist() + np.atleast_1d( in_row22 ).tolist()
        row3 = [ 'Init' , subject_id ] + np.atleast_1d( in_row31 ).tolist() + np.atleast_1d( in_row32 ).tolist()

        with open( out_file, 'a' ) as f:
            w = csv.writer( f )
            w.writerow( row1 )
            w.writerow( row2 )
            w.writerow( row3 )

        return out_file

    def _build_paths( data_dir, subject_id ):
        import os.path as op
        in_file = op.join( data_dir, subject_id, 'spect.nii.gz' )
        in_rois = [ op.join( data_dir, subject_id, 'roi_%s.nii.gz' % t ) for t in [ 'c_l', 'p_l', 'c_r', 'p_r' ] ]
        return (in_file,in_rois)

    inputnode = pe.Node( niu.IdentityInterface(
                fields=[ 'data_dir', 'subject_id', 'template' ] ),
                name='inputnode' )
    ds = pe.Node( niu.Function( input_names=['data_dir','subject_id'], output_names=['in_file','in_rois'], function=_build_paths),
                  name='GrabSubject' )

    outputnode = pe.Node( niu.IdentityInterface(
                 fields=[ 'out_pysbr', 
                          'out_ants',
                          'out_pysbr_rois',
                          'out_ants_rois',
                          'out_init_rois',
                          'out_stats' ] ),
                 name='outputnode' )

    # Create ANTS non-linear registration workflow
    ants_reg = ants_normalization_2( name='IBR2')

    # Create PySBR non-linear registration workflow
    pysbr_reg = normalization( evaluate=True, name='PySBR' )

    select = pe.Node( niu.Select( index=[0] ), name='GetFirst')
    selectROI = pe.Node( niu.Select( index=[0] ), name='GetFirstROI')

    resample1 = pe.MapNode( fs.MRIConvert(out_datatype='float',resample_type='cubic' ), iterfield='in_file', name='Resample0' )
    resample2 = pe.MapNode( fs.MRIConvert(out_datatype='float',resample_type='cubic' ), iterfield='in_file', name='Resample1' )
    resample3 = pe.MapNode( fs.MRIConvert(out_datatype='float',resample_type='cubic' ), iterfield='in_file', name='Resample2' )

    # Connect to overlap measurement
    overlap1 = pe.Node( misc.FuzzyOverlap(), name='OverlapPySBR' )
    overlap2 = pe.Node( misc.FuzzyOverlap(), name='OverlapANTS' )
    overlap3 = pe.Node( misc.FuzzyOverlap(), name='OverlapInit' )

    # Connect to csvwriter
    tocsv = pe.Node( niu.Function( input_names=['in_row11','in_row21','in_row12','in_row22', 'in_row31', 'in_row32', 'subject_id'], output_names=['out_file'], function=_append_csv ), name='OutputStats' )

    pipeline.connect([
                         ( inputnode,    pysbr_reg, [ ('template','inputnode.template')])
                        ,( inputnode,     ants_reg, [ ('template','inputnode.template')])
                        ,( inputnode,           ds, [ ('data_dir','data_dir'),('subject_id','subject_id')])
                        ,( ds,           pysbr_reg, [ ('in_file','inputnode.in_file' ) ])
                        ,( ds,            ants_reg, [ ('in_file','inputnode.in_file' ) ])
                        ,( pysbr_reg,       select, [ ('outputnode.out_init_transform', 'inlist' ) ])
                        ,( select,        ants_reg, [ ('out', 'inputnode.initial_transforms' ) ] )
                        ,( pysbr_reg,    resample1, [ ('outputnode.out_rois','in_file' ) ])
                        ,( ants_reg,     resample2, [ ('outputnode.out_rois','in_file' ) ])
                        ,( pysbr_reg,    resample3, [ ('outputnode.out_init_rois','in_file' ) ])
                        ,( ds,           selectROI, [ ('in_rois','inlist' ) ])
                        ,( selectROI,    resample1, [ ('out','reslice_like' ) ])
                        ,( selectROI,    resample2, [ ('out','reslice_like' ) ])
                        ,( selectROI,    resample3, [ ('out','reslice_like' ) ])
                        ,( resample1,     overlap1, [ ('out_file','in_ref') ])
                        ,( resample2,     overlap2, [ ('out_file','in_ref') ])
                        ,( resample3,     overlap3, [ ('out_file','in_ref') ])
                        ,( ds,            overlap1, [ ('in_rois','in_tst' ) ])
                        ,( ds,            overlap2, [ ('in_rois','in_tst' ) ])
                        ,( ds,            overlap3, [ ('in_rois','in_tst' ) ])
                        ,( pysbr_reg,   outputnode, [ ('outputnode.out_registered','out_pysbr'),('outputnode.out_rois','out_pysbr_rois'),
                                                      ('outputnode.out_init_rois', 'out_init_rois' ) ])
                        ,( ants_reg,    outputnode, [ ('outputnode.out_registered','out_ants'), ('outputnode.out_rois','out_ants_rois') ])
                        ,( overlap1,         tocsv, [ ('dice','in_row11'),('class_fdi','in_row12') ])
                        ,( overlap2,         tocsv, [ ('dice','in_row21'),('class_fdi','in_row22') ])
                        ,( overlap3,         tocsv, [ ('dice','in_row31'),('class_fdi','in_row32') ])
                        ,( inputnode,        tocsv, [ ('subject_id','subject_id') ])
                        ,( tocsv,       outputnode, [ ('out_file','out_stats') ])
                    ])
    return pipeline


if __name__== '__main__':
    from argparse import ArgumentParser
    from argparse import RawTextHelpFormatter
    from os import getcwd
    from shutil import copyfileobj
    from pysbr.tools import figures as fig

    parser = ArgumentParser(description='Run registration for one case', 
                            formatter_class=RawTextHelpFormatter)
    
    g_input = parser.add_argument_group('Input')
    
    g_input.add_argument( '-S', '--subjects_dir', action='store',
                          default=os.getenv( 'PYSBR_SUBJECTS_DIR', '/media/mnemea/MINDt-Quantidopa/PySBR-SimulatedData/' ),
                          help='directory where subjects should be found' )

    g_input.add_argument( '-s', '--subject', action='store',
                          default='S*', help='subject id or pattern' )

    g_input.add_argument( "-T", "--template",  action="store",  choices=['simulated', 'normal'], 
                          help="template to be used", default='simulated')

    g_input.add_argument( '-w', '--work_dir', action='store', default=os.getcwd(),
                          help='directory to store intermediate results' )

    g_input.add_argument( '-N', '--name', action='store', default='PythonInNeuroscience2-Simulated',
                          help='default workflow name, it will create a new folder' )

    
    g_output = parser.add_argument_group( 'Output' )

    g_output.add_argument( '-o', '--out_csv', action='store', default='results.csv',
                           help='output csv file name collecting all evaluation results' )

    g_output.add_argument( '-O', '--out_pdf', action='store', default='results.pdf',
                           help='output pdf file with evaluation results plot' )

    options = parser.parse_args()

    if not op.exists( options.work_dir ):
        os.makedirs( options.work_dir )
  
    subjects = [ op.basename( sub ) for sub in glob.glob( op.join( options.subjects_dir, options.subject ) ) ] 

    # Limit execution (only debugging)
    #subjects = [ subjects[0] ]

    wf = pe.Workflow( name=options.name )
    wf.base_dir = options.work_dir


    infosource = pe.Node( niu.IdentityInterface(fields=['subject_id']),
                          name="infosource")

    infosource.iterables = [('subject_id', subjects)]
  
    ev = evaluation_workflow()
    ev.inputs.inputnode.template = 'simulated'
    ev.inputs.inputnode.data_dir = options.subjects_dir
    wf.connect([
                 ( infosource, ev, [ ('subject_id','inputnode.subject_id') ])
               ])

    try:
        wf.write_graph( graph2use='hierarchical', format='pdf', simple_form=True )
    except RuntimeError as e:
        print e

    try:
        wf.run()
    except RuntimeError as e:
        print e

    csv_paths = []
    out_csv = op.join( options.work_dir, options.name, options.out_csv )

    open( out_csv, 'w').close()

    for root, dirs, files, in os.walk( op.join( options.work_dir, options.name ) ):
        for name in files:
            fname, ext = op.splitext( name )
            if ext == '.csv':
                csv_paths+= [ op.join( root, name ) ]

    if len(csv_paths)>0:
        with open( out_csv, 'wb' ) as of:
            for fname in csv_paths:
                with open( fname ) as f:
                    copyfileobj( f, of )


        grade_list = [ 'Simulated' ]

            
        out_pdf = options.out_pdf
        if op.basename( out_pdf ) == out_pdf:
            out_pdf = op.join( options.work_dir, options.name, out_pdf )

        out = fig.plot_boxplot( out_csv, out_file=out_pdf, grade_names=grade_list, nograde=True, legend=False )
        print 'PDF file written: ' + out


