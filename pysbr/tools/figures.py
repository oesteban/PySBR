#!/usr/bin/env python
# coding: utf-8
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

"""
figures.py - Generates the boxplots for the Frontiers paper

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



def plot_boxplot( in_file, out_file=None, grade_names=None, method_list=[ 'PySBR','ANTS' ], method_names=['PySBR','IBR'], roi_names=[ 'All', 'LCAU', 'LPUT', 'RCAU', 'RPUT' ], nograde=False, lims=[ -60.0, 40.0 ], legend=True, blocksize=(4,6) ):
    import numpy as np
    import os, os.path as op
    import pandas as pd
    import matplotlib
    from matplotlib import pyplot as plt
    from matplotlib.patches import Polygon
    from matplotlib.lines import Line2D

    if out_file is None:
        fname,fext = op.splitext( op.basename( in_file ) )
        out_file = op.abspath( fname + '.pdf' )

    matplotlib.rcParams.update({'font.size': 22})

    if nograde:
        dt = { 'method': np.str, 'subject': np.str, 'totalOverlap': np.float32, 'ov1': np.float32,'ov2': np.float32,'ov3': np.float32,'ov4': np.float32 }
        datacols = range(2,7)
    else:
        dt = { 'method': np.str, 'subject': np.str, 'grade': np.uint8,'totalOverlap': np.float32, 'ov1': np.float32,'ov2': np.float32,'ov3': np.float32,'ov4': np.float32 }
        datacols = range(3,8)

    df = pd.read_csv( in_file, header=None, dtype=dt)

    nSubplots = 1
    nMethods = len(method_list)

    if not grade_names is None:
        nSubplots = len( grade_names )

    boxColors = [ 'darkblue', 'deepskyblue', 'dodgerblue', 'limegreen', 'green' ]

    f, axs = plt.subplots( 1,nSubplots, sharex=True, sharey=True)
    axs = np.atleast_1d( axs )
   
    if legend: 
        f.set_size_inches( ( blocksize[0]*nSubplots + 4, blocksize[1]+2 ) )
        plt.subplots_adjust(wspace=0, left=0.07, right=0.99, top=0.90, bottom=0.25)
    else:
        f.set_size_inches( ( blocksize[0]*nSubplots + 4, blocksize[1] ) )
        plt.subplots_adjust(wspace=0, left=0.16, right=0.99, top=0.90, bottom=0.07 )
        


    plt.setp( [a.get_yticklabels() for a in f.axes[1:]], visible=False )
    plt.hold(True)
  
    overlaps = np.array( [ np.squeeze( np.array( df[col], dtype=np.float32 ) * 100) for col in datacols ] )
    nCols = np.shape(overlaps)[0]
    poss = np.array( range(1, nCols+1 ) )
    boxPolygons = []
    xticks_pos = [ 3, 10 ]

    axs[0].set_ylabel('Overlap improvement (%)')
    axs[0].set_ylim( lims[0], lims[1])

    for grade,ax in enumerate( np.atleast_1d( axs ) ):
        if not grade_names is None:
            ax.set_xlabel(grade_names[grade])
            ax.xaxis.set_label_position('top')

        ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
        ax.set_axisbelow(True)
        
        positions = [ j + i*(nCols+2) for i in range(len(method_list)) for j in range(1,nCols+1) ]  
        alldata = []
        
        for mid,method in enumerate(method_list):
            if nograde:
                idxs = np.squeeze( np.argwhere( (df[0]==method) ) )
            else:
                idxs = np.squeeze( np.argwhere( (df[0]==method) & (df[2]==grade) ) )

            d = overlaps[:,idxs].copy()
            
            initvals = []
            for i in idxs:
                subj_id = df[1][i]
                if nograde:
                    init_id = np.squeeze( np.argwhere( (df[0]=='Init') & (df[1]==subj_id) ) )
                else:
                    init_id = np.squeeze( np.argwhere( (df[0]=='Init') & (df[2]==grade) & (df[1]==subj_id) ) )

                initvals.append( overlaps[:,init_id] )
            
           
            initvals = np.array( initvals, dtype=np.float32 ).T
            axData = np.squeeze( np.array( d - initvals ) ).T
            medians = np.median( axData, axis=0 )
            alldata.append( axData.copy() )
        
        
        alldata = np.hstack( tuple( alldata ) )   
        bp = ax.boxplot(alldata, notch=0, sym='+', vert=True, whis=1.5, widths=0.9, positions=positions )
        
        plt.xlim( positions[0]-1, positions[-1]+1 )
        plt.xticks( xticks_pos, method_names  )   
        plt.setp(bp['boxes'], color='black')
        plt.setp(bp['medians'], color='black')
        plt.setp(bp['whiskers'], color='black')
        plt.setp(bp['fliers'], color='red', marker='+',markersize=15)
        
            
        upperLabels = []
        for i,box in enumerate( bp['boxes'] ):
            boxX = box.get_xdata().copy()
            boxY = box.get_ydata().copy()
            boxCoords = zip( boxX, boxY )
            
            boxPolygon = Polygon( boxCoords, facecolor=boxColors[i%nCols], label=roi_names[i%nCols] )
            boxPolygons.append(boxPolygon)
            ax.add_patch(boxPolygon)
            
            # Draw median lines
            med = bp['medians'][i]
            medX = med.get_xdata().copy()
            medY = med.get_ydata().copy()
            
            medians[i%nCols] = medY[0]
            
            ax.plot( [ np.average( med.get_xdata() ) ], [ np.average( alldata[:,i] ) ], color='w', marker='*', markersize=15, markeredgecolor='k')
        
            upperLabels+= [str(np.round(s, 2)) for s in medians] 
        
                            
            #pos = np.arange(len(method_list))+1
            #upperLabels = [str(np.round(s, 2)) for s in medians]
        
            #for tick,label in enumerate( ax.get_xticklabels() ):
            #    k = tick % len(method_list)
            #    ax.text(pos[tick], lims[1]-(lims[1]*0.15), upperLabels[tick],
            #            horizontalalignment='center', size='small', weight='semibold',
            #            color=boxColors[k])
        
    star = Line2D( [ 0.0 ], [ 0.0 ], color='w', marker='*', markersize=20, markeredgecolor='k')
    boxes = boxPolygons[0:len(roi_names)]
    boxes.append( star )
      
    names = [ a for a in roi_names ]
    names.append( 'ROI average' )

    if legend:        
        f.legend( boxes, names , ncol=(len(boxes)), loc=8, numpoints=1, bbox_to_anchor=( 0.5, 0.01 ) )

    plt.savefig( out_file, dpi=1200 )
    return out_file

