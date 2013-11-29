#!/usr/bin/env python
# coding: utf-8
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

"""
polar.py: This file gathers the necessary functions and methods
          for polar coordinates methods in PySBR


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


import numpy as np
import nibabel as nb

class Registration:
    def set_reference( self, im ):
        self.im1 = (im-im.mean())/im.std()
    
    def set_moving( self, im ):
        self.im2 = (im-im.mean())/im.std()
    
    def xcorr( self, p=[0,0] ):
        from scipy.signal import correlate2d
        from scipy.ndimage import affine_transform

        im1 = self.im1.copy()
        im2 = self.im2.copy()
        im1_tf = affine_transform(im1, np.identity(2) ,mode='wrap',offset=p )
        result = correlate2d( (im1_tf - im1_tf.mean()), (im2 - im2.mean()), mode='full', boundary='wrap' )
        return np.sum( result ) / result.size
        
    def msd( self, p=[0,0] ):
        from scipy.ndimage import affine_transform

        im1 = self.im1.copy()
        im2 = self.im2.copy()
        im1_tf = affine_transform(im1, np.identity(2) ,mode='wrap',offset=p )
        r = (im1_tf.reshape(-1) - im2.reshape(-1))**2
        return r.sum()/r.size
        
    def optimize( self, p0=[0,0], iterprint=False ):
        from scipy.optimize import minimize

        printp = None
        if iterprint: 
            def printp( xk ):
                print xk
        
        self.res = minimize( self.metric, p0, method='Powell', tol=1e-5, options={'disp': False, 'maxiter': 200}, callback=printp)
        return self.res
       
    def __init__(self):
        self.im1 = None
        self.im2 = None
        self.res = None
        self.mfunc = 'msd'
        
    def __init__( self, ref, mov, mfunc='msd' ):
        self.res = None
        self.set_reference( ref )
        self.set_moving( mov )
        self.mfunc = mfunc
        
    def metric( self, p ):
        if self.mfunc == 'xcorr':
            return self.xcorr( p )
        else:
            return self.msd( p )

    def get_result( self ):
        return self.res

def polar_rep( points, center, feat=None, res=[ 200, 400 ], interp='nearest' ):
    import numpy as np
    import matplotlib.mlab as mlab
    import scipy.ndimage as si
    from scipy.interpolate import griddata
    from collections import Counter
    from collections import defaultdict
    from itertools import chain
    
    points = points - center
    
    # Map coordinates
    R = []
    theta = []
    phi = []   
    for p in points:
        r = np.linalg.norm( p )
        R.append(r)
        theta.append( np.arccos( p[2]/r ) )
        phi.append( np.arctan2( p[1],p[0] ) )

    R = np.array( R )
    theta = np.around( np.array( theta ), decimals=3 )
    phi = np.around( np.array( phi ), decimals=3 )

    # Switch feature in case it is present
    if not feat is None:
        assert( R.size == np.shape( feat )[0] )
        R = [ tuple( val ) for val in feat ]
        # Remove duplicates (2nd, and following)
        pos = np.array( zip( theta, phi ) )
        dt = np.dtype( [ ('a', pos.dtype ), ('b', pos.dtype ) ] )
        posy = pos.view( dtype=dt ).squeeze()
        u,idx,inv= np.unique( posy, return_index=True, return_inverse=True )
        theta = theta[ idx ]
        phi = phi[ idx ]
        R = np.take( R, idx, axis=0 )
    else:       
        # Remove (closest) duplicates
        D = defaultdict( list )
        for i, item in enumerate( zip(theta,phi) ):
            D[item].append(i)
    
        dlist = []
        delcount = 0
        for k,v in D.items():
            if len(v)>1:
                Rs = [ R[i] for i in v ]
                maxid = np.argmax( Rs, axis=0 )
                dlist.append( [ i for i in v if i!=maxid ] )
                delcount+=1
   
        if delcount>0:
            # Flatten list of lists and delete elements
            dlist = np.fromiter( chain.from_iterable(dlist), dtype='int' )
            R = np.delete( R, dlist )
            theta = np.delete( theta, dlist )
            phi = np.delete( phi, dlist )


    # Expand data to be periodical before fitting
    if R.ndim > 1:
        expR = np.tile( R, (9,1) )
    else:
        expR = np.tile( R, 9 )

    expTheta = []
    expPhi = []
    a = np.linspace( -1.0, 1.0, 3 )
    incs = np.transpose([np.tile(a, len(a)), np.repeat(a, len(a))]) * np.array( ( np.pi, 2*np.pi ) )

    for q in incs:
        expTheta.append( (np.array( theta.copy())+ q[0]).tolist() )
        expPhi.append( (np.array( phi.copy())+ q[0]).tolist() )

    expTheta = np.reshape( expTheta, -1 ).tolist()
    expPhi = np.reshape( expPhi, -1).tolist()

    assert( np.shape( expR )[0] == len(expTheta) )
    assert( np.shape( expR )[0] == len(expPhi) )

    xs = np.linspace( 0.0, np.pi - (np.pi/res[0]), res[0] )
    ys = np.linspace( -np.pi, np.pi - (2*np.pi/res[1]), res[1] )
   
    X,Y = np.meshgrid( xs, ys )
    grid = np.transpose([np.tile(xs, len(ys)), np.repeat(ys, len(xs))])

    pointcloud = np.array( zip( expTheta, expPhi ) )
    Z = griddata( pointcloud, expR, grid , method=interp )
    Z = np.squeeze( Z.reshape( (X.shape[0],X.shape[1],-1) ) )
    Z = np.ma.masked_invalid( Z )
    return (X,Y,Z)
