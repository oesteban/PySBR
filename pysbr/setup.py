#!/usr/bin/env python
# coding: utf-8
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

"""

setup.py - Installation script

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

from distutils.core import setup
from distutils.extension import Extension

import numpy;
import subprocess

depvars = {}
depvars['include_dirs'] = [numpy.get_include()]
depvars['extra_compile_args'] = ['-DPY_ARRAY_UNIQUE_SYMBOL=sbr_ARRAY_API', '-std=c++0x']

extension = Extension('csbr', ['tools/findspot.cc', 'tools/thinning_3d.cc', 'tools/borderdistance.cc', 'tools/sbr.cc'], **depvars)
                      
setup(name='csbr',
      version='0.1',
      description='Functions for Shape Based Registration in Python',
      author='Gert Wollny and Oscar Esteban',
      author_email=['gw.fossdev@gmail.com ', 'code@oscaresteban.es'],
      url='https://github.com/oesteban/PySBR', 
      requires=['numpy'], 
      ext_modules = [extension]
      )





