/** 
 * @file       borderdistance.cc
 * @brief      
 * @Author     Gert Wollny (gw.fossdev@gmail.com)
 * @date       September, 2013
 * @ingroup    PySBR
 *
 * 
 * 
 */

/* 
 * Copyright (c) 2013, code@oscaresteban.es (Oscar Esteban), and 
 *                     gw.fossdev@gmail.com (Gert Wollny)
 *                     with Biomedical Image Technology, UPM (BIT-UPM)
 * All rights reserved.
 * This file is part of PySBR.
 *
 */


#include "sbr.hh"
#include <numpy/noprefix.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <climits>

#include <queue>
#include <stdexcept>
#include <iostream>
#include <algorithm>
#include <stdio.h>

using std::queue; 
using std::priority_queue; 
using std::vector; 
using std::ostream; 
using std::transform; 
using std::cout; 
using std::invalid_argument; 

#ifdef __GNUC__
#define UNUSED(x) x __attribute__((unused))
#else 
#define UNUSED(x) x
#endif 


#if NPY_API_VERSION < 0x00000007
#define NPY_ARRAY_WRITEABLE NPY_WRITEABLE
#define NPY_ARRAY_ALIGNED NPY_ALIGNED
#endif 


static PyObject *do_get_distance_stats(PyArrayObject *image, PyObject *pyspot)
{
        int spot[3]; 
        if (!PyArg_ParseTuple(pyspot, "iii", &spot[0], &spot[1], &spot[2])) 
                throw invalid_argument("Unable to understand the spot tuple");
        
        if (image->descr->type_num != NPY_BOOL) 
                throw invalid_argument("input must be a binary mask");

        const TMirrorPythonArray<Bool> data(image); 
        npy_intp *dims = data.get_dims(); 
        

        double mean = 0.0; 
        double var = 0.0; 
        int n = 0; 
        double gravcenter[3] = {0,0,0}; 
        for (int z = 0; z < dims[2]; ++z) {
                double dz = z - spot[2]; 
		for (int y = 0; y < dims[1]; ++y) {
                        double dy = y - spot[1]; 
			for (int x = 0; x < dims[0]; ++x) {
                                double dx = x - spot[0];
                                if (data(x,y,z)) {
                                        double v = sqrt(dx * dx + dy * dy + dz * dz);
                                        mean += v; 
                                        var += v*v; 
                                        gravcenter[0] += x; 
                                        gravcenter[1] += y; 
                                        gravcenter[2] += z; 
                                        ++n; 
                                }
                        }
                }
        }
        // collected values

        double dgs = 0.0; 
        if (n > 0) {
                gravcenter[0] /= n; 
                gravcenter[1] /= n; 
                gravcenter[2] /= n;
                
                gravcenter[0] -= spot[0]; 
                gravcenter[1] -= spot[1]; 
                gravcenter[2] -= spot[2];
                
                // distance border-gravity center 
                dgs = sqrt(gravcenter[0]*gravcenter[0] + 
                           gravcenter[1]*gravcenter[1] + 
                           gravcenter[2]*gravcenter[2]); 
        }

        // mean border-spot distance 
        if (n > 0) 
                mean /= n; 
        if (n > 1) 
                var = sqrt((var - n * mean * mean) / (n - 1));

        return Py_BuildValue("ddd", dgs, mean, var); 

}


extern PyObject *Py_get_distance_stats(PyObject *UNUSED(self), PyObject *args)
{
	PyArrayObject *image = NULL; 
        PyObject *pyspot = NULL; 

        
	if (!PyArg_ParseTuple(args,  "O!O", &PyArray_Type, &image, &pyspot))
		return NULL;

        if (!PyTuple_Check(pyspot) || (PyTuple_Size(pyspot) != PyArray_NDIM(image))) {
                throw invalid_argument("The spot coordinate ust be a tuple type and it must be "
                                       "of the same size like the image dimenstion.");
        }
        
	return do_get_distance_stats(image, pyspot); 
}
