/**  
 * @file       thinning_3d.cc
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


/** The implementation is based on 
    
 */

#include <Python.h>
#if PY_MAJOR_VERSION >= 3
#define IS_PYTHON3
#endif

/* this is imported in the main module */
#define NO_IMPORT_ARRAY

#include <numpy/arrayobject.h>
#include <numpy/noprefix.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include <stdio.h>

#if NPY_API_VERSION < 0x00000007
#define NPY_ARRAY_WRITEABLE NPY_WRITEABLE
#define NPY_ARRAY_ALIGNED NPY_ALIGNED
#endif 

#ifdef __GNUC__
#define UNUSED(x) x __attribute__((unused))
#else 
#define UNUSED(x) x
#endif 


static const int eulerLUT[128] = {
	  1,  -1,  -1,   1, 
	 -3,  -1,  -1,   1, 
	 -1,   1,   1,  -1, 
	  3,   1,   1,  -1, 
	 -3,  -1,   3,   1, 
	  1,  -1,   3,   1, 
	 -1,   1,   1,  -1, 
	  3,   1,   1,  -1, 

	 -3,   3,  -1,   1, 
	  1,   3,  -1,   1, 
	 -1,   1,   1,  -1, 
	  3,   1,   1,  -1, 
	  1,   3,   3,   1, 
	  5,   3,   3,   1, 
	 -1,   1,   1,  -1, 
	  3,   1,   1,  -1, 

	 -7,  -1,  -1,   1, 
	 -3,  -1,  -1,   1, 
	 -1,   1,   1,  -1, 
	  3,   1,   1,  -1, 
	 -3,  -1,   3,   1, 
	  1,  -1,   3,   1, 
	 -1,   1,   1,  -1, 
	  3,   1,   1,  -1, 

	 -3,   3,  -1,   1, 
	  1,   3,  -1,   1, 
	 -1,   1,   1,  -1, 
	  3,   1,   1,  -1, 
	  1,   3,   3,   1, 
	  5,   3,   3,   1, 
	 -1,   1,   1,  -1, 
	  3,   1,   1,  -1
}; 
 

static void octree_labeling(int octant, int label, int *cube); 

static int not_Euler_invariant(int *neighbors); 

static int is_simple_point(int *neighborhood); 

static int get_neighbourhood(const Bool *base, int *neighbourhood, npy_intp *n, npy_intp *dims, npy_intp *strides) 
{

	// read the neigborhood of the current pixel 
	memset(neighbourhood, 0, 27 * sizeof(int));
	int nneighbours = -1; 
	
	for (int nn2 = -1; nn2 < 2; ++nn2) {
		int i2 = n[2] + nn2; 
		if (i2 < 0 || i2 >= dims[2]) 
			continue; 
		int ofs2 = i2 * strides[2]; 
		
		for (int nn1 = -1; nn1 < 2; ++nn1) {
			int i1 = n[1] + nn1; 
			if (i1 < 0 || i1 >= dims[1]) 
				continue; 
			int ofs1 = i1 * strides[1] + ofs2; 
			for (int nn0 = -1; nn0 < 2; ++nn0) {

				int i0 = n[0] + nn0; 
				if (i0 < 0 || i0 >= dims[0]) 
					continue; 
				Bool val = base[ofs1 + i0 * strides[0]]; 
				
				neighbourhood[nn0  + 3 * nn1 + 9 * nn2 + 13] = val; 
				if (val)
					++nneighbours;
				
				
			}
		}
	}
	
	return nneighbours; 
}

typedef struct _pixel_queue pixel_queue; 
struct _pixel_queue {
	pixel_queue *next; 
	npy_intp coord[3]; 
	Bool *ptr; 
}; 

static pixel_queue *queue_push(pixel_queue *q, npy_intp *c, Bool *p)
{
	pixel_queue *r = (pixel_queue *)malloc(sizeof(pixel_queue)); 
	r->next = q; 
	memcpy(r->coord, c, 3 * sizeof(npy_intp)); 
	r->ptr = p; 
	return r; 
}

static void queue_top(pixel_queue *q, npy_intp *c, Bool **p)
{
	assert(q); 
	assert(p); 
	assert(c);
	
	*p = q->ptr; 
	memcpy(c, q->coord,  3 * sizeof(npy_intp)); 
}

static pixel_queue *queue_pop(pixel_queue *q)
{
	assert(q); 
	pixel_queue *r = q->next; 
	free(q); 
	return r; 
}

static void queue_free(pixel_queue *q)
{
	while (q) {
		pixel_queue *r = q; 
		q = q->next; 
		free(r);
	}
}

static inline int queue_empty(pixel_queue *q)
{
	return q == NULL; 
}

static const char pixels[28] = "****B*****N*WXO*S*****U****";

static int thinning3d(PyArrayObject *data, int maxiter)
{
	Bool *base_ptr; 
	int ndims; 
	npy_intp *dims; 
	npy_intp *strides; 

	int neighbourhood[27];  
	int pixel_changed; 
	
	const int border_indices[6] = {10, 16, 22, 4, 14, 12 }; 

	
	assert(data); 
	ndims = PyArray_NDIM(data); 
	assert(data->flags & NPY_ARRAY_WRITEABLE); 
	assert(data->flags & NPY_ARRAY_ALIGNED); 
	assert(data->descr->type_num == NPY_BOOL); 
	
	if (ndims != 3) 
		return -1; 

	
	dims = PyArray_DIMS(data); 
	strides = PyArray_STRIDES(data); 
	base_ptr = (Bool*)PyArray_DATA(data);
	
	int iter = 0; 
	do {
		pixel_changed = 0; 
		for (int border = 0; border < 6; ++border) {
			int candidate_pixels = 0; 
			int local_changed_pixels = 0; 
			pixel_queue *q = NULL;
			int current_border_index = border_indices[border]; 
			npy_intp n[3]; 

			for (n[2] = 0; n[2] < dims[2]; ++n[2])
				for (n[1] = 0; n[1] < dims[1]; ++n[1])
					for (n[0] = 0; n[0] < dims[0]; ++n[0]) {
						int nneighbours; 
						Bool *ptr = base_ptr + (n[0] * strides[0] + n[1] * strides[1] + n[2] * strides[2]); 
									
						if (!*ptr)
							continue; 
						
						nneighbours = get_neighbourhood(base_ptr, neighbourhood, n, dims, strides); 
							
						// This is not a border pixel? 
						if (neighbourhood[current_border_index]) 
							continue;

						// current point has only one or no neighbour == not deletable 
						if (nneighbours < 2)
							continue;
						
						int not_euler = not_Euler_invariant(neighbourhood); 
						if (not_euler)
							continue; 

						// does deletion change connectivity? 
						if (!is_simple_point(neighbourhood)) 
							continue;  
						
						q = queue_push(q, n, ptr); 
						++candidate_pixels; 
					}
			while (!queue_empty(q)) {
				Bool *pixel; 
				queue_top(q, n, &pixel); 
				q = queue_pop(q);
				
				assert(*pixel != 0); 
				
                                get_neighbourhood(base_ptr, neighbourhood, n, dims, strides); 
                                
                                *pixel = 0; 
                                
                                
                                // check if pixel really can be deleted
                                if (!is_simple_point(neighbourhood)) {
                                        *pixel = 1; 
                                }else {
                                        ++local_changed_pixels; 
                                        ++pixel_changed;
                                }
			}
			queue_free(q);
		} // for borders

		++iter; 
	} while (pixel_changed > 0 && (maxiter > 0 && iter < maxiter) ); 
	

	return 0; 
}

PyObject *Py_thinning3d(PyObject *UNUSED(self), PyObject *args)
{
	PyArrayObject *image = NULL; 
        int maxiter = 0; 

	if (!PyArg_ParseTuple(args,"O!|i", &PyArray_Type, &image, &maxiter))
		return NULL;

	assert(image); 
		

	int result = thinning3d(image, maxiter); 
	return Py_BuildValue("i", result); 
}

static int is_simple_point(int *neighbourhood)
{
	int cube[26];
	int i; 
	// copy all save the center pixel
	for(i = 0; i < 13; i++ )
		cube[i] = neighbourhood[i];
	for(i = 14; i < 27; i++ ) 
		cube[i-1] = neighbourhood[i];

	// set initial label
	int label = 2;
	
	for(i = 0; i < 26; i++ ) {
		if( cube[i]!=1 )
			continue; 
		// voxel has not been labelled yet
		// start recursion with any octant that contains the point i
		switch( i ){
		case 0:
		case 1:
		case 3:
		case 4:
		case 9:
		case 10:
		case 12:
			octree_labeling(1, label, cube );
			break;
		case 2:
		case 5:
		case 11:
		case 13:
			octree_labeling(2, label, cube );
			break;
		case 6:
		case 7:
		case 14:
		case 15:
			octree_labeling(3, label, cube );
			break;
		case 8:
		case 16:
			octree_labeling(4, label, cube );
			break;
		case 17:
		case 18:
		case 20:
		case 21:
			octree_labeling(5, label, cube );
			break;
		case 19:
		case 22:
			octree_labeling(6, label, cube );
			break;
		case 23:
		case 24:
			octree_labeling(7, label, cube );
			break;
		case 25:
			octree_labeling(8, label, cube );
			break;
		}
		++label;
		
		if( label-2 >= 2 ) {
			return 0;
		}
	}
	return 1;
}

/** 
 * Check for Euler invariance. (see [Lee94])
 */
static int not_Euler_invariant(int *neighbors)
{
	// calculate Euler characteristic for each octant and sum up
	int EulerChar = 0;
	unsigned char n;
	// Octant SWU
	n = 0;
	if( neighbors[24] )
		n |=  64;
	if( neighbors[25] )
		n |=  32;
	if( neighbors[15] )
		n |=  16;
	if( neighbors[16] )
		n |=   8;
	if( neighbors[21] )
		n |=   4;
	if( neighbors[22] )
		n |=   2;
	if( neighbors[12] )
		n |=   1;
	EulerChar += eulerLUT[n];
// Octant SEU
	n = 0;
	if( neighbors[26] )
		n |=  64;
	if( neighbors[23] )
		n |=  32;
	if( neighbors[17] )
		n |=  16;
	if( neighbors[14] )
		n |=   8;
	if( neighbors[25] )
		n |=   4;
	if( neighbors[22] )
		n |=   2;
	if( neighbors[16] )
		n |=   1;
	EulerChar += eulerLUT[n];
// Octant NWU
	n = 0;
	if( neighbors[18] )
		n |=  64;
	if( neighbors[21] )
		n |=  32;
	if( neighbors[9] )
		n |=  16;
	if( neighbors[12] )
		n |=   8;
	if( neighbors[19] )
		n |=   4;
	if( neighbors[22] )
		n |=   2;
	if( neighbors[10] )
		n |=   1;
	EulerChar += eulerLUT[n];
// Octant NEU
	n = 0;
	if( neighbors[20] )
		n |=  64;
	if( neighbors[23] )
		n |=  32;
	if( neighbors[19] )
		n |=  16;
	if( neighbors[22] )
		n |=   8;
	if( neighbors[11] )
		n |=   4;
	if( neighbors[14] )
		n |=   2;
	if( neighbors[10] )
		n |=   1;
	EulerChar += eulerLUT[n];
// Octant SWB
	n = 0;
	if( neighbors[6] )
		n |=  64;
	if( neighbors[15] )
		n |=  32;
	if( neighbors[7] )
		n |=  16;
	if( neighbors[16] )
		n |=   8;
	if( neighbors[3] )
		n |=   4;
	if( neighbors[12] )
		n |=   2;
	if( neighbors[4] )
		n |=   1;
	EulerChar += eulerLUT[n];
// Octant SEB
	n = 0;
	if( neighbors[8] )
		n |=  64;
	if( neighbors[7] )
		n |=  32;
	if( neighbors[17] )
		n |=  16;
	if( neighbors[16] )
		n |=   8;
	if( neighbors[5] )
		n |=   4;
	if( neighbors[4] )
		n |=   2;
	if( neighbors[14] )
		n |=   1;
	EulerChar += eulerLUT[n];
// Octant NWB
	n = 0;
	if( neighbors[0] )
		n |=  64;
	if( neighbors[9] )
		n |=  32;
	if( neighbors[3] )
		n |=  16;
	if( neighbors[12] )
		n |=   8;
	if( neighbors[1] )
		n |=   4;
	if( neighbors[10] )
		n |=   2;
	if( neighbors[4] )
		n |=   1;
	EulerChar += eulerLUT[n];
// Octant NEB
	n = 0;
	if( neighbors[2] )
		n |=  64;
	if( neighbors[1] )
		n |=  32;
	if( neighbors[11] )
		n |=  16;
	if( neighbors[10] )
		n |=   8;
	if( neighbors[5] )
		n |=   4;
	if( neighbors[4] )
		n |=   2;
	if( neighbors[14] )
		n |=   1;
	EulerChar += eulerLUT[n];
	return EulerChar; 
}

static void octree_labeling(int octant, int label, int *cube)
{
	// check if there are points in the octant with value 1
	if( octant==1 ){
		// set points in this octant to current label
		// and recursive labeling of adjacent octants
		if( cube[0] == 1 )
			cube[0] = label;
		if( cube[1] == 1 ){
			cube[1] = label;        
			octree_labeling( 2, label, cube);
		}
		if( cube[3] == 1 )  {
			cube[3] = label;        
			octree_labeling( 3, label, cube);
		}
		if( cube[4] == 1 )  {
			cube[4] = label;        
			octree_labeling( 2, label, cube);
			octree_labeling( 3, label, cube);
			octree_labeling( 4, label, cube);
		}
		if( cube[9] == 1 )   {
			cube[9] = label;        
			octree_labeling( 5, label, cube);
		}
		if( cube[10] == 1 ) {
			cube[10] = label;        
			octree_labeling( 2, label, cube);
			octree_labeling( 5, label, cube);
			octree_labeling( 6, label, cube);
		}
		if( cube[12] == 1 )  {
			cube[12] = label;        
			octree_labeling( 3, label, cube);
			octree_labeling( 5, label, cube);
			octree_labeling( 7, label, cube);
		}
	}
	
	if( octant==2 ) {
		if( cube[1] == 1 ) {
			cube[1] = label;
			octree_labeling( 1, label, cube);
		}
		if( cube[4] == 1 ){
			cube[4] = label;        
			octree_labeling( 1, label, cube);
			octree_labeling( 3, label, cube);
			octree_labeling( 4, label, cube);
		}
		if( cube[10] == 1 )  {
			cube[10] = label;        
			octree_labeling( 1, label, cube);
			octree_labeling( 5, label, cube);
			octree_labeling( 6, label, cube);
		}
		if( cube[2] == 1 )
			cube[2] = label;        
		if( cube[5] == 1 )  {
			cube[5] = label;        
			octree_labeling( 4, label, cube);
		}
		if( cube[11] == 1 ) {
			cube[11] = label;        
			octree_labeling( 6, label, cube);
		}
		if( cube[13] == 1 )  {
			cube[13] = label;        
			octree_labeling( 4, label, cube);
			octree_labeling( 6, label, cube);
			octree_labeling( 8, label, cube);
		}
	}

	if( octant==3 ) {
		
		if( cube[3] == 1 ) {
			cube[3] = label;        
			octree_labeling( 1, label, cube);
		}
		if( cube[4] == 1 ) {
			cube[4] = label;        
			octree_labeling( 1, label, cube);
			octree_labeling( 2, label, cube);
			octree_labeling( 4, label, cube);
		}
		if( cube[12] == 1 ) {
			cube[12] = label;        
			octree_labeling( 1, label, cube);
			octree_labeling( 5, label, cube);
			octree_labeling( 7, label, cube);
		}
		if( cube[6] == 1 )
			cube[6] = label;        
		if( cube[7] == 1 )  {
			cube[7] = label;        
			octree_labeling( 4, label, cube);
		}
		if( cube[14] == 1 ){
			cube[14] = label;        
			octree_labeling( 7, label, cube);
		}
		if( cube[15] == 1 )
		{
			cube[15] = label;        
			octree_labeling( 4, label, cube);
			octree_labeling( 7, label, cube);
			octree_labeling( 8, label, cube);
		}
	}
	if( octant==4 )
	{
		if( cube[4] == 1 )
		{
			cube[4] = label;        
			octree_labeling( 1, label, cube);
			octree_labeling( 2, label, cube);
			octree_labeling( 3, label, cube);
		}
		if( cube[5] == 1 )
		{
			cube[5] = label;        
			octree_labeling( 2, label, cube);
		}
		if( cube[13] == 1 )
		{
			cube[13] = label;        
			octree_labeling( 2, label, cube);
			octree_labeling( 6, label, cube);
			octree_labeling( 8, label, cube);
		}
		if( cube[7] == 1 )
		{
			cube[7] = label;        
			octree_labeling( 3, label, cube);
		}
		if( cube[15] == 1 )
		{
			cube[15] = label;        
			octree_labeling( 3, label, cube);
			octree_labeling( 7, label, cube);
			octree_labeling( 8, label, cube);
		}
		if( cube[8] == 1 )
			cube[8] = label;        
		if( cube[16] == 1 )
		{
			cube[16] = label;        
			octree_labeling( 8, label, cube);
		}
	}
	if( octant==5 )
	{
		if( cube[9] == 1 )
		{
			cube[9] = label;        
			octree_labeling( 1, label, cube);
		}
		if( cube[10] == 1 )
		{
			cube[10] = label;        
			octree_labeling( 1, label, cube);
			octree_labeling( 2, label, cube);
			octree_labeling( 6, label, cube);
		}
		if( cube[12] == 1 )
		{
			cube[12] = label;        
			octree_labeling( 1, label, cube);
			octree_labeling( 3, label, cube);
			octree_labeling( 7, label, cube);
		}
		if( cube[17] == 1 )
			cube[17] = label;        
		if( cube[18] == 1 )
		{
			cube[18] = label;        
			octree_labeling( 6, label, cube);
		}
		if( cube[20] == 1 )
		{
			cube[20] = label;        
			octree_labeling( 7, label, cube);
		}
		if( cube[21] == 1 )
		{
			cube[21] = label;        
			octree_labeling( 6, label, cube);
			octree_labeling( 7, label, cube);
			octree_labeling( 8, label, cube);
		}
	}
	if( octant==6 )
	{
		if( cube[10] == 1 )
		{
			cube[10] = label;        
			octree_labeling( 1, label, cube);
			octree_labeling( 2, label, cube);
			octree_labeling( 5, label, cube);
		}
		if( cube[11] == 1 )
		{
			cube[11] = label;        
			octree_labeling( 2, label, cube);
		}
		if( cube[13] == 1 )
		{
			cube[13] = label;        
			octree_labeling( 2, label, cube);
			octree_labeling( 4, label, cube);
			octree_labeling( 8, label, cube);
		}
		if( cube[18] == 1 )
		{
			cube[18] = label;        
			octree_labeling( 5, label, cube);
		}
		if( cube[21] == 1 )
		{
			cube[21] = label;        
			octree_labeling( 5, label, cube);
			octree_labeling( 7, label, cube);
			octree_labeling( 8, label, cube);
		}
		if( cube[19] == 1 )
			cube[19] = label;        
		if( cube[22] == 1 )
		{
			cube[22] = label;        
			octree_labeling( 8, label, cube);
		}
	}
	if( octant==7 )
	{
		if( cube[12] == 1 )
		{
			cube[12] = label;        
			octree_labeling( 1, label, cube);
			octree_labeling( 3, label, cube);
			octree_labeling( 5, label, cube);
		}
		if( cube[14] == 1 )
		{
			cube[14] = label;        
			octree_labeling( 3, label, cube);
		}
		if( cube[15] == 1 )
		{
			cube[15] = label;        
			octree_labeling( 3, label, cube);
			octree_labeling( 4, label, cube);
			octree_labeling( 8, label, cube);
		}
		if( cube[20] == 1 )
		{
			cube[20] = label;        
			octree_labeling( 5, label, cube);
		}
		if( cube[21] == 1 )
		{
			cube[21] = label;        
			octree_labeling( 5, label, cube);
			octree_labeling( 6, label, cube);
			octree_labeling( 8, label, cube);
		}
		if( cube[23] == 1 )
			cube[23] = label;        
		if( cube[24] == 1 )
		{
			cube[24] = label;        
			octree_labeling( 8, label, cube);
		}
	}
	if( octant==8 )
	{
		if( cube[13] == 1 )
		{
			cube[13] = label;        
			octree_labeling( 2, label, cube);
			octree_labeling( 4, label, cube);
			octree_labeling( 6, label, cube);
		}
		if( cube[15] == 1 )
		{
			cube[15] = label;        
			octree_labeling( 3, label, cube);
			octree_labeling( 4, label, cube);
			octree_labeling( 7, label, cube);
		}
		if( cube[16] == 1 )
		{
			cube[16] = label;        
			octree_labeling( 4, label, cube);
		}
		if( cube[21] == 1 )
		{
			cube[21] = label;        
			octree_labeling( 5, label, cube);
			octree_labeling( 6, label, cube);
			octree_labeling( 7, label, cube);
		}
		if( cube[22] == 1 )
		{
			cube[22] = label;        
			octree_labeling( 6, label, cube);
		}
		if( cube[24] == 1 )
		{
			cube[24] = label;        
			octree_labeling( 7, label, cube);
		}
		if( cube[25] == 1 )
			cube[25] = label;        
	} 

}


