/**  
 * @file       sbr.cc
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


#include <Python.h>
#if PY_MAJOR_VERSION >= 3
#define IS_PYTHON3
#endif

#include <numpy/arrayobject.h>
#include <sstream>
#include <stdexcept>

#ifdef __GNUC__
#define UNUSED(x) x __attribute__((unused))
#else 
#define UNUSED(x) x
#endif 

static PyObject *SbrError; 

extern PyObject *Py_thinning3d(PyObject *UNUSED(self), PyObject *args); 
extern PyObject *Py_findspots3d(PyObject *UNUSED(self), PyObject *args, PyObject *kwdict);
extern PyObject *Py_get_distance_stats(PyObject *UNUSED(self), PyObject *args); 


PyObject *Py_findspots3d_call(PyObject *self, PyObject *args, PyObject *kwdict)
{
	std::ostringstream msg;			\
	try {
		return Py_findspots3d(self, args, kwdict); 
	} catch (std::runtime_error& x) {				
                msg << "csbr runtime error:'" << x.what() << "'";
        }                                                              
        catch (std::invalid_argument& x) {
                msg << "csbr invalid argument:'" << x.what() << "'";   
        }                                                              
        catch (std::exception& x) {
                msg << "csbr exception: '" << x.what() << "'";
        } 
        catch (...) { 
                msg << "csbr: unknown error";
        }                                   
        PyErr_SetString(SbrError, msg.str().c_str());
        return NULL;
}

PyObject *Py_get_distance_stats_call(PyObject *self, PyObject *args)
{
	std::ostringstream msg;			\
	try {
		return Py_get_distance_stats(self, args); 
	} catch (std::runtime_error& x) {				
                msg << "csbr runtime error:'" << x.what() << "'";
        }                                                              
        catch (std::invalid_argument& x) {
                msg << "csbr invalid argument:'" << x.what() << "'";   
        }                                                              
        catch (std::exception& x) {
                msg << "csbr exception: '" << x.what() << "'";
        } 
        catch (...) { 
                msg << "csbr: unknown error";
        }                                   
        PyErr_SetString(SbrError, msg.str().c_str());
        return NULL;
}


static struct PyMethodDef csbr_methods[]={
	{ "thinning3d", (PyCFunction)Py_thinning3d, METH_VARARGS, "run morphological thinning on a "
	  "given binary input image. Parameter:\n"
	  "   image: the input image which will be overwritten." },
	{ "findspots3d", (PyCFunction)Py_findspots3d_call, METH_VARARGS| METH_KEYWORDS, "find two hostspots in the  "
	  "given image and return the label masks of the areas. Parameters: \n"
          "   image: the image to search the two spots in\n"
	  "   rseedlimit: relative intensity (from mean image intensity) for possible hotspot seed pixel.\n"
	  "   rspotmeanlimit: relative intensity for region growing (from mean masked area image intensity)\n"
	  "   rspotbridgelimit: relative intensity for region growing (from minimum intensity of the pixels"
          "   that form the bridge between the hotspots)\n"
	  "   maxvolume: Maxium volume of a hotspot area in pixels\n"
	  "   boundinbox: tuple ((x1,y1,z1),(x2,y2,z2)) that defines a bounding box to restrict the search for spots\n"
          "   uphill: factor that describes in region growing by which factor the intensity may increase to still add "
          "to the region, should be >= 1.0\n"
	  "Returns a tuple (mask, nlabels, (seeds), (volumes)) with mask: the label mask image containing "
          "the hot areas; lower label number corresponds to the hot spot with higher intensity, nlabels: the "
          "number of labeled areas, seeds: the hotspot area seeds points as tuples, and volumes: the volumes of "
          "the segmentate regions." },
        { "shapestats", (PyCFunction)Py_get_distance_stats_call, METH_VARARGS, 
          "use a hotspot and a shape boundary given as "
          "binary mask to evaluate statistics over the distance between the hostpot and the shape. "
          "Parameters are\n"
          "   mask: the boundary mask as binary image of three dimensions\n"
          "   spot: the coordinates of the hotsport given as tuple, e.g. the one resulting from calling findspots3d."
          "Returns: a tuple with the values (distance spot-boundary gravity center,  mean distance spot-to-boundary, "
          "variation of the distance spot-to-boundary)"},
	{ NULL, NULL, 0, NULL} /* stop mark */
};




#ifdef IS_PYTHON3	



static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "csbr",
        NULL,
        0,
        csbr_methods,
        NULL,
        NULL,
        NULL,
        NULL
};

#endif

#ifdef IS_PYTHON3	
PyObject * PyInit_csbr(void)
{
	PyObject *m,*d ;

	m = PyModule_Create(&moduledef);


	/* initialize exception object */
	d=PyModule_GetDict(m) ; /* get module dictionary */
	SbrError=PyErr_NewException("csbr.error",NULL,NULL) ;
	PyDict_SetItemString(d,"error",SbrError) ;
	
	import_array();
	
	if (PyErr_Occurred()) /* something went wrong ?*/
		Py_FatalError("can't initialize module csbr") ;
        return m; 

}



#else 

PyMODINIT_FUNC initcsbr() 
{
	PyObject *m,*d ;

	m=Py_InitModule("csbr", csbr_methods);

	/* initialize exception object */
	d=PyModule_GetDict(m) ; /* get module dictionary */
	SbrError=PyErr_NewException(const_cast<char*>("csbr.error"),NULL,NULL) ;
	PyDict_SetItemString(d,"error",SbrError) ;
	
	import_array();
	
	if (PyErr_Occurred()) /* something went wrong ?*/
		Py_FatalError("can't initialize module csbr") ;
}
#endif 
