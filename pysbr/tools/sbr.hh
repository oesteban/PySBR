/** 
 * @file       sbr.hh
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

#ifndef sbr_hh
#define sbr_hh

#include <Python.h>
#if PY_MAJOR_VERSION >= 3
#define IS_PYTHON3
#endif

#define NO_IMPORT_ARRAY

#include <cassert>
#include <ostream>
#include <algorithm>
#include <vector>
#include <numpy/arrayobject.h>

template <class T>
class TMirrorPythonArray {
public: 
        TMirrorPythonArray(PyArrayObject *array); 

        const T operator ()(int x , int y, int z) const; 
        
        T& operator ()(int x , int y, int z); 

       
	npy_intp *get_dims() const; 

        void clear(); 
private: 
        int m_ndim; 
	npy_intp *m_dims; 
        std::vector<npy_intp> m_strides; 
	T *m_base;
};

template <class T>
TMirrorPythonArray<T>::TMirrorPythonArray(PyArrayObject *array):
        m_ndim(PyArray_NDIM(array)), 
        m_dims(PyArray_DIMS(array)),
        m_strides(m_ndim),
        m_base(reinterpret_cast<T *>(PyArray_DATA(array)))
{
        npy_intp *strides = PyArray_STRIDES(array); 
        std::transform(strides, strides + m_ndim, m_strides.begin(), 
                       [](npy_intp i){return i / sizeof(T);}); 
}

template <class T>
const T TMirrorPythonArray<T>::operator ()(int x , int y, int z) const
{
        assert(m_ndim == 3); 
        return m_base[x*m_strides[0] + y * m_strides[1] + z* m_strides[2]]; 
}


template <class T>
T& TMirrorPythonArray<T>::operator ()(int x , int y, int z)
{
        assert(m_ndim == 3); 
        return m_base[x*m_strides[0] + y * m_strides[1] + z* m_strides[2]]; 
}

template <class T>
void TMirrorPythonArray<T>::clear()
{
        bzero(m_base, m_dims[0] * m_dims[1] * m_dims[2] * sizeof(T)); 
}

template <class T>
npy_intp *TMirrorPythonArray<T>::get_dims() const
{
        return m_dims;
}


template <class T>
class TMirrorPythonArrayConst {
public: 
        TMirrorPythonArrayConst(PyArrayObject *array); 

        const T operator ()(int x , int y, int z) const; 
        
        T& operator ()(int x , int y, int z); 

       
	npy_intp *get_dims() const; 

        void clear(); 
private: 
        int m_ndim; 
	npy_intp *m_dims; 
        std::vector<npy_intp> m_strides; 
	const T *m_base;
};

template <class T>
TMirrorPythonArrayConst<T>::TMirrorPythonArrayConst(PyArrayObject *array):
        m_ndim(PyArray_NDIM(array)), 
        m_dims(PyArray_DIMS(array)),
        m_strides(m_ndim),
        m_base(reinterpret_cast<T *>(PyArray_DATA(array)))
{
        npy_intp *strides = PyArray_STRIDES(array); 
        std::transform(strides, strides + m_ndim, m_strides.begin(), 
                       [](npy_intp i){return i / sizeof(T);}); 
}

template <class T>
const T TMirrorPythonArrayConst<T>::operator ()(int x , int y, int z) const
{
        assert(m_ndim == 3); 
        return m_base[x*m_strides[0] + y * m_strides[1] + z* m_strides[2]]; 
}


template <class T>
void TMirrorPythonArrayConst<T>::clear()
{
        bzero(m_base, m_dims[0] * m_dims[1] * m_dims[2] * sizeof(T)); 
}

template <class T>
npy_intp *TMirrorPythonArrayConst<T>::get_dims() const
{
        return m_dims;
}


#endif 
