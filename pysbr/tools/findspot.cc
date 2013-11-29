/**  
 * @file       findspot.cc
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

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <climits>

#include <queue>
#include <stdexcept>
#include <iostream>
#include <algorithm>
#include <stdio.h>
#include <map>
#include <set>
#include <vector>
#include <cmath>


#include "sbr.hh"
#include <numpy/noprefix.h>

using namespace std; 

#ifdef __GNUC__
#define UNUSED(x) x __attribute__((unused))
#else 
#define UNUSED(x) x
#endif 


#if NPY_API_VERSION < 0x00000007
#define NPY_ARRAY_WRITEABLE NPY_WRITEABLE
#define NPY_ARRAY_ALIGNED NPY_ALIGNED
#endif 

struct CandidatePixel {
        int x, y, z; 
        double value; 
        int label;  
        
        CandidatePixel():
                x(0), y(0),z(0),value(0), label(-1){}
        CandidatePixel(int _x, int _y, int _z, double v, int l):
                x(_x), y(_y), z(_z), value(v), label(l)
        {
        }
        
        void print(std::ostream& os) const {
                os << x << "," << y << "," << z <<"=(" << value << ")" << "[" << label << "]";  
        }
}; 


struct Pixel {
        int x, y, z; 
        double value; 
        
        Pixel():
                x(0), y(0),z(0),value(0){}
        Pixel(int _x, int _y, int _z, double v):
                x(_x), y(_y), z(_z), value(v)
        {
        }
        
        void print(std::ostream& os) const {
                os << x << "," << y << "," << z <<"=(" << value << ")";  
        }
        
        friend bool operator == (const Pixel& lhs, const Pixel& rhs); 

        unsigned delta (const Pixel& other) const {
                if (*this == other) 
                        return 0.0; 
                return abs(x - other.x) + abs(y - other.y) + abs(z - other.z); 
        }
}; 

enum EPivotDir { pd_x, pd_y, pd_z}; 

ostream& operator << (ostream& os, const Pixel& p) 
{
        p.print(os); 
        return os; 
}

/* Use the reverse order to remove the smallest elememts first*/
bool operator < (const Pixel& lhs, const Pixel& rhs) 
{
	return lhs.x < rhs.x || (lhs.x ==  rhs.x && ((lhs.y < rhs.y)|| (lhs.y == rhs.y && lhs.z < rhs.z))); 
}

bool operator == (const Pixel& lhs, const Pixel& rhs) 
{
	return (lhs.x ==  rhs.x && lhs.y == rhs.y && lhs.z == rhs.z); 
}


struct Path {
        double distance; 
        double distance_remain_estimate; 
        double min_intensity; 
        set<Pixel> path; 
        Pixel last_added_pixel; 
        bool end_pixel; 
        
        Path(const Pixel& start_pixel, double dhr); 
        
        bool will_loop(const Pixel& pixel)const; 
        
        Path add_pixel(const Pixel& pixel, double dhr) const;
}; 

ostream& operator << (ostream& os, const Path& p) 
{
        os << "["; 
        for_each(p.path.begin(), p.path.end(), [&os](const Pixel& pix) {
                        os << 10000*pix.z + 100*pix.y + pix.x << " "; 
                }); 
        os << "]"; 
        return os; 
}

struct path_compare {
        bool operator () (const Path& lhs, const Path& rhs) {
                return lhs.distance + lhs.distance_remain_estimate > 
                        rhs.distance  + rhs.distance_remain_estimate ; 
        }; 
}; 

Path::Path(const Pixel& start_pixel, double dhr):
        distance(0.0), 
        distance_remain_estimate(dhr), 
        last_added_pixel(start_pixel), 
        
        end_pixel(false)
{
        path.insert(start_pixel);
        min_intensity = start_pixel.value; 
}

bool Path::will_loop(const Pixel& pixel)const
{
        return (path.find(pixel) != path.end()); 
}

Path Path::add_pixel(const Pixel& pixel, double dhr)const
{
        Path retval = *this;     
        retval.path.insert(pixel);
        retval.distance += fabs(retval.last_added_pixel.value - pixel.value); 
        retval.distance_remain_estimate = dhr; 
        retval.last_added_pixel = pixel; 
        retval.end_pixel = (dhr == 0.0); 
        retval.min_intensity = min_intensity < pixel.value ? min_intensity : pixel.value; 
        return retval; 
}

template <typename T> 
float hotspot_bridge_intensity_low_by_path(const vector<CandidatePixel>& seed_pixels, 
                                           const TMirrorPythonArrayConst<T>& data, double maxpath, EPivotDir pivot_dir_sec)
{
	if (seed_pixels.size() < 2) 
		return 0.0; 



        clog << "Maxpath=" << maxpath << "\n"; 
        priority_queue<Path, vector<Path>, path_compare> path_queue; 
        
        Pixel start(seed_pixels[0].x, seed_pixels[0].y, seed_pixels[0].z, seed_pixels[0].value); 
        Pixel end(seed_pixels[1].x, seed_pixels[1].y, seed_pixels[1].z, seed_pixels[1].value); 
        

        int min_x = 0; 
        int max_x = data.get_dims()[0]; 
        int min_y = 0; 
        int max_y = data.get_dims()[1]; 
        int min_z = 0; 
        int max_z = data.get_dims()[2]; 


        double step_max = max (start.value, end.value); 

        switch (pivot_dir_sec) {
        case pd_x: 
                if (start.x > end.x) 
                        swap(start,end); 
                min_x = start.x; 
                max_x = end.x; 
                step_max = maxpath / (max_x - min_x); 
                break; 
        case pd_y: 
                if (start.y > end.y) 
                        swap(start,end); 
                min_y = start.y; 
                max_y = end.y; 
                step_max = maxpath / (max_y - min_y); 
                break; 
        case pd_z: 
                if (start.z > end.z) 
                        swap(start,end); 
                min_z = start.z; 
                max_z = end.z; 
                step_max = maxpath / (max_z - min_z); 
                break; 
                // default: n/a
        }; 

        step_max *= 0.5; 
        
        int max_path_steps = start.delta(end); 
        clog << "max-hops=" << max_path_steps << " , avg-step="<< step_max << "\n";

        path_queue.push(Path(start,  step_max * max_path_steps));         

        maxpath *= 4; 
        while (!path_queue.empty()) {

                Path p = path_queue.top(); 
                path_queue.pop(); 


                // found the shortest path
                if (p.end_pixel) {
                        clog << "Shortest path " << p.distance << " min intensity " << p.min_intensity << "\n"; 
                        return p.min_intensity;
                }


                Pixel pix = p.last_added_pixel; 
                if (p.path.size() > max_path_steps * 3u) 
                        continue; 
                
                if (pix.x - 1 >= min_x && ! (pivot_dir_sec == pd_x)) {
                        --pix.x; 
                        
                        if (!p.will_loop(pix)){
                                pix.value = data(pix.x, pix.y, pix.z); 
                                if (pix.value > 0) {
                                        auto new_p = p.add_pixel(pix, pix.delta(end) * step_max); 
                                        if (new_p.distance <= maxpath)
                                                path_queue.push(new_p); 
                                }
                        }
                        ++pix.x; 
                }
                if (pix.x + 1 <= max_x) {
                        ++pix.x; 
                        
                        if (!p.will_loop(pix)){
                                pix.value = data(pix.x, pix.y, pix.z); 
                                if (pix.value > 0) {
                                        auto new_p = p.add_pixel(pix, pix.delta(end) * step_max); 
                                        if (new_p.distance <= maxpath)
                                                path_queue.push(new_p);                                 
                                }
                        }
                        --pix.x; 
                }
                
                if (pix.y - 1 >= min_y && ! (pivot_dir_sec == pd_y)) {
                        --pix.y; 
                        
                        if (!p.will_loop(pix)){
                                pix.value = data(pix.x, pix.y, pix.z); 
                                if (pix.value > 0) {
                                        auto new_p = p.add_pixel(pix, pix.delta(end) * step_max); 
                                        if (new_p.distance <= maxpath)
                                                path_queue.push(new_p); 
                                }
                        }
                        ++pix.y; 
                }

                if (pix.y + 1 <= max_y) {
                        ++pix.y; 

                        if (!p.will_loop(pix)){
                                pix.value = data(pix.x, pix.y, pix.z); 
                                if (pix.value > 0) {
                                        auto new_p = p.add_pixel(pix, pix.delta(end) * step_max); 
                                        if (new_p.distance <= maxpath)
                                                path_queue.push(new_p); 
                                }
                        }
                        --pix.y; 
                }
                
                if (pix.z - 1 >= min_z && ! (pivot_dir_sec == pd_z)) {
                        --pix.z; 

                        if (!p.will_loop(pix)){
                                pix.value = data(pix.x, pix.y, pix.z); 
                                if (pix.value > 0) {
                                        auto new_p = p.add_pixel(pix, pix.delta(end) * step_max); 
                                        if (new_p.distance <= maxpath)
                                                path_queue.push(new_p); 
                                }
                        }
                        ++pix.z; 
                }
                if (pix.z + 1 <= max_z) {
                        ++pix.z; 

                        if (!p.will_loop(pix)){
                                pix.value = data(pix.x, pix.y, pix.z); 
                                if (pix.value > 0) {
                                        auto new_p = p.add_pixel(pix, pix.delta(end) * step_max); 
                                        if (new_p.distance <= maxpath)
                                                path_queue.push(new_p); 
                                }
                        }
                        --pix.z; 
                }
                
        }
        return 0.0; 
}

struct SFindPSpotParams {
        float seed_range_limit; 
        float spot_range_limit_from_mean; 
        float spot_range_limit_from_bridge; 
        int max_volume_per_spot; 
        float uphill_factor; 
        PyObject *boundingbox; 
}; 


ostream& operator << (ostream& os, const CandidatePixel& c) {
        c.print(os); 
        return os; 
}; 

/* Use the reverse order to remove the smallest elememts first*/
bool operator < (const CandidatePixel& lhs, const CandidatePixel& rhs) 
{
	return lhs.value < rhs.value; 
}

struct less_reverse {
	bool operator () (const CandidatePixel& lhs, const CandidatePixel& rhs) 
	{
		return !(lhs < rhs);
	}
}; 
 

/*
  Returns 1 of the center pixel is the maximum (or equal to the maximum) in neighborhood, 0 otherwise
*/
template <typename T> 
bool center_is_neighbourhood_max(const TMirrorPythonArrayConst<T>& data, CandidatePixel& p) 
{
	p.value = data(p.x, p.y, p.z); ;
	for (int nn2 = -1; nn2 < 2; ++nn2) {
		int i2 = p.z + nn2; 
		for (int nn1 = -1; nn1 < 2; ++nn1) {
			int i1 = p.y + nn1; 
                        for (int nn0 = -1; nn0 < 2; ++nn0) {
				int i0 = p.x + nn0; 
                                T val = data(i0, i1, i2); 
				if (p.value < val) {
					return false; 
                                }
			}
		}
	}
	return true; 
}

template <typename T> 
void add_6_neighbourhood(const CandidatePixel& point, const TMirrorPythonArrayConst<T>& data, priority_queue<CandidatePixel>& q, float uphill_factor)
{
        // this is the multiplyer we allow for moving up-hill 
        CandidatePixel p(point); 

        const npy_intp *dims = data.get_dims();         

        p.x = point.x - 1; 
        p.value = data(p.x, p.y, p.z); 
        T fuzzy = static_cast<T>(point.value * uphill_factor); 
	if (p.x >= 0 && p.value <= fuzzy)
		q.push(p); 
	
        p.x = point.x + 1; 
        p.value = data(p.x, p.y, p.z); 
	if (p.x < dims[0] && p.value <= fuzzy)
		q.push(p); 

        p.x = point.x;
	p.y = point.y - 1; 
        p.value = data(p.x, p.y, p.z); 
	if (p.y >= 0 && p.value <= fuzzy)
		q.push(p); 
	
	p.y = point.y + 1; 
        p.value = data(p.x, p.y, p.z); 
	if (p.y < dims[1] && p.value <= fuzzy)
		q.push(p); 
        p.y = point.y;

	p.z = point.z - 1; 
        p.value = data(p.x, p.y, p.z); 
	if (p.z >= 0 && p.value <= fuzzy)
		q.push(p); 
	
	p.z = point.z + 1; 
        p.value = data(p.x, p.y, p.z); 
	if (p.z < dims[2] && p.value <= fuzzy)
		q.push(p); 

}

template <typename T>
vector<int> grow_region(TMirrorPythonArray<Bool>& result, const vector<CandidatePixel>& p, 
                 const TMirrorPythonArrayConst<T>& data, float low_thresh, int max_volume, float uphill_factor)
{
	priority_queue<CandidatePixel> q; 
        for_each(p.begin(), p.end(), [&q](const CandidatePixel& pixel){	q.push(pixel);});  

        vector<int> maxv(p.size(), max_volume); 

        clog << "low thresh=:"<< low_thresh << "\n"; 
	while (!q.empty()) {
		CandidatePixel point = q.top(); 
		q.pop();

                if (!maxv[point.label-1]) 
                        continue; 

		unsigned char& pixel = result(point.x, point.y, point.z); 
		if (!pixel) {
                        if ( data(point.x, point.y, point.z) < low_thresh)
                                continue; 
			pixel = point.label; 
                        --maxv[point.label-1]; 
			add_6_neighbourhood(point, data, q, uphill_factor);
		}
	}
        for_each(maxv.begin(), maxv.end(), [max_volume](int& x) { x = max_volume - x;}); 

        clog << "Finished growing, volumes are:\n"; 
        for (unsigned i = 0; i < maxv.size(); ++i) 
                clog << "  " << p[i] << ": " << maxv[i] << "\n"; 
        clog << "\n";
                        
        return maxv; 
}

static void get_boundinbox (PyObject *pybox, npy_intp *box_min, npy_intp *box_max)
{
        clog << "Read boundingbox\n"; 
        if (!PyArg_ParseTuple(pybox, "(iii)(iii)", &box_min[0], &box_min[1], &box_min[2], 
                              &box_max[0], &box_max[1], &box_max[2])) {
                throw invalid_argument("bounding box not specified correctly, expect ((x1,y1,z1),(x2,y2,z2))"); 
        }
        clog << "Get boundingbox (" << box_min[0]<< ", " << box_min[1]<< ", " << box_min[2] << ")-("
             << box_max[0]<< ", " << box_max[1]<< ", " << box_max[2] << "\n"; 
}



template <typename T>
pair<float,float>  hotspot_bridge_intensity_low_x(CandidatePixel& sp0, CandidatePixel& sp1, const TMirrorPythonArrayConst<T>& data)
{
	if (sp0.x > sp1.x)
                swap(sp0, sp1); 

        int nx = sp1.x - sp0.x; 
        
        if (nx == 0) 
                throw runtime_error("findspot3d: impossible cause that the hotspots have the same coordinate in pivot direction (x)"); 
        
        float dy = float(sp1.y - sp0.y) / nx; 
        float dz = float(sp1.z - sp0.z) / nx;  
	
	float low_value = sp0.value; 
        
        float y = sp0.y; 
        float z = sp0.z;
        int x = sp0.x; 
        float delta_sum = 0; 
        float intensity = sp0.value; 
	for (int i = 0; i < nx; ++i, y += dy, z += dz) {
                int iy = static_cast<int>(floor(y + 0.5)); 
                int iz = static_cast<int>(floor(z + 0.5)); 
                auto v = data(x+i,iy,iz); 
                delta_sum += fabs(intensity - v); 
                intensity = v; 
                if (v < low_value) {
                        low_value = v; 
                }
        }
        delta_sum += fabs(intensity - sp1.value); 
	return make_pair(low_value, delta_sum); 
}				       

template <typename T>
pair<float,float>  hotspot_bridge_intensity_low_z(CandidatePixel& sp0, CandidatePixel& sp1, const TMirrorPythonArrayConst<T>& data)
{
	if (sp0.z > sp1.z)
                swap(sp0, sp1); 

        int nz = sp1.z - sp0.z; 
        
        if (nz == 0) 
		throw runtime_error("findspot3d: impossible cause that the hotspots have the same coordinate in pivot direction (z)"); 
        
        float dy = float(sp1.y - sp0.y) / nz; 
        float dx = float(sp1.x - sp0.x) / nz;  
	
	float low_value = sp0.value; 
        
        float y = sp0.y; 
        float x = sp0.x;
        int z = sp0.z; 
        float delta_sum = 0; 
        float intensity = sp0.value; 
        
	for (int i = 0; i < nz; ++i, y += dy, x += dx) {
                int iy = static_cast<int>(floor(y + 0.5)); 
                int ix = static_cast<int>(floor(x + 0.5)); 
                auto v = data(ix,iy,z+i); 
                delta_sum += fabs(intensity - v); 
                intensity = v; 
                if (v < low_value) {
                        low_value = v; 
                }
        }
        delta_sum += fabs(intensity - sp1.value); 
	return make_pair(low_value, delta_sum); 
}				       

template <typename T>
pair<float,float>  hotspot_bridge_intensity_low_y(CandidatePixel& sp0, CandidatePixel& sp1, const TMirrorPythonArrayConst<T>& data)
{
	if (sp0.y > sp1.y)
                swap(sp0, sp1); 

        int ny = sp1.y - sp0.y; 
        
        if (ny == 0) 
		throw runtime_error("findspot3d: impossible cause that the hotspots have the same coordinate in pivot direction (y)"); 
        
        float dz = float(sp1.z - sp0.z) / ny; 
        float dx = float(sp1.x - sp0.x) / ny;  
	
	float low_value = sp0.value; 

        float z = sp0.z; 
        float x = sp0.x;
        int y = sp0.y; 

        float delta_sum = 0; 
        float intensity = sp0.value; 
	for (int i = 0; i < ny; ++i, z += dz, x += dx) {
                int iz = static_cast<int>(floor(z + 0.5)); 
                int ix = static_cast<int>(floor(x + 0.5)); 
                auto v = data(ix,i+y,iz); 
                delta_sum += fabs(intensity - v); 
                intensity = v; 
                if (v < low_value) {
                        low_value = v; 
                }
        }
        delta_sum += fabs(intensity - sp1.value); 
	return make_pair(low_value, delta_sum);
}				       

template <typename T>
pair<float,float> hotspot_bridge_intensity_low(const vector<CandidatePixel>& seed_pixels, const TMirrorPythonArrayConst<T>& data, EPivotDir pd)
{
	if (seed_pixels.size() < 2) 
		return make_pair(0.0, 0.0); 
                                 
	auto sp0 = seed_pixels[0]; 
	auto sp1 = seed_pixels[1];
	
	switch (pd) {
	case pd_x: return hotspot_bridge_intensity_low_x(sp0, sp1, data); 
	case pd_y: return hotspot_bridge_intensity_low_y(sp0, sp1, data); 
	case pd_z: return hotspot_bridge_intensity_low_z(sp0, sp1, data); 
	}
        assert(0);  
}



#define COMPARE_SWAP(A,B,C) \
	if (p1. A < p2. A || (p1. A == p2. A && (p1. B < p2. B || (p1. B == p2. B && p1. C < p2. C)))){ \
		swap(p1,p2);						\
		}

void swap_pixels_if_needed(CandidatePixel& p1, CandidatePixel& p2, EPivotDir pivot_dir_prim, EPivotDir pivot_dir_sec)
{
	switch (pivot_dir_prim) {
	case pd_x: {
		if (pivot_dir_sec == pd_y) {
			COMPARE_SWAP(x,y,z); 
		}else {
			COMPARE_SWAP(x,z,y); 
		}
	}break; 
	case pd_y: {
		if (pivot_dir_sec == pd_x){
			COMPARE_SWAP(y,x,z); 
		}else {
			COMPARE_SWAP(y,z,x); 
		}
	}break; 
	case pd_z: {
		if (pivot_dir_sec == pd_x) {
			COMPARE_SWAP(z,x,y); 
		}else {
			COMPARE_SWAP(z,y,x); 
		}
	}break; 
	}
}
#undef COMPARE_SWAP

template <typename T> 
PyObject *do_findspotsT(PyArrayObject *image, const SFindPSpotParams& params) 
{
        const TMirrorPythonArrayConst<T> data(image); 

        npy_intp *dims = data.get_dims(); 
	priority_queue<CandidatePixel, vector<CandidatePixel>, less_reverse> pixels; 
        
	double max_value = 0; 
	CandidatePixel p(0,0,0,0, 0); 
        double mean = 0.0; 
        double npixels = 0;

        npy_intp box_min[3] = {dims[0] / 8, dims[1] / 8, dims[2] / 8}; 
        npy_intp box_max[3] = {dims[0] - box_min[0], dims[1]  - box_min[1], dims[2] - box_min[2]}; 

        if (params.boundingbox != NULL) 
                get_boundinbox (params.boundingbox, box_min, box_max); 
 
	for (p.z = box_min[2]; p.z < box_max[2] -1; ++p.z) {
		for (p.y = box_min[1]; p.y < box_max[1] - 1; ++p.y)
			for (p.x = box_min[0]; p.x < box_max[0] - 1; ++p.x) {
                                if (center_is_neighbourhood_max(data, p)) {
					if (max_value < p.value)
						max_value = p.value; 
					
					pixels.push(p); 
					if (pixels.size() > 100) {
						pixels.pop(); 
					}
				}
                                if (p.value > 0 ) {
                                        mean += p.value; 
                                        ++npixels; 
                                }
			}
        }
        mean /= npixels; 
	// set cutoff for interesting pixels, should be a parameter 
	max_value *= params.seed_range_limit; 
	
	// revert to order of test pixels to start at the highest intensities. 
	priority_queue<CandidatePixel> pixels2;
	while (!pixels.empty()) {
		CandidatePixel p = pixels.top(); 
		if (max_value < p.value) {
			pixels2.push(p); 
                }
		pixels.pop(); 
	}

	
	// create result image as bitmask 
	PyArrayObject* out_array = reinterpret_cast<PyArrayObject*>(PyArray_SimpleNew(3, dims, NPY_UBYTE));
        TMirrorPythonArray<unsigned char> result(out_array); 
        result.clear(); 
	


        // this is a simple region growing, doing the region growing in parallel would be better 

        bool left_cx = 0; 
        bool left_cy = 0; 
        bool left_cz = 0; 
        
        int center_x = dims[0] / 2; 
        int center_y = dims[1] / 2; 
        int center_z = dims[2] / 2; 

        vector<CandidatePixel> seed_pixels; 
	while (!pixels2.empty() && seed_pixels.size() < 2) {
		CandidatePixel p = pixels2.top(); 
		pixels2.pop(); 

		if (p.x < box_min[0] || p.x > box_max[0]) 
			continue; 

		if (p.y < box_min[1] || p.y > box_max[1]) 
			continue; 
		
		if (p.z < box_min[2] || p.z > box_max[2])
			continue; 

                p.label = seed_pixels.size() + 1; 

                // first hotspot, so set quadrant 
                if (seed_pixels.empty()) {
                        left_cx = p.x < center_x;  
                        left_cy = p.y < center_y;  
                        left_cz = p.z < center_z;  
                        clog << p << "\n"; 
                        seed_pixels.push_back(p); 
                }else {
                        // add second pixel only if at least one coordinate is in the 
                        // other side of the input area 
                        if ( (left_cx != (p.x < center_x)) ||
                             (left_cy != (p.y < center_y)) ||
                             (left_cz != (p.z < center_z))) {
                                seed_pixels.push_back(p); 
                        }
                }
        }

	float low_thresh = mean * params.spot_range_limit_from_mean; 
	
	// more then two points are never searched 
	if (seed_pixels.size() == 2)  {

		auto p1 = seed_pixels[0]; 
		auto p2 = seed_pixels[1]; 
		

		std::multimap <int,  EPivotDir> pivot_map; 

		pivot_map.insert(make_pair(abs(p1.x - p2.x), pd_x)); 
		pivot_map.insert(make_pair(abs(p1.y - p2.y), pd_y)); 
		pivot_map.insert(make_pair(abs(p1.z - p2.z), pd_z)); 

		auto ip = pivot_map.begin(); 
		++ip; 
		EPivotDir pivot_dir_sec = ip->second; 
		++ip; 
		EPivotDir pivot_dir_prim = ip->second; 




		
		swap_pixels_if_needed(p1, p2, pivot_dir_prim, pivot_dir_sec); 
		
		seed_pixels[0] = p1; 
		seed_pixels[1] = p2; 
		
		clog << "spot1 sorted = " << p1 << "\n"; 
		clog << "spot2 sorted = " << p2 << "\n"; 

		for (unsigned i = 0; i < seed_pixels.size(); ++i) 
			seed_pixels[i].label = i+1; 
		
	
		float bridge_low = 0.0; 
		
		if (params.spot_range_limit_from_bridge > 0.0 ) {
                        auto bridges = hotspot_bridge_intensity_low(seed_pixels, data, pivot_dir_prim); 

                        bridge_low = bridges.first; 
#if 0
			float new_bridge_low = hotspot_bridge_intensity_low_by_path(seed_pixels, data, bridges.second, pivot_dir_prim); 

                        if (bridge_low < new_bridge_low) 
                                bridge_low = new_bridge_low; 
#endif 
			bridge_low *= params.spot_range_limit_from_bridge; 
                }
		
		if (bridge_low > low_thresh) 
			low_thresh = bridge_low; 
	}

        auto volumes = grow_region(result, seed_pixels, data, low_thresh, 
                                   params.max_volume_per_spot, params.uphill_factor);  


        PyObject *pyresult = 0; 
        
	if (seed_pixels.size() > 1) 
		pyresult = Py_BuildValue("Oi((iii)(iii))(ii)", out_array, seed_pixels.size(), 
				     seed_pixels[0].x, seed_pixels[0].y, seed_pixels[0].z, 
				     seed_pixels[1].x, seed_pixels[1].y, seed_pixels[1].z, 
                                     volumes[0], volumes[1]
                                     ); 
	else if (seed_pixels.size() == 1) 
		pyresult = Py_BuildValue("Oi((iii))(i)", out_array, seed_pixels.size(), 
				     seed_pixels[0].x, seed_pixels[0].y, seed_pixels[0].z, 
                                     volumes[0]); 
	else 
		pyresult = Py_BuildValue("Oi((iii))(i)", out_array, 0, -1, -1, -1, 0); 
        
        return pyresult; 
}

PyObject *do_findspots(PyArrayObject *image, const SFindPSpotParams& params) 
{
	int ndims = PyArray_NDIM(image); 
	
	if (ndims != 3) 
		throw invalid_argument("expect 3 dimensional image"); 
        
	switch (image->descr->type_num) {
	case NPY_BYTE:   return do_findspotsT<signed char>(image, params);
	case NPY_UBYTE:  return do_findspotsT<unsigned char>(image, params);
	case NPY_SHORT:  return do_findspotsT<signed short>(image, params);
	case NPY_USHORT: return do_findspotsT<unsigned short>(image, params);
	case NPY_INT:    return do_findspotsT<signed int>(image, params);
	case NPY_UINT:	 return do_findspotsT<unsigned int>(image, params);
	case NPY_LONG:   return do_findspotsT<signed long>(image, params);
	case NPY_ULONG:  return do_findspotsT<unsigned long>(image, params);
	case NPY_FLOAT:  return do_findspotsT<float>(image, params);
	case NPY_DOUBLE: return do_findspotsT<double>(image, params);
	default: 
                throw invalid_argument("Input pixel type not supported."); 
	}; 
}


extern PyObject *Py_findspots3d(PyObject *UNUSED(self), PyObject *args, PyObject *kwdict)
{
	PyArrayObject *image = NULL; 
        
        SFindPSpotParams params = {
                0.7, 
                1.5, 
                1.1, 
                INT_MAX, 
                1.1, 
                NULL
        }; 

	static const char *kwlist[] = {"image", "rseedlimit","rspotmeanlimit", "rspotbridgelimit", 
                                       "maxvolume", "boundingbox", "uphill", NULL};


	if (!PyArg_ParseTupleAndKeywords(args, kwdict, "O!|fffiOf", 
                                         const_cast<char**>(kwlist), 
                                         &PyArray_Type, &image, &params.seed_range_limit, 
                                         &params.spot_range_limit_from_mean, 
                                         &params.spot_range_limit_from_bridge, 
                                         &params.max_volume_per_spot, 
                                         &params.boundingbox, 
                                         &params.uphill_factor
                                         ))
		return NULL;

	assert(image); 
        
        if (!params.boundingbox) 
                clog << "no bounding box given, limit boundaries by 1/8\n"; 
        else 
                clog << "got a bounding box parameter\n"; 
        
        if (params.boundingbox && PyTuple_Size(params.boundingbox) != 2) {
                throw invalid_argument("The bounding box must be given as a tuple ((x1,y1,z1), (x2,y2,z2))");
        }
        
	return do_findspots(image, params); 
}
