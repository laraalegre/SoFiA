/* =================================================================== */
/* This module provides time-critical and memory-critical statistical  */
/* functions implemented in plain C99. They can be called from within  */
/* Python using the ctypes module after compilation into a shared ob-  */
/* ject library named statistics.so.                                   */
/* =================================================================== */
/* Compilation: gcc -O3 -fPIC -shared -o statistics.so statistics.c    */
/* =================================================================== */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#define RMS_STD 0
#define RMS_MAD 1
#define RMS_GAUSS 2

#define RMS_NEG 0
#define RMS_ALL 1
#define RMS_POS 2



/* ===================== */
/* Function declarations */
/* ===================== */

/* Check for NaN */
unsigned int check_nan(const float *data, const size_t size, const unsigned int byte_order);
/* Replace NaN */
void replace_nan(float *data, const float *mask, const size_t size, float value, const unsigned int byte_order_data, const unsigned int byte_order_mask);
/* Set mask based on threshold */
void set_mask(unsigned char *mask, const float *data, const size_t size, const float threshold, const unsigned int byte_order);
/* Measure RMS */
double standard_deviation(const float *data, const size_t size, const size_t cadence, const unsigned int flux_range, const unsigned int byte_order);
/* Summation */
double sum(const float *data, const size_t size);
/* Moment map generation */
double *moment(const float *data, const size_t nx, const size_t ny, const size_t nz, const unsigned int mom, const double *mom0, const double *mom1);
/* Uniform filter */
void uniform_filter_1d(float *data, const size_t nx, const size_t ny, const size_t nz, const size_t width, const unsigned int edge_mode);
/* Memory release */
void free_memory(double *data);
/* Check native byte order of machine */
inline unsigned int native_byte_order(void);
/* Swap byte order */
inline void swap_byte_order(float *value);



/* =================================== */
/* FUNCTION: Check for presence of NaN */
/* =================================== */

unsigned int check_nan(const float *data, const size_t size, const unsigned int byte_order)
{
	const unsigned int swap_needed = (byte_order != native_byte_order());
	
	for(size_t i = size; i--;)
	{
		float value = data[i];
		if(swap_needed) swap_byte_order(&value);
		if(isnan(value)) return 1;
	}
	
	return 0;
}



/* ================================ */
/* FUNCTION: Replace NaN with value */
/* ================================ */

void replace_nan(float *data, const float *mask, const size_t size, float value, const unsigned int byte_order_data, const unsigned int byte_order_mask)
{
	const unsigned int swap_needed_data = (byte_order_data != native_byte_order());
	const unsigned int swap_needed_mask = (byte_order_mask != native_byte_order());
	
	for(size_t i = size; i--;)
	{
		float mask_value = mask[i];
		if(swap_needed_mask) swap_byte_order(&mask_value);
		
		if(isnan(mask_value))
		{
			if(swap_needed_data) swap_byte_order(&value);
			data[i] = value;
		}
	}
	
	return;
}



/* ===================================== */
/* FUNCTION: Set mask based on threshold */
/* ===================================== */

void set_mask(unsigned char *mask, const float *data, const size_t size, const float threshold, const unsigned int byte_order)
{
	const unsigned int swap_needed = (byte_order != native_byte_order());
	
	for(size_t i = size; i--;)
	{
		float value = data[i];
		if(swap_needed) swap_byte_order(&value);
		if(fabs(value) >= threshold) mask[i] = 1;
	}
	
	return;
}



/* ==================================== */
/* FUNCTION: Measure standard deviation */
/* ==================================== */

double standard_deviation(const float *data, const size_t size, const size_t cadence, const unsigned int flux_range, const unsigned int byte_order)
{
	double result = 0.0;
	size_t counter = 0;
	const unsigned int swap_needed = (byte_order != native_byte_order());
	
	for(size_t i = size - size % cadence; i -= cadence;)
	{
		float value = data[i];
		if(swap_needed) swap_byte_order(&value);
		if(isnan(value)) continue;
		
		if(value < 0.0)
		{
			if(flux_range <= RMS_ALL)
			{
				result += value * value;
				++counter;
			}
		}
		else
		{
			if(flux_range >= RMS_ALL)
			{
				result += value * value;
				++counter;
			}
		}
	}
	
	return sqrt(result / counter);
}



/* ============================================= */
/* FUNCTION: Calculate the sum over a data array */
/* ============================================= */

double sum(const float *data, const size_t size)
{
	double sum = 0.0;
	unsigned int flag = 0;
	size_t i = size;
	
	for(;i--;) {
		if(!isnan(data[i])) {
			sum += data[i];
			flag = 1;
		}
	}
	
	if(flag) return sum;
	return NAN;
}



/* =============================================== */
/* FUNCTION: Calculate the moment 0 of a 3-D array */
/* =============================================== */

double *moment(const float *data, const size_t nx, const size_t ny, const size_t nz, const unsigned int mom, const double *mom0, const double *mom1)
{
	size_t x = nx;
	size_t y = ny;
	size_t z = nz;
	size_t index_2d;
	size_t index_3d;
	double offset;
	unsigned int flag;
	
	double *moment_map = (double*)calloc(nx * ny, sizeof(double));
	
	if(moment_map == NULL) {
		fprintf(stderr, "ERROR: Failed to allocate memory for moment map.\n");
		exit(1);
	}
	
	for(;x--;) {
		for(;y--;) {
			flag = 0;
			index_2d = y + ny * x;
			
			for(;z--;) {
				index_3d = z + nz * index_2d;
				
				if(!isnan(data[index_3d])) {
					flag = 1;
					switch(mom) {
						case 2:
							offset = mom1[index_2d] - z;
							moment_map[index_2d] += offset * offset * data[index_3d];
							break;
						case 1:
							moment_map[index_2d] += data[index_3d] * z;
							break;
						default:
							moment_map[index_2d] += data[index_3d];
							break;
					}
				}
			}
			
			if(flag) {
				switch(mom) {
					case 2:
						moment_map[index_2d] = sqrt(moment_map[index_2d] / mom0[index_2d]);
						break;
					case 1:
						moment_map[index_2d] /= mom0[index_2d];
						break;
					default:
						break;
				}
			}
			else {
				moment_map[index_2d] = NAN;
			}
		}
	}
	
	return moment_map;
}



/* ===================================== */
/* FUNCTION: Uniform filter along z axis */
/* ===================================== */

/* WARNING: This function is untested and cannot handle NaN and inf! */
void uniform_filter_1d(float *data, const size_t nx, const size_t ny, const size_t nz, const size_t width, const unsigned int edge_mode)
{
	size_t x, y, z, n;
	size_t index_2d, index_3d;
	size_t radius = width / 2;
	float value;
	
	/* Allocate memory for spectrum */
	float *spectrum = (float*)malloc(nz * sizeof(float));
	
	if(spectrum == NULL)
	{
		fprintf(stderr, "ERROR: Failed to allocate memory for convolution.\n");
		exit(1);
	}
	
	/* Loop over the image plane */
	for(x = nx; x--;)
	{
		for(y = ny; y--;)
		{
			/* Determine pixel index */
			index_2d = y + ny * z;
			
			/* Copy data into spectrum */
			for(z = nz; z--;)
			{
				index_3d = x + nx * index_2d;
				spectrum[z] = data[index_3d];
			}
			
			value = 0.0;
			
			/* Apply filter */
			for(z = radius; z < nz - radius; ++z)
			{
				index_3d = z + nz * index_2d;
				value += spectrum[z] / (float)width;
				
				for(n = 1; n <= radius; ++n)
				{
					value += spectrum[z - n] / (float)width;
					value += spectrum[z + n] / (float)width;
				}
				
				data[index_3d] = value;
			}
			
			/* Mask unfiltered channels with NaN */
			for(n = 0; n < radius; ++n)
			{
				index_3d = n + nz * index_2d;
				data[index_3d] = NAN;
				index_3d = (nz - n - 1) + nz * index_2d;
				data[index_3d] = NAN;
			}
		}
	}
	
	/* Release memory again */
	free(spectrum);
	
	return;
}



/* ============================ */
/* FUNCTION: De-allocate memory */
/* ============================ */

void free_memory(double *data)
{
	if(data != NULL) free(data);
	return;
}



/* ============================================ */
/* FUNCTION: Check native byte order of machine */
/* ============================================ */

inline unsigned int native_byte_order(void)
{
	// Returns 0 on little-endian and 1 on big-endian machines
	// NOTE: This check will only work on systems
	//       where sizeof(long) > sizeof(char) = 1
	
	long n = 1;
	return *(char *)&n != 1;
}



/* ================================== */
/* FUNCTION: Swap byte order of float */
/* ================================== */

inline void swap_byte_order(float *value)
{
	uint32_t tmp;
	memcpy(&tmp, value, 4);
	tmp = __builtin_bswap32(tmp);
	memcpy(value, &tmp, 4);
	return;
}
