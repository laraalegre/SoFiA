// ===================================================================
// This module provides time-critical and memory-critical statistical
// functions implemented in plain C99. They can be called from within
// Python using the ctypes module after compilation into a shared ob-
// ject library named statistics.so.
// ===================================================================
// Compilation: gcc -std=c99 -O3 -fPIC -shared -o statistics.so statistics.c
// ===================================================================

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#define RMS_STD 0
#define RMS_MAD 1
#define RMS_GAUSS 2

#define RMS_NEG -1
#define RMS_ALL  0
#define RMS_POS  1


#define loop_desc(I,N) for(size_t (I) = (N); (I)--;)


typedef float data_t;



/* ===================== */
/* Function declarations */
/* ===================== */

// Check for NaN
unsigned int check_nan(const data_t *data, const size_t size);
// Set mask based on threshold
void set_mask(unsigned char *mask, const data_t *data, const size_t size, const data_t threshold);
// Standard deviation
double stddev(const data_t *data, const size_t size, const size_t cadence, const int flux_range, data_t value);
// Median
data_t median(data_t *data, const size_t size);
// Median absolute deviation
data_t mad(data_t *data, const size_t size, data_t value);
// Summation
double sum(const data_t *data, const size_t size, const unsigned int mean);
// Kahan sum
double kahan_sum(const data_t *data, const size_t size, const unsigned int mean);
// Moment map generation
double *moment(const data_t *data, const size_t nx, const size_t ny, const size_t nz, const unsigned int mom, const double *mom0, const double *mom1);
// Uniform filter
void uniform_filter_1d(data_t *data, const size_t nx, const size_t ny, const size_t nz, const size_t width, const unsigned int edge_mode);
// Partial sorting of n-th element
data_t nth_element(data_t *data, const size_t size, const size_t n);
// Maximum and minimum
data_t max(const data_t *data, const size_t size);
data_t min(const data_t *data, const size_t size);
// Memory release
void free_memory(double *data);
// Check native byte order of machine
inline unsigned int native_byte_order(void);
// Swap byte order
inline void swap_byte_order_32(float *value);
inline void swap_byte_order_64(double *value);
// Check for NaN
inline unsigned int is_nan(const data_t *value);



// -------------------------
// Check if any NaN in array
// -------------------------

unsigned int check_nan(const data_t *data, const size_t size)
{
	const data_t *ptr = data + size;
	while(ptr --> data) if(is_nan(ptr)) return 1;
	return 0;
}



// ----------------------------------
// Maximum and minimum value in array
// ----------------------------------

data_t max(const data_t *data, const size_t size)
{
	data_t result = NAN;
	const data_t *ptr = data + size;
	
	while(ptr --> data)
	{
		if(is_nan(ptr)) continue;
		if(is_nan(&result) || *ptr > result) result = *ptr;
	}
	
	return result;
}

data_t min(const data_t *data, const size_t size)
{
	data_t result = NAN;
	const data_t *ptr = data + size;
	
	while(ptr --> data)
	{
		if(is_nan(ptr)) continue;
		if(is_nan(&result) || *ptr < result) result = *ptr;
	}
	
	return result;
}



// ---------------------------
// Summation of array elements
// ---------------------------

double sum(const data_t *data, const size_t size, const unsigned int mean)
{
	const data_t *ptr = data + size;
	double result = 0.0;
	size_t counter = 0;
	
	while(ptr --> data)
	{
		if(!is_nan(ptr))
		{
			result += *ptr;
			++counter;
		}
	}
	
	if(mean && counter) return result / (double)counter;
	return result;
}



// ---------------------------------
// Kahan summation of array elements
// ---------------------------------

double kahan_sum(const data_t *data, const size_t size, const unsigned int mean)
{
	const data_t *ptr = data + size;
	double result = 0.0;
	double error = 0.0;
	size_t counter = 0;
	
	while(ptr --> data)
	{
		if(!is_nan(ptr))
		{
			double y = *ptr - error;
			double t = result + y;
			error = (t - result) - y;
			result = t;
			++counter;
		}
	}
	
	if(mean && counter) return result / (double)counter;
	return result;
}



// --------------------------------------
// Standard deviation about value or mean
// --------------------------------------

double stddev(const data_t *data, const size_t size, const size_t cadence, const int flux_range, data_t value)
{
	double result = 0.0;
	size_t counter = 0;
	const data_t *ptr = data + size;
	
	// Calculate mean if no value specified
	if(is_nan(&value)) value = sum(data, size, 1);
	
	while(data < ptr)
	{
		ptr -= cadence;
		if(is_nan(ptr)) continue;
		
		if(*ptr < 0.0)
		{
			if(flux_range <= RMS_ALL)
			{
				result += (*ptr - value) * (*ptr - value);
				++counter;
			}
		}
		else
		{
			if(flux_range >= RMS_ALL)
			{
				result += (*ptr - value) * (*ptr - value);
				++counter;
			}
		}
	}
	
	return sqrt(result / counter);
}



// ---------------
// Median of array
// ---------------

// WARNING: Not NaN-safe!
data_t median(data_t *data, const size_t size)
{
	/* Exact median */
	const size_t n = size / 2;
	const data_t value = nth_element(data, size, n);
	if(size & 1U) return value;
	return (value + max(data, n)) / 2.0;
	
	/* Approximate median (marginally faster) */
	/*return nth_element(data, size, size / 2);*/
}



// -------------------------
// Median absolute deviation
// -------------------------

data_t mad(data_t *data, const size_t size, data_t value)
{
	data_t *ptr = data + size;
	if(is_nan(&value)) value = median(data, size);
	while(ptr --> data) *ptr = fabs(*ptr - value);
	return 1.4826 * median(data, size);
}



// ------------------------------
// N-th smallest element in array
// ------------------------------

// WARNING: Not NaN-safe!
data_t nth_element(data_t *data, const size_t size, const size_t n)
{
	data_t *l = data;
	data_t *m = data + size - 1;
	data_t *ptr = data + n;
	
	while(l < m)
	{
		data_t value = *ptr;
		data_t *i = l;
		data_t *j = m;
		
		do
		{
			while(*i < value) ++i;
			while(value < *j) --j;
			
			if(i <= j)
			{
				data_t tmp = *i;
				*i = *j;
				*j = tmp;
				++i;
				--j;
			}
		} while(i <= j);
		
		if(j < ptr) l = i;
		if(ptr < i) m = j;
	}
	
	return *ptr;
}



// ---------------------------
// Set mask based on threshold
// ---------------------------

void set_mask(unsigned char *mask, const data_t *data, const size_t size, const data_t threshold)
{
	const data_t *ptr_data = data + size;
	unsigned char *ptr_mask = mask + size;
	
	while(ptr_data --> data)
	{
		--ptr_mask;
		if(fabs(*ptr_data) >= threshold) *ptr_mask = 1U;
	}
	
	return;
}



/* =============================================== */
/* FUNCTION: Calculate the moment 0 of a 3-D array */
/* =============================================== */

// ALERT: This needs to be updated and sped up!
double *moment(const data_t *data, const size_t nx, const size_t ny, const size_t nz, const unsigned int mom, const double *mom0, const double *mom1)
{
	double *moment_map = (double*)calloc(nx * ny, sizeof(double));
	
	if(moment_map == NULL)
	{
		fprintf(stderr, "ERROR: Failed to allocate memory for moment map.\n");
		exit(1);
	}
	
	loop_desc(x, nx)
	{
		loop_desc(y, ny)
		{
			unsigned int flag = 0;
			size_t index_2d = y + ny * x;
			
			loop_desc(z, nz)
			{
				size_t index_3d = z + nz * index_2d;
				
				if(!isnan(data[index_3d]))
				{
					flag = 1;
					double offset;
					
					switch(mom)
					{
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
			
			if(flag)
			{
				switch(mom)
				{
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
			else
			{
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
void uniform_filter_1d(data_t *data, const size_t nx, const size_t ny, const size_t nz, const size_t width, const unsigned int edge_mode)
{
	size_t radius = width / 2;
	
	/* Allocate memory for spectrum */
	data_t *spectrum = (data_t*)malloc(nz * sizeof(data_t));
	
	if(spectrum == NULL)
	{
		fprintf(stderr, "ERROR: Failed to allocate memory for convolution.\n");
		exit(1);
	}
	
	/* Loop over the image plane */
	loop_desc(x, nx)
	{
		loop_desc(y, ny)
		{
			/* Determine pixel index */
			size_t index_2d = y + ny * x;
			
			/* Copy data into spectrum */
			loop_desc(z, nz)
			{
				size_t index_3d = z + nz * index_2d;
				spectrum[z] = data[index_3d];
			}
			
			data_t value = 0.0;
			
			/* Apply filter */
			for(size_t z = radius; z < nz - radius; ++z)
			{
				size_t index_3d = z + nz * index_2d;
				value += spectrum[z] / (data_t)width;
				
				for(size_t n = 1; n <= radius; ++n)
				{
					value += spectrum[z - n] / (data_t)width;
					value += spectrum[z + n] / (data_t)width;
				}
				
				data[index_3d] = value;
			}
			
			/* Mask unfiltered channels with NaN */
			loop_desc(n, radius)
			{
				size_t index_3d = n + nz * index_2d;
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



/* ========================= */
/* FUNCTION: Swap byte order */
/* ========================= */

inline void swap_byte_order_32(float *value)
{
	uint32_t tmp;
	memcpy(&tmp, value, 4);
	tmp = __builtin_bswap32(tmp);
	memcpy(value, &tmp, 4);
	return;
}

inline void swap_byte_order_64(double *value)
{
	uint64_t tmp;
	memcpy(&tmp, value, 8);
	tmp = __builtin_bswap64(tmp);
	memcpy(value, &tmp, 8);
	return;
}


/* ============ */
/* Check if NaN */
/* ============ */

inline unsigned int is_nan(const data_t *value)
{
	return *value != *value;
}
