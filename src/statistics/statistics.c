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
#include <math.h>



/* ===================== */
/* Function declarations */
/* ===================== */

/* Summation */
double sum(const double *data, const size_t size);
/* Moment map generation */
double *moment(const float *data, const size_t nx, const size_t ny, const size_t nz, const unsigned int mom, const double *mom0, const double *mom1);
/* Memory release */
void free_memory(double *data);



/* ============================================= */
/* FUNCTION: Calculate the sum over a data array */
/* ============================================= */

double sum(const double *data, const size_t size)
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



/* ============================ */
/* FUNCTION: De-allocate memory */
/* ============================ */

void free_memory(double *data)
{
	if(data != NULL) free(data);
	return;
}
