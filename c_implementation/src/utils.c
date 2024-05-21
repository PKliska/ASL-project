#include "utils.h"
#include <stdlib.h>
#include <string.h>
#include <immintrin.h>
#include <assert.h>

double* zero_array(size_t len){
    size_t num_bytes = len*sizeof(double);
    double *a = (double*)(aligned_alloc(32, num_bytes));
    memset(a, 0, num_bytes);
    return a;
}
double* aligned_zero_array(size_t len){
    size_t num_bytes = len*sizeof(double);
    assert(num_bytes % 32 == 0);
    double *a = (double*)(aligned_alloc(32, num_bytes));
    memset(a, 0, num_bytes);
    return a;
}
double* copy_array(double* a, size_t len){
    size_t num_bytes = len*sizeof(double);
    double *b = _mm_malloc(num_bytes, 32);
    memcpy(b, a, num_bytes);
    return b;
}
void write_matrix(double* m, size_t nx, size_t ny, FILE* fp){
    for(size_t i=0;i<nx;i++){
        for(size_t j=0; j < (ny-1); j++){
            // print (ny-1) times with comma at the end
            // if ya wann know why ".17g" check out: https://stackoverflow.com/a/21162120/5338761
            fprintf(fp, "%.6g,", m[i*ny+j]);
        }
	// print last element withouth comma at the end (cuz there are
	// no elems after it)
	fprintf(fp, "%.6g\n", m[i*ny + ny-1]);
    }
}
