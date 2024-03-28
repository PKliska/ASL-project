#include "utils.h"
#include <stdlib.h>
#include <string.h>

double* zero_array(size_t len){
    size_t num_bytes = len*sizeof(double);
    double *a = malloc(num_bytes);
    memset(a, 0, num_bytes);
    return a;
}
double* copy_array(double* a, size_t len){
    size_t num_bytes = len*sizeof(double);
    double *b = malloc(num_bytes);
    memcpy(b, a, num_bytes);
    return b;
}
void write_matrix(double* m, size_t nx, size_t ny, FILE* fp){
    for(size_t i=0;i<nx;i++){
        for(size_t j=0; j < (ny-1); j++){
	    // print (ny-1) times with comma at the end
            fprintf(fp, "%lf,", m[i*ny+j]);
        }
	// print last element withouth comma at the end (cuz there are
	// no elems after it)
	fprintf(fp, "%lf\n", m[i*ny + ny-1]);
    }
}
