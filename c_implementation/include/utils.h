#ifndef UTILS_H
#define UTILS_H
#include <stdio.h>

/* Returns a pointer to an array of doubles of size *len* initialized to
 * zeroes. Caller should call free on the array after it is done.
 * */
double* zero_array(size_t len);
/* Returns a pointer to an array of doubles of initialized to first *len*
 * values of the array *a*. UB if len longer then length of a.
 * Caller should call free on the array after it is done.
 * */
double* copy_array(double* a, size_t len);
/* Outputs the nx x ny matrix stored in row first format at *m* to *fp*.
* */
void write_matrix(double* m, size_t nx, size_t ny, FILE* fp);
#endif
