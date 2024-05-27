#ifndef UTILS_H
#define UTILS_H
#include <stdio.h>
#ifdef NDEBUG
    #define DPRINTF(format, ...) 
#else
    #define DPRINTF(format, ...)         \
{                                        \
    fprintf(stderr, format, __VA_ARGS__);\
    fflush(stderr);                      \
}
#endif


struct trapeze{
    unsigned t0, t1;
    unsigned x0, x1, y0, y1;
    int dx0, dx1, dy0, dy1;
};

/* Returns a pointer to an array of doubles of size *len* initialized to
 * zeroes. Caller should call free on the array after it is done.
 * */
double* zero_array(size_t len);
/* Returns a pointer to an array of doubles of size *len* initialized to
 * zeroes. Caller should call free on the array after it is done.
 * Array is 32-byte aligned.
 * */
double* aligned_zero_array(size_t len);
/* Returns a pointer to an array of doubles of initialized to first *len*
 * values of the array *a*. UB if len longer then length of a.
 * Caller should call free on the array after it is done.
 * */
double* copy_array(double* a, size_t len);
/* Outputs the nx x ny matrix stored in row first format at *m* to *fp*.
* */
void write_matrix(double* m, size_t nx, size_t ny, FILE* fp);

#define SWAP(T, a, b){ \
    T t = (a); \
    (a) = (b); \
    (b) = t; \
}
#endif
