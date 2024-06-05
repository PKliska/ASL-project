#ifndef SKEWED_SIMULATION_H
#define SKEWED_SIMULATION_H
#include <stdio.h>
#include "trapeze_macros.h"
#include "faster_math_simulation.h"

typedef faster_math_simulation skewed_simulation;

skewed_simulation* new_skewed_simulation(
    size_t dimension, double size, double rho, double nu);

void init_skewed_simulation(skewed_simulation* sim,
    size_t dimension, double size, double rho, double nu);

void advance_skewed_simulation(skewed_simulation* sim,
                                 unsigned int steps, unsigned int pit,
                                 double dt);

void write_skewed_simulation(skewed_simulation* sim,
                               FILE* fp);

void destroy_skewed_simulation(skewed_simulation* sim);
void deinit_skewed_simulation(skewed_simulation* sim);
#define DO_TIME_BLOCK(b, p0, p1, d,                                           \
                      t0, t1,                                                 \
                      x_blocks, last_block_x,                                 \
                      y_blocks, last_block_y)                                 \
{                                                                             \
    DO_TRAPEZE_TOP_LEFT(b, p0, p1, d,                                         \
               t0, t1,                                                        \
               1, 0, SKEWING_BLOCK_SIZE_X, -1,                                \
               1, 0, SKEWING_BLOCK_SIZE_Y, -1);                               \
    for(int bj=1;bj<y_blocks-1;bj++){                                         \
        int y0 = bj*SKEWING_BLOCK_SIZE_Y;                                     \
        int y1 = (bj+1)*SKEWING_BLOCK_SIZE_Y;                                 \
        DO_TRAPEZE_TOP(b, p0, p1, d,                                          \
                   t0, t1,                                                    \
                   1, 0, SKEWING_BLOCK_SIZE_X, -1,                            \
                   y0, -1, y1, -1);                                           \
    }                                                                         \
    DO_TRAPEZE_TOP_RIGHT(b, p0, p1, d,                                        \
               t0, t1,                                                        \
               1, 0, SKEWING_BLOCK_SIZE_X, -1,                                \
               d-last_block_y, -1, d-1, 0);                                   \
    for(int bi=1;bi<x_blocks - 1;bi++){                                       \
        int x0 = bi*SKEWING_BLOCK_SIZE_X;                                     \
        int x1 = (bi+1)*SKEWING_BLOCK_SIZE_X;                                 \
        DO_TRAPEZE_LEFT(b, p0, p1, d,                                         \
                   t0, t1,                                                    \
                   x0, -1, x1, -1,                                            \
                   1, 0, SKEWING_BLOCK_SIZE_Y, -1);                           \
        for(int bj=1;bj<y_blocks - 1;bj++){                                   \
            int y0 = bj*SKEWING_BLOCK_SIZE_Y;                                 \
            int y1 = (bj+1)*SKEWING_BLOCK_SIZE_Y;                             \
            DO_TRAPEZE_MID(b, p0, p1, d,                                      \
                       t0, t1,                                                \
                       x0, -1, x1, -1,                                        \
                       y0, -1, y1, -1);                                       \
        }                                                                     \
        DO_TRAPEZE_RIGHT(b, p0, p1, d,                                        \
                   t0, t1,                                                    \
                   x0, -1, x1, -1,                                            \
                   d-last_block_y, -1, d-1, 0);                               \
    }                                                                         \
    DO_TRAPEZE_BOTTOM_LEFT(b, p0, p1, d,                                      \
               t0, t1,                                                        \
               d-last_block_x, -1, d-1, 0,                                    \
               1, 0, SKEWING_BLOCK_SIZE_Y, -1);                               \
    for(int bj=1;bj<y_blocks - 1;bj++){                                       \
        int y0 = bj*SKEWING_BLOCK_SIZE_Y;                                     \
        int y1 = (bj+1)*SKEWING_BLOCK_SIZE_Y;                                 \
        DO_TRAPEZE_BOTTOM(b, p0, p1, d,                                       \
                   t0, t1,                                                    \
                   d-last_block_x, -1, d-1, 0,                                \
                   y0, -1, y1, -1);                                           \
    }                                                                         \
    DO_TRAPEZE_BOTTOM_RIGHT(b, p0, p1, d,                                     \
               t0, t1,                                                        \
               d-last_block_x, -1, d-1, 0,                                    \
               d-last_block_y, -1, d-1, 0);                                   \
}
#endif
