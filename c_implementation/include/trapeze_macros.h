#include "utils.h"

#define TRAPEZE_LOOP_UTD for(int tx=tx0;tx<tx1;tx++)
#define TRAPEZE_LOOP_DTU for(int tx=tx1-1;tx>=tx0;tx--)
#define TRAPEZE_LOOP_LTR for(int ty=ty0;ty<ty1;ty++)
#define TRAPEZE_LOOP_RTL for(int ty=ty1-1;ty>=ty0;ty--)
#define AFTER_LEFT new[tx*d + 0] = new[tx*d + 1]
#define AFTER_RIGHT new[tx*d + d - 1] = new[tx*d + d - 2]
#define AFTER_TOP memcpy(new + 0*d+ty0, new + 1*d+ty0, (ty1-ty0)*sizeof(double))
#define AFTER_BOTTOM memset(new + (d-1)*d+ty0, 0, (ty1-ty0)*sizeof(double))
#define EMPTY_AFTER 

#define DO_TRAPEZE_GENERAL(loop1, loop2,                                      \
                           after_loop1, after_loop2,                          \
                           b, p0, p1, d,                                      \
                           t0, t1,                                            \
                           x0, dx0, x1, dx1,                                  \
                           y0, dy0, y1, dy1)                                  \
{                                                                             \
    double *restrict old, *restrict new;                                      \
    if(t0 % 2 == 0){                                                          \
        new = p0;                                                             \
        old = p1;                                                             \
    }else{                                                                    \
        new = p1;                                                             \
        old = p0;                                                             \
    }                                                                         \
    int tx0 = x0, tx1 = x1;                                                   \
    int ty0 = y0, ty1 = y1;                                                   \
    for(int t=t0;t<t1;t++){                                                   \
        loop1{                                                                \
            loop2{                                                            \
                new[tx*d + ty] = (old[(tx+1)*d + ty] + old[(tx-1)*d + ty]     \
                                + old[tx*d + (ty+1)] + old[tx*d + (ty-1)]     \
                                - b[tx*d+ty]) * 0.25;                         \
            }                                                                 \
            after_loop2;                                                      \
        }                                                                     \
        after_loop1;                                                          \
        tx0 += dx0; tx1 += dx1; ty0 += dy0; ty1 += dy1;                       \
        SWAP(double*, new, old);                                              \
    }                                                                         \
}

#define DO_TRAPEZE_TOP_LEFT(b, p0, p1, d,                                     \
                           t0, t1,                                            \
                           x0, dx0, x1, dx1,                                  \
                           y0, dy0, y1, dy1)                                  \
DO_TRAPEZE_GENERAL(TRAPEZE_LOOP_DTU, TRAPEZE_LOOP_RTL,                        \
                   AFTER_TOP, AFTER_LEFT,                                     \
                   b, p0, p1, d,                                              \
                   t0, t1,                                                    \
                   x0, dx0, x1, dx1,                                          \
                   y0, dy0, y1, dy1)                                          \

#define DO_TRAPEZE_TOP(b, p0, p1, d,                                          \
                           t0, t1,                                            \
                           x0, dx0, x1, dx1,                                  \
                           y0, dy0, y1, dy1)                                  \
DO_TRAPEZE_GENERAL(TRAPEZE_LOOP_DTU, TRAPEZE_LOOP_LTR,                        \
                   AFTER_TOP, EMPTY_AFTER,                                    \
                   b, p0, p1, d,                                              \
                   t0, t1,                                                    \
                   x0, dx0, x1, dx1,                                          \
                   y0, dy0, y1, dy1)                                          \


#define DO_TRAPEZE_TOP_RIGHT(b, p0, p1, d,                                    \
                           t0, t1,                                            \
                           x0, dx0, x1, dx1,                                  \
                           y0, dy0, y1, dy1)                                  \
DO_TRAPEZE_GENERAL(TRAPEZE_LOOP_DTU, TRAPEZE_LOOP_LTR,                        \
                   AFTER_TOP, AFTER_RIGHT,                                    \
                   b, p0, p1, d,                                              \
                   t0, t1,                                                    \
                   x0, dx0, x1, dx1,                                          \
                   y0, dy0, y1, dy1)                                          \

#define DO_TRAPEZE_LEFT(b, p0, p1, d,                                         \
                           t0, t1,                                            \
                           x0, dx0, x1, dx1,                                  \
                           y0, dy0, y1, dy1)                                  \
DO_TRAPEZE_GENERAL(TRAPEZE_LOOP_UTD, TRAPEZE_LOOP_RTL,                        \
                   EMPTY_AFTER, AFTER_LEFT,                                   \
                   b, p0, p1, d,                                              \
                   t0, t1,                                                    \
                   x0, dx0, x1, dx1,                                          \
                   y0, dy0, y1, dy1)                                          \


#define DO_TRAPEZE_MID(b, p0, p1, d,                                          \
                           t0, t1,                                            \
                           x0, dx0, x1, dx1,                                  \
                           y0, dy0, y1, dy1)                                  \
DO_TRAPEZE_GENERAL(TRAPEZE_LOOP_UTD, TRAPEZE_LOOP_LTR,                        \
                   EMPTY_AFTER, EMPTY_AFTER,                                  \
                   b, p0, p1, d,                                              \
                   t0, t1,                                                    \
                   x0, dx0, x1, dx1,                                          \
                   y0, dy0, y1, dy1)                                          \

#define DO_TRAPEZE_RIGHT(b, p0, p1, d,                                        \
                           t0, t1,                                            \
                           x0, dx0, x1, dx1,                                  \
                           y0, dy0, y1, dy1)                                  \
DO_TRAPEZE_GENERAL(TRAPEZE_LOOP_UTD, TRAPEZE_LOOP_LTR,                        \
                   EMPTY_AFTER, AFTER_RIGHT,                                  \
                   b, p0, p1, d,                                              \
                   t0, t1,                                                    \
                   x0, dx0, x1, dx1,                                          \
                   y0, dy0, y1, dy1)                                          \

#define DO_TRAPEZE_BOTTOM_LEFT(b, p0, p1, d,                                  \
                           t0, t1,                                            \
                           x0, dx0, x1, dx1,                                  \
                           y0, dy0, y1, dy1)                                  \
DO_TRAPEZE_GENERAL(TRAPEZE_LOOP_UTD, TRAPEZE_LOOP_RTL,                        \
                   AFTER_BOTTOM, AFTER_LEFT,                                   \
                   b, p0, p1, d,                                              \
                   t0, t1,                                                    \
                   x0, dx0, x1, dx1,                                          \
                   y0, dy0, y1, dy1)                                          \


#define DO_TRAPEZE_BOTTOM(b, p0, p1, d,                                       \
                           t0, t1,                                            \
                           x0, dx0, x1, dx1,                                  \
                           y0, dy0, y1, dy1)                                  \
DO_TRAPEZE_GENERAL(TRAPEZE_LOOP_UTD, TRAPEZE_LOOP_LTR,                        \
                   AFTER_BOTTOM, EMPTY_AFTER,                                 \
                   b, p0, p1, d,                                              \
                   t0, t1,                                                    \
                   x0, dx0, x1, dx1,                                          \
                   y0, dy0, y1, dy1)                                          \

#define DO_TRAPEZE_BOTTOM_RIGHT(b, p0, p1, d,                                 \
                           t0, t1,                                            \
                           x0, dx0, x1, dx1,                                  \
                           y0, dy0, y1, dy1)                                  \
DO_TRAPEZE_GENERAL(TRAPEZE_LOOP_UTD, TRAPEZE_LOOP_LTR,                        \
                   AFTER_BOTTOM, AFTER_RIGHT,                                 \
                   b, p0, p1, d,                                              \
                   t0, t1,                                                    \
                   x0, dx0, x1, dx1,                                          \
                   y0, dy0, y1, dy1)                                          \

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
