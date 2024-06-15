#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include "uskewed_simulation.h"
#include "faster_math_simulation.h"
#include "utils.h"
#include "quarter_trapeze_macros.h"

#ifndef USKEWING_BLOCK_SIZE_X
#error "USKEWING_BLOCK_SIZE_X was not defined, setting to default (36)"
#endif
#ifndef USKEWING_BLOCK_SIZE_Y
#error "USKEWING_BLOCK_SIZE_Y was not defined, setting to default (36)"
#endif
#ifndef USKEWING_TIMESTEPS
#error "USKEWING_TIMESTEPS was not defined, setting to default (10)"
#endif
#if USKEWING_BLOCK_SIZE_X <= USKEWING_TIMESTEPS \
 || USKEWING_BLOCK_SIZE_Y <= USKEWING_TIMESTEPS
#error "USKEWING_BLOCK_SIZE must be greater than USKEWING_TIMESTEPS"
#endif



static const struct simulation_vtable_ uskewed_SIMULATION_VTABLE[] = {{
    .advance=(void (*)(struct simulation *, unsigned int, unsigned int, double))advance_uskewed_simulation,
    .write=(void (*)(const struct simulation *, FILE *))write_uskewed_simulation,
    .destroy=(void (*)(struct simulation *))destroy_uskewed_simulation
}};


uskewed_simulation* new_uskewed_simulation(
    size_t dimension, double size, double rho, double nu){
    uskewed_simulation* sim = malloc(sizeof(*sim));
    init_uskewed_simulation(sim, dimension, size, rho, nu);
    return sim;
}

void init_uskewed_simulation(uskewed_simulation* sim,
    size_t dimension, double size, double rho, double nu){
    init_faster_math_simulation(sim, dimension, size, rho, nu);
    assert((dimension + USKEWING_BLOCK_SIZE_X - 1) / USKEWING_BLOCK_SIZE_X >= 2);
    assert((dimension + USKEWING_BLOCK_SIZE_Y - 1) / USKEWING_BLOCK_SIZE_Y >= 2);
    assert(dimension % 8 == 0);
    assert(USKEWING_BLOCK_SIZE_Y % 8 == 0);
    sim->base.vtable_ = uskewed_SIMULATION_VTABLE;
}


static double sq(const double x){
    return x*x;
}

#define UCOMPUTE_A(old, d, x, y) old[((x)+1)*(d) + (y)] + old[(x)*(d) + (y)+1]
#define UCOMPUTE_B(old, d, x, y) old[((x)-1)*(d) + (y)] + old[(x)*(d) + (y)-1]
#define UCOMPUTE_C(a, b) (a) + (b)
#define UCOMPUTE_D(c, b, d, x, y) (c)*0.25 - (b)[(x)*(d)+(y)]

__attribute__((always_inline))
static inline void uskewed_trapeze_mid(double *restrict b, double *restrict p0, double *restrict p1,
                                size_t d, unsigned t0, unsigned t1,
                                int x0, int x1,
                                int y0, int y1){
    double *restrict old, *restrict new;
    if(t0 % 2 == 0){
        new = p0;
        old = p1;
    }else{
        new = p1;
        old = p0;
    }
    int tx0 = x0, tx1 = x1;
    int ty0 = y0, ty1 = y1;
    for(unsigned t=t0;t<(unsigned)t1;t++){
        for(int tx=tx0;tx<tx1;tx+=8){
        for(int ty=ty0;ty<ty1;ty+=4){
		#pragma clang loop unroll(full)
		for(int i=0;i<8;i++){
                double a0 = UCOMPUTE_A(old, d, tx+i, ty+0);
                double b0 = UCOMPUTE_B(old, d, tx+i, ty+0);
                double a1 = UCOMPUTE_A(old, d, tx+i, ty+1);
                double b1 = UCOMPUTE_B(old, d, tx+i, ty+1);
                double a2 = UCOMPUTE_A(old, d, tx+i, ty+2);
                double b2 = UCOMPUTE_B(old, d, tx+i, ty+2);
                double a3 = UCOMPUTE_A(old, d, tx+i, ty+3);
                double b3 = UCOMPUTE_B(old, d, tx+i, ty+3);
                double c0 = UCOMPUTE_C(a0, b0);
                double c1 = UCOMPUTE_C(a1, b1);
                double c2 = UCOMPUTE_C(a2, b2);
                double c3 = UCOMPUTE_C(a3, b3);
                new[(tx+i)*d + ty+0] = UCOMPUTE_D(c0, b, d, tx+i, ty+0);
                new[(tx+i)*d + ty+1] = UCOMPUTE_D(c1, b, d, tx+i, ty+1);
                new[(tx+i)*d + ty+2] = UCOMPUTE_D(c2, b, d, tx+i, ty+2);
                new[(tx+i)*d + ty+3] = UCOMPUTE_D(c3, b, d, tx+i, ty+3);
		}
            }
        }
        tx0--; tx1--; ty0--; ty1--;
        SWAP(double*, new, old);
    }
}

static void uskewed_time_block(double *restrict b, double *restrict p0, double *restrict p1,
                               size_t d, unsigned t0, unsigned t1,
                               int x_blocks, int last_block_x,
                               int y_blocks, int last_block_y){
    DO_TRAPEZE_TOP_LEFT(b, p0, p1, d,
               t0, t1,
               1, 0, USKEWING_BLOCK_SIZE_X, -1,
               1, 0, USKEWING_BLOCK_SIZE_Y, -1);
    for(int bj=1;bj<y_blocks-1;bj++){
        int y0 = bj*USKEWING_BLOCK_SIZE_Y;
        int y1 = (bj+1)*USKEWING_BLOCK_SIZE_Y;
        DO_TRAPEZE_TOP(b, p0, p1, d,
                   t0, t1,
                   1, 0, USKEWING_BLOCK_SIZE_X, -1,
                   y0, -1, y1, -1);
    }
    DO_TRAPEZE_TOP_RIGHT(b, p0, p1, d,
               t0, t1,
               1, 0, USKEWING_BLOCK_SIZE_X, -1,
               d-last_block_y, -1, d-1, 0);
    for(int bi=1;bi<x_blocks - 1;bi++){
        int x0 = bi*USKEWING_BLOCK_SIZE_X;
        int x1 = (bi+1)*USKEWING_BLOCK_SIZE_X;
        DO_TRAPEZE_LEFT(b, p0, p1, d,
                   t0, t1,
                   x0, -1, x1, -1,
                   1, 0, USKEWING_BLOCK_SIZE_Y, -1);
        for(int bj=1;bj<y_blocks - 1;bj++){
            int y0 = bj*USKEWING_BLOCK_SIZE_Y;
            int y1 = (bj+1)*USKEWING_BLOCK_SIZE_Y;


            uskewed_trapeze_mid(b, p0, p1, d,
                       t0, t1,
                       x0, x1,
                       y0, y1);
        }
        DO_TRAPEZE_RIGHT(b, p0, p1, d,
                   t0, t1,
                   x0, -1, x1, -1,
                   d-last_block_y, -1, d-1, 0);
    }
    DO_TRAPEZE_BOTTOM_LEFT(b, p0, p1, d,
               t0, t1,
               d-last_block_x, -1, d-1, 0,
               1, 0, USKEWING_BLOCK_SIZE_Y, -1);
    for(int bj=1;bj<y_blocks - 1;bj++){
        int y0 = bj*USKEWING_BLOCK_SIZE_Y;
        int y1 = (bj+1)*USKEWING_BLOCK_SIZE_Y;
        DO_TRAPEZE_BOTTOM(b, p0, p1, d,
                   t0, t1,
                   d-last_block_x, -1, d-1, 0,
                   y0, -1, y1, -1);
    }
    DO_TRAPEZE_BOTTOM_RIGHT(b, p0, p1, d,
               t0, t1,
               d-last_block_x, -1, d-1, 0,
               d-last_block_y, -1, d-1, 0);
}
static void uskewed_pressure_poisson(uskewed_simulation* sim,
                 unsigned int pit){
    const size_t d = sim->d;
    double *restrict b = sim->b;
    double *restrict p = sim->p;
    double *restrict pn = sim->pn;
    const int time_blocks = (pit+USKEWING_TIMESTEPS-1) / USKEWING_TIMESTEPS;
    const int last_time = (pit-1) % USKEWING_TIMESTEPS + 1;
    const int x_blocks = (d+USKEWING_BLOCK_SIZE_X-1) / USKEWING_BLOCK_SIZE_X;
    const int y_blocks = (d+USKEWING_BLOCK_SIZE_Y-1) / USKEWING_BLOCK_SIZE_Y;
    const int last_block_x = (d-1) % USKEWING_BLOCK_SIZE_X + 1;
    const int last_block_y = (d-1) % USKEWING_BLOCK_SIZE_Y + 1;
    int bt;
    for(bt=0;bt + 1 < time_blocks - 1; bt += 2){
        const int t0 = bt * USKEWING_TIMESTEPS + 1;
        const int t1 = (bt+1) * USKEWING_TIMESTEPS + 1;
        const int t2 = (bt+2) * USKEWING_TIMESTEPS + 1;
        uskewed_time_block(b, p, pn, d,
                      t0, t1,
                      x_blocks, last_block_x,
                      y_blocks, last_block_y);
        uskewed_time_block(b, p, pn, d,
                      t1, t2,
                      x_blocks, last_block_x,
                      y_blocks, last_block_y);
    }
    if(bt < time_blocks - 1){
        const int t0 = bt * USKEWING_TIMESTEPS + 1;
        const int t1 = (bt+1) * USKEWING_TIMESTEPS + 1;
        uskewed_time_block(b, p, pn, d,
                      t0, t1,
                      x_blocks, last_block_x,
                      y_blocks, last_block_y);
    }
    uskewed_time_block(b, p, pn, d,
                  pit-last_time+1, pit+1,
                  x_blocks, last_block_x,
                  y_blocks, last_block_y);
    if(pit % 2 == 1) SWAP(double*, sim->p, sim->pn);
    sim->p[0] = sim->p[1];
    sim->p[d-1] = sim->p[d-2];
    sim->p[(d-1)*d+0] = sim->p[(d-1)*d+1];
    sim->p[(d-1)*d+d-1] = sim->p[(d-1)*d+d-2];
}
// 22 * (d-2)(d-2)*pit + 2

/* Advance the simulation sim by one step of size dt using pit iterations
 * for the calculation of pressure. Use the array b for the calculation
 * of the intermediate matrix b.
*/
static void step_uskewed_simulation(
    uskewed_simulation* sim, unsigned int pit, double dt){

    faster_math_build_up_quarter_b(sim, dt);
    uskewed_pressure_poisson(sim, pit);
    const size_t d = sim->d;
    const double ds = sim->size / (d - 1);

    const double rho = sim->rho;
    const double nu = sim->nu;

    const double dtds = dt*(d - 1)/(sim->size);
    //const double dtds = dt*(d - 1)/(2*sim->size);
    const double dtsqds = dt*(d - 1)*(d - 1)/(sim->size*sim->size);
    //const double dtsqds = dt*(d - 1)*(d - 1)/(2*sim->size);
    const double dt2rhods = dt*(d - 1)/(2*sim->size*rho);
    //Swap u and un
    double* tmp = sim->u;
    sim->u = sim->un;
    sim->un = tmp;
    //Swap v and vn
    double* tmp2 = sim->v;
    sim->v = sim->vn;
    sim->vn = tmp2;
    double *restrict u = sim->u;
    double *restrict v = sim->v;
    double *restrict p = sim->p;
    double *restrict un = sim->un;
    double *restrict vn = sim->vn;
    for(size_t i=1;i<d-1;i++){
        for(size_t j=1;j<d-1;j++){
            const double un_here  = un[d*i     + j  ];
            const double un_left  = un[d*i     + j-1];
            const double un_right = un[d*i     + j+1];
            const double un_below = un[d*(i+1) + j  ];
            const double un_above = un[d*(i-1) + j  ];
            const double vn_here  = vn[d*i     + j  ];
            const double vn_left  = vn[d*i     + j-1];
            const double vn_right = vn[d*i     + j+1];
            const double vn_below = vn[d*(i+1) + j  ];
            const double vn_above = vn[d*(i-1) + j  ];
            const double p_left   =  p[d*i     + j-1];
            const double p_right  =  p[d*i     + j+1];
            const double p_below  =  p[d*(i+1) + j  ];
            const double p_above  =  p[d*(i-1) + j  ];
            u[d*i+j] = un_here
                      -dtds * (un_here * (un_here - un_left) +
                               vn_here * (un_here - un_above)) -
                    dt2rhods * (p_right - p_left) +
                    nu * dtsqds *
                (un_right + un_left + un_below + un_above - 4*un_here);


            v[d*i+j] = vn_here
                      -dtds * (un_here * (vn_here - vn_left) +
                               vn_here * (vn_here - vn_above))
                      -dt2rhods * (p_below - p_above)
                      +nu * dtsqds *
                       (vn_right + vn_left + vn_below + vn_above
                        - 4*vn_here);

            }
        }
        for(size_t i=0;i<d;i++){
            u[d*i + 0   ] = 0;
            u[d*i + d-1] = 0;
            v[d*i + 0   ] = 0;
            v[d*i + d-1] = 0;
        }
        for(size_t j=0;j<d;j++){
            u[d*0      + j] = 0;
            u[d*(d-1) + j] = 1; // cavity lid
            v[d*0      + j] = 0;
            v[d*(d-1) + j] = 0;
        }
} // flops ?

/* Advance the simulation sim by steps steps of size dt, using pit
 * iterations for the calculation of pressure.
*/
void advance_uskewed_simulation(uskewed_simulation* sim,
                 unsigned int steps,
                 unsigned int pit, double dt){
    for(unsigned int i=0;i<steps;i++){
        step_uskewed_simulation(sim, pit, dt);
    }
} // Flops ?


void write_uskewed_simulation(uskewed_simulation* sim,
                   FILE* fp){
    write_preallocated_simulation(sim, fp);
}

void destroy_uskewed_simulation(uskewed_simulation* sim){
    deinit_uskewed_simulation(sim);
    free(sim);
}

void deinit_uskewed_simulation(uskewed_simulation *sim){
    deinit_preallocated_simulation(sim);
}
