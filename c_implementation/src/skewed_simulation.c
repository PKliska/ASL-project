#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include "skewed_simulation.h"
#include "faster_math_simulation.h"
#include "utils.h"
#include "trapeze_macros.h"

#ifndef SKEWING_BLOCK_SIZE
#warning "SKEWING_BLOCK_SIZE was not defined, setting to default (16)"
#define SKEWING_BLOCK_SIZE (36)
#endif
#ifndef SKEWING_TIMESTEPS
#warning "SKEWING_TIMESTEPS was not defined, setting to default (10)"
#define SKEWING_TIMESTEPS (10)
#endif
#if SKEWING_BLOCK_SIZE <= SKEWING_TIMESTEPS
#error "SKEWING_BLOCK_SIZE must be greater than SKEWING_TIMESTEPS"
#endif



static const struct simulation_vtable_ SKEWED_SIMULATION_VTABLE[] = {{
    .advance=(void (*)(struct simulation *, unsigned int, unsigned int, double))advance_skewed_simulation,
    .write=(void (*)(const struct simulation *, FILE *))write_skewed_simulation,
    .destroy=(void (*)(struct simulation *))destroy_skewed_simulation
}};


skewed_simulation* new_skewed_simulation(
    size_t dimension, double size, double rho, double nu){
    skewed_simulation* sim = malloc(sizeof(*sim)); 
    init_skewed_simulation(sim, dimension, size, rho, nu);
    return sim;
}

void init_skewed_simulation(skewed_simulation* sim,
    size_t dimension, double size, double rho, double nu){
    init_faster_math_simulation(sim, dimension, size, rho, nu);
    assert(dimension % SKEWING_BLOCK_SIZE == 0);
    assert(dimension / SKEWING_BLOCK_SIZE >= 3);
    sim->base.vtable_ = SKEWED_SIMULATION_VTABLE;
}


static double sq(const double x){
    return x*x;
}

static void skewed_pressure_poisson(skewed_simulation* sim,
                 unsigned int pit){
    const size_t d = sim->d;
    double *restrict b = sim->b;
    double *restrict p = sim->p;
    double *restrict pn = sim->pn;
    for(int bt=0;bt < 4; bt+=2){
        int t0 = bt * SKEWING_TIMESTEPS + 1;
        int t1 = (bt + 1) * SKEWING_TIMESTEPS + 1;
        int t2 = (bt + 1) * SKEWING_TIMESTEPS + 1;
        int t3 = (bt + 2) * SKEWING_TIMESTEPS + 1;

        //handle first block of first row
        DO_TRAPEZE_TOP_LEFT(b, p, pn, d,
                   t0, t1,
                   0, 0, SKEWING_BLOCK_SIZE, -1,
                   0, 0, SKEWING_BLOCK_SIZE, -1);
        for(int bj=1;bj<d/SKEWING_BLOCK_SIZE-1;bj++){
            int y0 = bj*SKEWING_BLOCK_SIZE;
            int y1 = (bj+1)*SKEWING_BLOCK_SIZE;
            DO_TRAPEZE_TOP(b, p, pn, d,
                       t0, t1,
                       0, 0, SKEWING_BLOCK_SIZE, -1,
                       y0, -1, y1, -1);
        }
        DO_TRAPEZE_TOP_RIGHT(b, p, pn, d,
                   t0, t1,
                   0, 0, SKEWING_BLOCK_SIZE, -1,
                   d-SKEWING_BLOCK_SIZE, -1, d, 0);
        for(int bi=1;bi<d/SKEWING_BLOCK_SIZE-1;bi++){
            int x0 = bi*SKEWING_BLOCK_SIZE;
            int x1 = (bi+1)*SKEWING_BLOCK_SIZE;
            DO_TRAPEZE_LEFT(b, p, pn, d,
                       t0, t1,
                       x0, -1, x1, -1,
                       0, 0, SKEWING_BLOCK_SIZE, -1);
            for(int bj=1;bj<d/SKEWING_BLOCK_SIZE-1;bj++){
                int y0 = bj*SKEWING_BLOCK_SIZE;
                int y1 = (bj+1)*SKEWING_BLOCK_SIZE;
                DO_TRAPEZE_MID(b, p, pn, d,
                           t0, t1,
                           x0, -1, x1, -1,
                           y0, -1, y1, -1);
            }
            DO_TRAPEZE_RIGHT(b, p, pn, d,
                       t0, t1,
                       x0, -1, x1, -1,
                       d-SKEWING_BLOCK_SIZE, -1, d, 0);
        }
        DO_TRAPEZE_BOTTOM_LEFT(b, p, pn, d,
                   t0, t1,
                   d-SKEWING_BLOCK_SIZE, -1, d, 0,
                   0, 0, SKEWING_BLOCK_SIZE, -1);
        for(int bj=1;bj<d/SKEWING_BLOCK_SIZE-1;bj++){
            int y0 = bj*SKEWING_BLOCK_SIZE;
            int y1 = (bj+1)*SKEWING_BLOCK_SIZE;
            DO_TRAPEZE_BOTTOM(b, p, pn, d,
                       t0, t1,
                       d-SKEWING_BLOCK_SIZE, -1, d, 0,
                       y0, -1, y1, -1);
        }
        DO_TRAPEZE_BOTTOM_RIGHT(b, p, pn, d,
                   t0, t1,
                   d-SKEWING_BLOCK_SIZE, -1, d, 0,
                   d-SKEWING_BLOCK_SIZE, -1, d, 0);

        //unroll once
        DO_TRAPEZE_TOP_LEFT(b, p, pn, d,
                   t2, t3,
                   0, 0, SKEWING_BLOCK_SIZE, -1,
                   0, 0, SKEWING_BLOCK_SIZE, -1);
        for(int bj=1;bj<d/SKEWING_BLOCK_SIZE-1;bj++){
            int y0 = bj*SKEWING_BLOCK_SIZE;
            int y1 = (bj+1)*SKEWING_BLOCK_SIZE;
            DO_TRAPEZE_TOP(b, p, pn, d,
                       t2, t3,
                       0, 0, SKEWING_BLOCK_SIZE, -1,
                       y0, -1, y1, -1);
        }
        DO_TRAPEZE_TOP_RIGHT(b, p, pn, d,
                   t2, t3,
                   0, 0, SKEWING_BLOCK_SIZE, -1,
                   d-SKEWING_BLOCK_SIZE, -1, d, 0);
        for(int bi=1;bi<d/SKEWING_BLOCK_SIZE-1;bi++){
            int x0 = bi*SKEWING_BLOCK_SIZE;
            int x1 = (bi+1)*SKEWING_BLOCK_SIZE;
            DO_TRAPEZE_LEFT(b, p, pn, d,
                       t2, t3,
                       x0, -1, x1, -1,
                       0, 0, SKEWING_BLOCK_SIZE, -1);
            for(int bj=1;bj<d/SKEWING_BLOCK_SIZE-1;bj++){
                int y0 = bj*SKEWING_BLOCK_SIZE;
                int y1 = (bj+1)*SKEWING_BLOCK_SIZE;
                DO_TRAPEZE_MID(b, p, pn, d,
                           t2, t3,
                           x0, -1, x1, -1,
                           y0, -1, y1, -1);
            }
            DO_TRAPEZE_RIGHT(b, p, pn, d,
                       t2, t3,
                       x0, -1, x1, -1,
                       d-SKEWING_BLOCK_SIZE, -1, d, 0);
        }
        DO_TRAPEZE_BOTTOM_LEFT(b, p, pn, d,
                   t2, t3,
                   d-SKEWING_BLOCK_SIZE, -1, d, 0,
                   0, 0, SKEWING_BLOCK_SIZE, -1);
        for(int bj=1;bj<d/SKEWING_BLOCK_SIZE-1;bj++){
            int y0 = bj*SKEWING_BLOCK_SIZE;
            int y1 = (bj+1)*SKEWING_BLOCK_SIZE;
            DO_TRAPEZE_BOTTOM(b, p, pn, d,
                       t2, t3,
                       d-SKEWING_BLOCK_SIZE, -1, d, 0,
                       y0, -1, y1, -1);
        }
        DO_TRAPEZE_BOTTOM_RIGHT(b, p, pn, d,
                   t2, t3,
                   d-SKEWING_BLOCK_SIZE, -1, d, 0,
                   d-SKEWING_BLOCK_SIZE, -1, d, 0);
    }
    int t9 = 4 * SKEWING_TIMESTEPS + 1;
    int t10 = 5 * SKEWING_TIMESTEPS + 1;
    DO_TRAPEZE_TOP_LEFT(b, p, pn, d,
                   t9, t10,
                   0, 0, SKEWING_BLOCK_SIZE, -1,
                   0, 0, SKEWING_BLOCK_SIZE, -1);
        for(int bj=1;bj<d/SKEWING_BLOCK_SIZE-1;bj++){
            int y0 = bj*SKEWING_BLOCK_SIZE;
            int y1 = (bj+1)*SKEWING_BLOCK_SIZE;
            DO_TRAPEZE_TOP(b, p, pn, d,
                       t9, t10,
                       0, 0, SKEWING_BLOCK_SIZE, -1,
                       y0, -1, y1, -1);
        }
        DO_TRAPEZE_TOP_RIGHT(b, p, pn, d,
                   t9, t10,
                   0, 0, SKEWING_BLOCK_SIZE, -1,
                   d-SKEWING_BLOCK_SIZE, -1, d, 0);
        for(int bi=1;bi<d/SKEWING_BLOCK_SIZE-1;bi++){
            int x0 = bi*SKEWING_BLOCK_SIZE;
            int x1 = (bi+1)*SKEWING_BLOCK_SIZE;
            DO_TRAPEZE_LEFT(b, p, pn, d,
                       t9, t10,
                       x0, -1, x1, -1,
                       0, 0, SKEWING_BLOCK_SIZE, -1);
            for(int bj=1;bj<d/SKEWING_BLOCK_SIZE-1;bj++){
                int y0 = bj*SKEWING_BLOCK_SIZE;
                int y1 = (bj+1)*SKEWING_BLOCK_SIZE;
                DO_TRAPEZE_MID(b, p, pn, d,
                           t9, t10,
                           x0, -1, x1, -1,
                           y0, -1, y1, -1);
            }
            DO_TRAPEZE_RIGHT(b, p, pn, d,
                       t9, t10,
                       x0, -1, x1, -1,
                       d-SKEWING_BLOCK_SIZE, -1, d, 0);
        }
        DO_TRAPEZE_BOTTOM_LEFT(b, p, pn, d,
                   t9, t10,
                   d-SKEWING_BLOCK_SIZE, -1, d, 0,
                   0, 0, SKEWING_BLOCK_SIZE, -1);
        for(int bj=1;bj<d/SKEWING_BLOCK_SIZE-1;bj++){
            int y0 = bj*SKEWING_BLOCK_SIZE;
            int y1 = (bj+1)*SKEWING_BLOCK_SIZE;
            DO_TRAPEZE_BOTTOM(b, p, pn, d,
                       t9, t10,
                       d-SKEWING_BLOCK_SIZE, -1, d, 0,
                       y0, -1, y1, -1);
        }
        DO_TRAPEZE_BOTTOM_RIGHT(b, p, pn, d,
                   t9, t10,
                   d-SKEWING_BLOCK_SIZE, -1, d, 0,
                   d-SKEWING_BLOCK_SIZE, -1, d, 0);
    //if(pit % 2 == 1) SWAP(double*, sim->p, sim->pn);
}
// 22 * (d-2)(d-2)*pit + 2

/* Advance the simulation sim by one step of size dt using pit iterations
 * for the calculation of pressure. Use the array b for the calculation
 * of the intermediate matrix b.
*/
static void step_skewed_simulation(
    skewed_simulation* sim, unsigned int pit, double dt){
    
    faster_math_build_up_b(sim, dt);
    skewed_pressure_poisson(sim, pit);
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
void advance_skewed_simulation(skewed_simulation* sim,
                 unsigned int steps,
                 unsigned int pit, double dt){
    assert(pit % SKEWING_TIMESTEPS == 0);
    for(unsigned int i=0;i<steps;i++){
        step_skewed_simulation(sim, pit, dt);
    }
} // Flops ?


void write_skewed_simulation(skewed_simulation* sim,
                   FILE* fp){
    write_preallocated_simulation(sim, fp);
}

void destroy_skewed_simulation(skewed_simulation* sim){
    deinit_skewed_simulation(sim);
    free(sim);
}

void deinit_skewed_simulation(skewed_simulation *sim){
    deinit_preallocated_simulation(sim);
}
