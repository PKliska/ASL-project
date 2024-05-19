#include <stdio.h>
#include <stdlib.h>
#include "utils.h"
#include "preallocated_simulation.h"
#include "simulation.h"
#include "trapeze_simulation.h"


static const struct simulation_vtable_ TRAPEZE_SIMULATION_VTABLE[] = {{
    .advance=(void (*)(struct simulation *, unsigned int, unsigned int, double))advance_trapeze_simulation,
    .write=(void (*)(const struct simulation *, FILE *))write_trapeze_simulation,
    .destroy=(void (*)(struct simulation *))destroy_trapeze_simulation
}};


trapeze_simulation* new_trapeze_simulation(
    size_t dimension, double size, double rho, double nu){
    trapeze_simulation* sim = malloc(sizeof(trapeze_simulation));
    sim->d = dimension;
    sim->rho = rho;
    sim->nu = nu;
    sim->size = size;
    size_t matrix_size = dimension*dimension;
    sim->u = zero_array(matrix_size);
    sim->un = zero_array(matrix_size);
    sim->v = zero_array(matrix_size);
    sim->vn = zero_array(matrix_size);
    //add 2 rows of padding
    sim->p = zero_array(matrix_size + 2*dimension);
    sim->pn = zero_array(matrix_size + 2*dimension);
    sim->b = zero_array(matrix_size + 2*dimension);
    sim->base.vtable_ = TRAPEZE_SIMULATION_VTABLE;
    return sim;
}
static double sq(const double x){
    return x*x;
}

static void build_up_b(const trapeze_simulation* sim,
               double dt){
    //ignore padding row
    double *restrict b = sim->b + sim->d;
    const size_t d = sim->d;
    const double rho = sim->rho;
    const double multiplier = (d - 1) / (2.0*sim->size);
    const double rdt = 1.0 / dt;
    const double sqds = sq(sim->size / (d - 1));
    // ? flops
    const double *restrict u = sim->u;
    const double *restrict v = sim->v;
    for(size_t i=1;i<d - 1;i++){
    for(size_t j=1;j<d - 1;j++){
        const double u_left =  u[d*i     + j-1];
        const double u_right = u[d*i     + j+1];
        const double u_below = u[d*(i+1) + j  ];
        const double u_above = u[d*(i-1) + j  ];
        const double v_left =  v[d*i     + j-1];
        const double v_right = v[d*i     + j+1];
        const double v_below = v[d*(i+1) + j  ];
        const double v_above = v[d*(i-1) + j  ];
        b[i*d+j] = sqds * rho *
            (rdt *
                ((u_right - u_left + v_below - v_above) * multiplier)
             - (sq(u_right - u_left)
                + 2*((u_below - u_above)*(v_right - v_left))
                + sq(v_below - v_above))*sq(multiplier));

    }
    }
}

static inline void  pressure_poisson_base_case(trapeze_simulation* sim,
                                              unsigned t,
                                              unsigned x0, unsigned x1,
                                              unsigned y0, unsigned y1){
    const unsigned d = sim->d;
    double *restrict p_old, *restrict p_new;
    double *restrict b = sim->b;
    if(t % 2 == 0){
        p_old = sim->pn;
        p_new = sim->p;
    }else{
        p_new = sim->pn;
        p_old = sim->p;
    }
    for(unsigned x=x0;x<x1;x++){
    for(unsigned y=y0;y<y1;y++){
        p_new[x*d+y] = (p_old[(x+1)*d + y] + p_old[(x-1)*d + y]
                      + p_old[x*d + (y+1)] + p_old[x*d + (y-1)]
                      - b[x*d+y]) * 0.25;
    }
    }
    if(y1 == d) for(unsigned x=x0;x<x1;x++) p_new[x*d + (d-1)] = p_new[x*d +(d-2)];
    if(y0 == 0) for(unsigned x=x0;x<x1;x++) p_new[x*d + 0] = p_new[x*d + 1];
    if(x0 == 1) for(unsigned y=y0;y<y1;y++) p_new[1*d + y] = p_new[2*d + y];
    if(x1 == d+1) for(unsigned y=y0;y<y1;y++) p_new[d*d + y] = 0;

}

// Reference
// Cache Oblivious Stencil Computations by Frigo and Strumpen 2005
static void pressure_poisson_walk(trapeze_simulation* sim,
                                  struct trapeze* trap){
    const unsigned t0 = trap->t0, t1 = trap->t1;
    const unsigned x0 = trap->x0, dx0 = trap->dx0, x1 = trap->x1, dx1 = trap->dx1;
    const unsigned y0 = trap->y0, dy0 = trap->dy0, y1 = trap->y1, dy1 = trap->dy1;
    const unsigned dt = t1 - t0;
    if(dt == 1){
        pressure_poisson_base_case(sim, t0, x0, x1, y0, y1);
        return;
    }
    //Try cut space along x
    if(2 * (x1 - x0) + (dx1 - dx0) * dt >= 4 * dt){
        unsigned xm = (2 * (x0+x1) + (2 + dx0 + dx1) * dt) / 4;
        
        trap->x1 = xm; trap->dx1 = -1;
        pressure_poisson_walk(sim, trap);
        trap->x1 = x1; trap->dx1 = dx1;
        
        trap->x0 = xm; trap->dx0 = -1;
        pressure_poisson_walk(sim, trap);
        trap->x0 = x0; trap->dx0 = dx0;
        return;
    }
    if(2 * (y1 - y0) + (dy1 - dy0) * dt >= 4 * dt){
        unsigned ym = (2 * (y0+y1) + (2 + dy0 + dy1) * dt) / 4;
        
        trap->y1 = ym; trap->dy1 = -1;
        pressure_poisson_walk(sim, trap);
        trap->y1 = y1; trap->dy1 = dy1;
        
        trap->y0 = ym; trap->dy0 = -1;
        pressure_poisson_walk(sim, trap);
        trap->y0 = y0; trap->dy0 = dy0;
        return;
    }
    unsigned s = dt / 2;
    
    trap->t1 = t0 + s;
    pressure_poisson_walk(sim, trap);
    trap->t1 = t1;

    trap->t0 = t0 + s;
    trap->x0 = x0 + dx0*s; trap->x1 = x1 + dx1*s;
    trap->y0 = y0 + dy0*s; trap->y1 = y1 + dy1*s;
    pressure_poisson_walk(sim, trap);
    trap->x0 = x0; trap->x1 = x1;
    trap->y0 = y0; trap->y1 = y1;
    trap->t0 = t0;
}

static void pressure_poisson(trapeze_simulation* sim,
                 unsigned pit){
    const unsigned d = sim->d;
    struct trapeze trap = {
        .x0 = 1, .dx0 = 0, .x1 = d+1, .dx1 = 0,
        .y0 = 0, .dy0 = 0, .y1 = d, .dy1 = 0,
        .t0 = 1, .t1 = pit+1
    };
    pressure_poisson_walk(sim, &trap);
    if(pit % 2 == 1){
        double* tmp;
        tmp = sim->p;
        sim->p = sim->pn;
        sim->pn = tmp;
    }
}
// 22 * (d-2)(d-2)*pit + 2

/* Advance the simulation sim by one step of size dt using pit iterations
 * for the calculation of pressure. Use the array b for the calculation
 * of the intermediate matrix b.
*/
static void step_trapeze_simulation(
    trapeze_simulation* sim, unsigned int pit, double dt){
    
    build_up_b(sim, dt);
    pressure_poisson(sim, pit);
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
    //ogmpre [adding row
    double *restrict p = sim->p + d;
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
void advance_trapeze_simulation(trapeze_simulation* sim,
                 unsigned int steps,
                 unsigned int pit, double dt){
    for(unsigned int i=0;i<steps;i++){
        step_trapeze_simulation(sim, pit, dt);
    }
} // Flops ?


void write_trapeze_simulation(trapeze_simulation* sim,
                   FILE* fp){
    //ignore padding row
    sim->p += sim->d;
    write_preallocated_simulation(sim, fp);
    sim->p -= sim->d;
}

void destroy_trapeze_simulation(trapeze_simulation* sim){
    destroy_preallocated_simulation(sim);
}
