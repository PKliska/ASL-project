#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "trapeze_simulation.h"
#include "utils.h"
#include "preallocated_simulation.h"
#include "simulation.h"
#include "cutoff_trapeze_simulation.h"

#if defined(TRAPEZE_CUTOFF_VOLUME) && defined(TRAPEZE_CUTOFF_AREA)
#error "At most one of TRAPEZE_CUTOFF_VOLUME, TRAPEZE_CUTOFF_AREA, can be set."
#elif defined(TRAPEZE_CUTOFF_VOLUME) + defined(TRAPEZE_CUTOFF_AREA) + defined(TRAPEZE_CUTOFF_MAX_AREA) == 0
#warning "No TRAPEZE_CUTOFF criterium set, using TRAPEZE_CUTOFF_VOLUME=(1<<18)"
#define TRAPEZE_CUTOFF_VOLUME (1<<18)
#endif

int should_cutoff(const struct trapeze* trap){
    const unsigned t0 = trap->t0, t1 = trap->t1;
    const unsigned x0 = trap->x0, dx0 = trap->dx0, x1 = trap->x1, dx1 = trap->dx1;
    const unsigned y0 = trap->y0, dy0 = trap->dy0, y1 = trap->y1, dy1 = trap->dy1;
    const unsigned dt = t1 - t0;
    const unsigned w = x1 - x0, h = y1 - y0;
    const unsigned dw = dx1 - dx0, dh = dy1 - dy0;
    const unsigned volume = dt * w * h + (h*dw + dh*dw) * (dt - 1)*dt/2 + dh*dw*(dt - 1)*dt*(2*dt - 1)/6;
    const unsigned area = w * h;
#ifdef TRAPEZE_CUTOFF_VOLUME
    return volume <= TRAPEZE_CUTOFF_VOLUME;
#elif defined(TRAPEZE_CUTOFF_AREA)
    return area <= TRAPEZE_CUTOFF_AREA;
#endif
}

static const struct simulation_vtable_ cutoff_trapeze_SIMULATION_VTABLE[] = {{
    .advance=(void (*)(struct simulation *, unsigned int, unsigned int, double))advance_cutoff_trapeze_simulation,
    .write=(void (*)(const struct simulation *, FILE *))write_cutoff_trapeze_simulation,
    .destroy=(void (*)(struct simulation *))destroy_cutoff_trapeze_simulation
}};


cutoff_trapeze_simulation* new_cutoff_trapeze_simulation(
    size_t dimension, double size, double rho, double nu){
    cutoff_trapeze_simulation* sim = malloc(sizeof(cutoff_trapeze_simulation));
    trapeze_simulation* base_sim = &sim->base;
    base_sim->d = dimension;
    base_sim->rho = rho;
    base_sim->nu = nu;
    base_sim->size = size;
    size_t matrix_size = dimension*dimension;
    base_sim->u = zero_array(matrix_size);
    base_sim->un = zero_array(matrix_size);
    base_sim->v = zero_array(matrix_size);
    base_sim->vn = zero_array(matrix_size);
    //add 2 rows of padding
    base_sim->p = zero_array(matrix_size + 2*dimension);
    base_sim->pn = zero_array(matrix_size + 2*dimension);
    base_sim->b = zero_array(matrix_size + 2*dimension);
    base_sim->base.vtable_ = cutoff_trapeze_SIMULATION_VTABLE;
    sim->poisson_len = 0;
    sim->poisson_cap = 64;
    sim->poisson_order = malloc(64 * sizeof(struct trapeze));
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

static inline void  pressure_poisson_do_trapeze(trapeze_simulation* sim,
                                                struct trapeze* trap){
    const unsigned t0 = trap->t0, t1 = trap->t1;
    const int dx0 = trap->dx0, dx1 = trap->dx1;
    const int dy0 = trap->dy0, dy1 = trap->dy1;
    unsigned x0 = trap->x0, x1 = trap->x1;
    unsigned y0 = trap->y0, y1 = trap->y1;
    
    const unsigned d = sim->d;
    double *restrict p_old, *restrict p_new;
    double *restrict b = sim->b;
    if(t0 % 2 == 0){
        p_old = sim->pn;
        p_new = sim->p;
    }else{
        p_new = sim->pn;
        p_old = sim->p;
    }
    for(unsigned t=t0;t<t1;t++){
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
        x0 += dx0;x1 += dx1;
        y0 += dy0;y1 += dy1;
        double *tmp;
        tmp = p_new;
        p_new = p_old;
        p_old = tmp;
    }
}

// Reference
// Cache Oblivious Stencil Computations by Frigo and Strumpen 2005
static void pressure_poisson_walk(cutoff_trapeze_simulation* sim,
                                  struct trapeze* trap){
    const unsigned t0 = trap->t0, t1 = trap->t1;
    const unsigned x0 = trap->x0, dx0 = trap->dx0, x1 = trap->x1, dx1 = trap->dx1;
    const unsigned y0 = trap->y0, dy0 = trap->dy0, y1 = trap->y1, dy1 = trap->dy1;
    const unsigned dt = t1 - t0;
    if(should_cutoff(trap)){
        if(sim->poisson_len == sim->poisson_cap){
            struct trapeze* tmp = malloc(2*sim->poisson_cap*sizeof(struct trapeze));
            memcpy(tmp, sim->poisson_order, sim->poisson_len*sizeof(struct trapeze));
            free(sim->poisson_order);
            sim->poisson_order = tmp;
            sim->poisson_cap *= 2;
        }
        sim->poisson_order[sim->poisson_len] = *trap;
        sim->poisson_len++;
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
static void precompute_poisson(cutoff_trapeze_simulation* sim, unsigned pit){
    sim->poisson_len = 0;
    const unsigned d = sim->base.d;
    struct trapeze trap = {
        .x0 = 1, .dx0 = 0, .x1 = d+1, .dx1 = 0,
        .y0 = 0, .dy0 = 0, .y1 = d, .dy1 = 0,
        .t0 = 1, .t1 = pit+1
    };
    pressure_poisson_walk(sim, &trap);
}

static void pressure_poisson(cutoff_trapeze_simulation* sim,
                 unsigned pit){
    for(size_t i = 0; i<sim->poisson_len;i++){
        pressure_poisson_do_trapeze(&sim->base, &sim->poisson_order[i]);
    }
    if(pit % 2 == 1){
        double* tmp;
        tmp = sim->base.p;
        sim->base.p = sim->base.pn;
        sim->base.pn = tmp;
    }
}
// 22 * (d-2)(d-2)*pit + 2

/* Advance the simulation sim by one step of size dt using pit iterations
 * for the calculation of pressure. Use the array b for the calculation
 * of the intermediate matrix b.
*/
static void step_cutoff_trapeze_simulation(
    cutoff_trapeze_simulation* sim, unsigned int pit, double dt){
    
    trapeze_simulation* base_sim = &sim->base;
    build_up_b(base_sim, dt);
    pressure_poisson(sim, pit);
    const size_t d = base_sim->d;
    const double ds = base_sim->size / (d - 1);

    const double rho = base_sim->rho;
    const double nu = base_sim->nu;

    const double dtds = dt*(d - 1)/(base_sim->size);
    //const double dtds = dt*(d - 1)/(2*base_sim->size);
    const double dtsqds = dt*(d - 1)*(d - 1)/(base_sim->size*base_sim->size);
    //const double dtsqds = dt*(d - 1)*(d - 1)/(2*base_sim->size);
    const double dt2rhods = dt*(d - 1)/(2*base_sim->size*rho);
    //Swap u and un
    double* tmp = base_sim->u;
    base_sim->u = base_sim->un;
    base_sim->un = tmp;
    //Swap v and vn
    double* tmp2 = base_sim->v;
    base_sim->v = base_sim->vn;
    base_sim->vn = tmp2;
    double *restrict u = base_sim->u;
    double *restrict v = base_sim->v;
    //ogmpre [adding row
    double *restrict p = base_sim->p + d;
    double *restrict un = base_sim->un;
    double *restrict vn = base_sim->vn;
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
void advance_cutoff_trapeze_simulation(cutoff_trapeze_simulation* sim,
                 unsigned int steps,
                 unsigned int pit, double dt){
    precompute_poisson(sim, pit);
    for(unsigned int i=0;i<steps;i++){
        step_cutoff_trapeze_simulation(sim, pit, dt);
    }
} // Flops ?


void write_cutoff_trapeze_simulation(cutoff_trapeze_simulation* sim,
                   FILE* fp){
    trapeze_simulation* base_sim = &sim->base;
    //ignore padding row
    base_sim->p += base_sim->d;
    write_preallocated_simulation(base_sim, fp);
    base_sim->p -= base_sim->d;
}

void destroy_cutoff_trapeze_simulation(cutoff_trapeze_simulation* sim){
    trapeze_simulation* base_sim = &sim->base;
    free(base_sim->u);
    free(base_sim->un);
    free(base_sim->v);
    free(base_sim->vn);
    free(base_sim->p);
    free(base_sim->pn);
    free(base_sim->b);
    free(sim->poisson_order);
    free(sim);
}
