#include <stdlib.h>
#include <stdio.h>
#include "utils.h"
#include "preallocated_simulation.h"


// @Pavel casting
static const struct simulation_vtable_ PREALLOCATED_SIMULATION_VTABLE_[] = {{
    .advance=(void (*)(struct simulation *, unsigned int, unsigned int, double))advance_preallocated_simulation,
    .write=(void (*)(const struct simulation *, FILE *))write_preallocated_simulation,
    .destroy=(void (*)(struct simulation *))destroy_preallocated_simulation
}};


struct preallocated_simulation* new_preallocated_simulation(
    size_t dimension, double size, double rho, double nu){
    struct preallocated_simulation *sim = malloc(
    sizeof(struct preallocated_simulation)
    );
    sim->base.vtable_ = PREALLOCATED_SIMULATION_VTABLE_;
    sim->d = dimension;
    sim->rho = rho;
    sim->nu = nu;
    sim->size = size;
    size_t matrix_size = dimension*dimension;
    sim->u = zero_array(matrix_size);
    sim->un = zero_array(matrix_size);
    sim->v = zero_array(matrix_size);
    sim->vn = zero_array(matrix_size);
    sim->p = zero_array(matrix_size);
    sim->pn = zero_array(matrix_size);
    sim->b = zero_array(matrix_size);
    return sim;
}
static double sq(const double x){
    return x*x;
}

static void build_up_b(const struct preallocated_simulation* sim,
               double dt){
    double *restrict b = sim->b;
    const size_t d = sim->d;
    const double rho = sim->rho;
    const double ds = sim->size / (d - 1);
    // 2 flops
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
        b[i*d+j] = rho * (1 / dt *
              ((u_right - u_left) / (2 * ds)
              +(v_below - v_above) / (2 * ds))
              -sq((u_right - u_left) / (2 * ds))
              -2*((u_below - u_above) / (2 * ds)
             *(v_right - v_left) / (2 * ds))
              -sq((v_below - v_above) / (2 * ds)));

    }
    }
} // 29(d-2)(d-2) + 2 flops

static void pressure_poisson(struct preallocated_simulation* sim,
                 unsigned int pit){
    double *restrict b = sim->b;
    const size_t d = sim->d;
    const double ds = sim->size / (d - 1);
    // 2 flops
    for(unsigned int q=0;q<pit;q++){
        //Swap p and pn
        double *tmp = sim->p;
        sim->p = sim->pn;
        sim->pn = tmp;

        double *restrict p = sim->p;
        double *restrict pn = sim->pn;
        for(size_t i=1;i<d-1;i++){
            for(size_t j=1;j<d-1;j++){
            const double pn_left =  pn[d*i     + j-1];
            const double pn_right = pn[d*i     + j+1];
            const double pn_below = pn[d*(i+1) + j  ];
            const double pn_above = pn[d*(i-1) + j  ];
            p[i*d+j] = ((pn_right + pn_left) * sq(ds)
                            +(pn_below + pn_above) * sq(ds)) /
                    (2 * (sq(ds) + sq(ds))) -
                    sq(ds) * sq(ds) / (2 * (sq(ds) + sq(ds))) *
                    b[i*d+j];
            }
        }
        for(size_t i=0;i<d;i++) p[d*i     + d-1] = p[d*i + d-2];
        for(size_t i=0;i<d;i++) p[d*i     + 0  ] = p[d*i + 1];
        for(size_t j=0;j<d;j++) p[d*0     + j  ] = p[d*1 + j];
        for(size_t j=0;j<d;j++) p[d*(d-1) + j  ] = 0;
    }
} // 22 * (d-2)(d-2)*pit + 2

/* Advance the simulation sim by one step of size dt using pit iterations
 * for the calculation of pressure. Use the array b for the calculation
 * of the intermediate matrix b.
*/
static void step_preallocated_simulation(
    struct preallocated_simulation* sim, unsigned int pit, double dt){
    
    build_up_b(sim, dt);
    pressure_poisson(sim, pit);
    const size_t d = sim->d;
    const double ds = sim->size / (d - 1);
    // 2 flops
    const double rho = sim->rho;
    const double nu = sim->nu;
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
            u[d*i+j] = (un_here -
                    un_here * dt / ds *
                    (un_here - un_left) -
                    vn_here * dt / ds *
                    (un_here - un_above) -
                    dt / (2 * rho * ds) * (p_right - p_left) +
                    nu * (dt / sq(ds) *
                (un_right - 2 * un_here + un_left) +
                dt / sq(ds) *
                (un_below - 2 * un_here + un_above)));
            v[d*i+j] = (vn_here -
                un_here * dt / ds *
                (vn_here - vn_left) -
                vn_here * dt / ds *
                (vn_here - vn_above) -
                dt / (2 * rho * ds) * (p_below - p_above) +
                    nu * (dt / sq(ds) *
                (vn_right - 2 * vn_here + vn_left) +
                dt / sq(ds) *
                (vn_below - 2 * vn_here + vn_above)));
            }
        }
        for(size_t i=0;i<d;i++){
            u[d*i + 0  ] = 0;
            u[d*i + d-1] = 0;
            v[d*i + 0  ] = 0;
            v[d*i + d-1] = 0;
        }
        for(size_t j=0;j<d;j++){
            u[d*0     + j] = 0;
            u[d*(d-1) + j] = 1; // cavity lid
            v[d*0     + j] = 0;
            v[d*(d-1) + j] = 0;
        }
} // flops 91 * (d-2)*(d-2) + 22 * (d-2)(d-2)*pit + 4

/* Advance the simulation sim by steps steps of size dt, using pit
 * iterations for the calculation of pressure.
*/
void advance_preallocated_simulation(struct preallocated_simulation* sim,
                 unsigned int steps,
                 unsigned int pit, double dt){
    for(unsigned int i=0;i<steps;i++){
        step_preallocated_simulation(sim, pit, dt);
    }
} // Flops (91 * (d-2)*(d-2) + 22 * (d-2)(d-2)*pit + 4) * steps


void write_preallocated_simulation(struct preallocated_simulation* sim,
                   FILE* fp){
    const size_t COLUMNS = sim->d;
    fprintf(fp, "%zu,%zu,%lf,%lf,", sim->d, sim->d, sim->rho, sim->nu);

    for(size_t i=4;i<COLUMNS;i++) fputc(',', fp);
    fputc('\n', fp);

    write_matrix(sim->u, sim->d, sim->d, fp);

    for(size_t i=0;i<COLUMNS;i++) fputc(',', fp);
    fputc('\n', fp);

    write_matrix(sim->v, sim->d, sim->d, fp);

    for(size_t i=0;i<COLUMNS;i++) fputc(',', fp);
    fputc('\n', fp);

    write_matrix(sim->p, sim->d, sim->d, fp);
}

void destroy_preallocated_simulation(struct preallocated_simulation* sim){
    free(sim->u);
    free(sim->un);
    free(sim->v);
    free(sim->vn);
    free(sim->p);
    free(sim->pn);
    free(sim->b);
    free(sim);
}
