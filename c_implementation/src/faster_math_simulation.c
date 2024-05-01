#include <stdlib.h>
#include <stdio.h>
#include "preallocated_simulation.h"
#include "faster_math_simulation.h"


static const struct simulation_vtable_ FASTER_MATH_SIMULATION_VTABLE[] = {{
    .advance=(void (*)(struct simulation *, unsigned int, unsigned int, double))advance_faster_math_simulation,
    .write=(void (*)(const struct simulation *, FILE *))write_faster_math_simulation,
    .destroy=(void (*)(struct simulation *))destroy_faster_math_simulation
}};


faster_math_simulation* new_faster_math_simulation(
    size_t nx, size_t ny, double rho, double nu){
    faster_math_simulation* sim = new_preallocated_simulation(nx, ny, rho, nu);
    sim->base.vtable_ = FASTER_MATH_SIMULATION_VTABLE;
    return sim;
}
static double sq(const double x){
    return x*x;
}

static void build_up_b(const faster_math_simulation* sim,
               double dt){
    double *restrict b = sim->b;
    const size_t nx = sim->nx;
    const size_t ny = sim->ny;
    const double rho = sim->rho;
    const double x_mult = (nx - 1) / 4.0;
    const double y_mult = (ny - 1) / 4.0;
    const double rdt = 1.0 / dt;
    // ? flops
    const double *restrict u = sim->u;
    const double *restrict v = sim->v;
    for(size_t i=1;i<nx - 1;i++){
    for(size_t j=1;j<ny - 1;j++){
        const double u_left =  u[ny*i     + j-1];
        const double u_right = u[ny*i     + j+1];
        const double u_below = u[ny*(i+1) + j  ];
        const double u_above = u[ny*(i-1) + j  ];
        const double v_left =  v[ny*i     + j-1];
        const double v_right = v[ny*i     + j+1];
        const double v_below = v[ny*(i+1) + j  ];
        const double v_above = v[ny*(i-1) + j  ];
        b[i*ny+j] = rho * (rdt *
              ((u_right - u_left) * x_mult
              +(v_below - v_above) * y_mult)
              -sq((u_right - u_left) * x_mult)
              -2*((u_below - u_above) * y_mult
             *(v_right - v_left) * x_mult)
              -sq((v_below - v_above) * y_mult));

    }
    }
} // ? flops

static void pressure_poisson(faster_math_simulation* sim,
                 unsigned int pit){
    double *restrict b = sim->b;
    const size_t nx = sim->nx;
    const size_t ny = sim->ny;
    const double sqdx = sq(2.0 / (nx - 1));
    const double sqdy = sq(2.0 / (ny - 1));
    const double sqdxsqdy = sqdx * sqdy;
    const double r2sqdxsqdy = ((long long) (nx-1)*(nx-1)*(ny-1)*(ny-1))
                                           /(8.0*
                                          ((nx-1)*(nx-1)+(ny-1)*(ny-1)));
    // ? flops
    for(unsigned int q=0;q<pit;q++){
        //Swap p and pn
        double *tmp = sim->p;
        sim->p = sim->pn;
        sim->pn = tmp;

        double *restrict p = sim->p;
        double *restrict pn = sim->pn;
        for(size_t i=1;i<nx-1;i++){
            for(size_t j=1;j<ny-1;j++){
            const double pn_left =  pn[ny*i     + j-1];
            const double pn_right = pn[ny*i     + j+1];
            const double pn_below = pn[ny*(i+1) + j  ];
            const double pn_above = pn[ny*(i-1) + j  ];
            p[i*ny+j] = ((pn_right + pn_left) * sqdy
                        +(pn_below + pn_above) * sqdx
                        - sqdxsqdy * b[i*ny + j])*r2sqdxsqdy;
            }
        }
        for(size_t i=0;i<nx;i++) p[ny*i      + ny-1] = p[ny*i + ny-2];
        for(size_t i=0;i<nx;i++) p[ny*i      + 0   ] = p[ny*i + 1];
        for(size_t j=0;j<ny;j++) p[ny*0      + j   ] = p[ny*1 + j];
        for(size_t j=0;j<ny;j++) p[ny*(nx-1) + j   ] = 0;
    }
} // 22 * (nx-2)(ny-2)*pit + 2

/* Advance the simulation sim by one step of size dt using pit iterations
 * for the calculation of pressure. Use the array b for the calculation
 * of the intermediate matrix b.
*/
static void step_faster_math_simulation(
    faster_math_simulation* sim, unsigned int pit, double dt){
    
    const size_t nx = sim->nx;
    const size_t ny = sim->ny;
    const double rho = sim->rho;
    const double nu = sim->nu;
    const double dtdx = dt*(nx - 1)/2.0;
    const double dtdy = dt*(ny - 1)/2.0;
    const double dtsqdx = dt*(nx - 1)*(nx - 1)/4.0;
    const double dtsqdy = dt*(ny - 1)*(ny - 1)/4.0;
    const double dt2rhodx = dt*(nx - 1)/(4.0*rho);
    const double dt2rhody = dt*(ny - 1)/(4.0*rho);
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
    build_up_b(sim, dt);
    pressure_poisson(sim, pit);
    for(size_t i=1;i<nx-1;i++){
        for(size_t j=1;j<ny-1;j++){
            const double un_here  = un[ny*i     + j  ];
            const double un_left  = un[ny*i     + j-1];
            const double un_right = un[ny*i     + j+1];
            const double un_below = un[ny*(i+1) + j  ];
            const double un_above = un[ny*(i-1) + j  ];
            const double vn_here  = vn[ny*i     + j  ];
            const double vn_left  = vn[ny*i     + j-1];
            const double vn_right = vn[ny*i     + j+1];
            const double vn_below = vn[ny*(i+1) + j  ];
            const double vn_above = vn[ny*(i-1) + j  ];
            const double p_left   =  p[ny*i     + j-1];
            const double p_right  =  p[ny*i     + j+1];
            const double p_below  =  p[ny*(i+1) + j  ];
            const double p_above  =  p[ny*(i-1) + j  ];
            u[ny*i+j] = (un_here -
                    un_here * dtdx *
                    (un_here - un_left) -
                    vn_here * dtdy *
                    (un_here - un_above) -
                    dt2rhodx * (p_right - p_left) +
                    nu * (dtsqdx *
                (un_right - 2 * un_here + un_left) +
                dtsqdy *
                (un_below - 2 * un_here + un_above)));
            v[ny*i+j] = (vn_here -
                un_here * dtdx *
                (vn_here - vn_left) -
                vn_here * dtdy *
                (vn_here - vn_above) -
                dt2rhody * (p_below - p_above) +
                    nu * (dtsqdx *
                (vn_right - 2 * vn_here + vn_left) +
                dtsqdy *
                (vn_below - 2 * vn_here + vn_above)));
            }
        }
        for(size_t i=0;i<nx;i++){
            u[ny*i + 0   ] = 0;
            u[ny*i + ny-1] = 0;
            v[ny*i + 0   ] = 0;
            v[ny*i + ny-1] = 0;
        }
        for(size_t j=0;j<ny;j++){
            u[ny*0      + j] = 0;
            u[ny*(nx-1) + j] = 1; // cavity lid
            v[ny*0      + j] = 0;
            v[ny*(nx-1) + j] = 0;
        }
} // flops ? 

/* Advance the simulation sim by steps steps of size dt, using pit
 * iterations for the calculation of pressure.
*/
void advance_faster_math_simulation(faster_math_simulation* sim,
                 unsigned int steps,
                 unsigned int pit, double dt){
    for(unsigned int i=0;i<steps;i++){
        step_faster_math_simulation(sim, pit, dt);
    }
} // Flops ?


void write_faster_math_simulation(faster_math_simulation* sim,
                   FILE* fp){
    write_preallocated_simulation(sim, fp);
}

void destroy_faster_math_simulation(faster_math_simulation* sim){
    destroy_preallocated_simulation(sim);
}
