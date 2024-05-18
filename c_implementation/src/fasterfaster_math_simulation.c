#include <stdlib.h>
#include <stdio.h>
#include "preallocated_simulation.h"
#include "fasterfaster_math_simulation.h"


static const struct simulation_vtable_ FASTERFASTER_MATH_SIMULATION_VTABLE[] = {{
    .advance=(void (*)(struct simulation *, unsigned int, unsigned int, double))advance_fasterfaster_math_simulation,
    .write=(void (*)(const struct simulation *, FILE *))write_fasterfaster_math_simulation,
    .destroy=(void (*)(struct simulation *))destroy_fasterfaster_math_simulation
}};

fasterfaster_math_simulation* new_fasterfaster_math_simulation(
    size_t dimension, double size, double rho, double nu){
    fasterfaster_math_simulation* sim = new_preallocated_simulation(dimension,
                                                              size,
                                                              rho, nu);
    sim->base.vtable_ = FASTERFASTER_MATH_SIMULATION_VTABLE;
    return sim;
}

static double sq(const double x){
    return x*x;
}

/* Advance the simulation sim by one step of size dt using pit iterations
 * for the calculation of pressure. Use the array b for the calculation
 * of the intermediate matrix b.
 */
static void step_fasterfaster_math_simulation(
    fasterfaster_math_simulation* sim, unsigned int pit, double dt){
    
    double *restrict b = sim->b;
    const size_t d = sim->d;
    const double sqds = sq(sim->size / (d - 1));
    const double ds = sim->size / (d - 1);
    const double rho = sim->rho;
    const double nu = sim->nu;
    const double dtds = dt*(d - 1)/(sim->size);
    const double dtsqds = dt*(d - 1)*(d - 1)/(sim->size*sim->size);
    const double dt2rhods = dt*(d - 1)/(2*sim->size*rho);
    
    // START BUILD_UP_B

    const double multiplier = (d - 1) / (2.0*sim->size);
    const double sqmultiplier = sq(multiplier);
    const double rdt = 1.0 / dt;

     double *restrict u = sim->u;
     double *restrict v = sim->v;
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
        b[i*d+j] = rho *
            ( rdt *


                ( (u_right - u_left + v_below - v_above) * multiplier )
                - ( (u_right - u_left) * (u_right - u_left)
                     + 2 * ( (u_below - u_above)*(v_right - v_left) )
                      + (v_below - v_above)*(v_below - v_above) )


                * sqmultiplier );


        b[i*d+j] = b[i*d+j] * sqds;
    }
    }
    // END BUILD_UP_B

    // START PRESSURE_POISSON
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

                p[i*d+j] = (pn_right + pn_left + pn_below + pn_above - b[i*d + j]) * 0.25;
            }
        }

        for(size_t i=0;i<d;i++) p[d*i     + d-1] = p[d*i + d-2];
        for(size_t i=0;i<d;i++) p[d*i     + 0  ] = p[d*i + 1];
        for(size_t j=0;j<d;j++) p[d*0     + j  ] = p[d*1 + j];
        for(size_t j=0;j<d;j++) p[d*(d-1) + j  ] = 0;
    }
    // END PRESSURE_POISSON

    //Swap u and un
    double* tmp = sim->u;
    sim->u = sim->un;
    sim->un = tmp;
    //Swap v and vn
    double* tmp2 = sim->v;
    sim->v = sim->vn;
    sim->vn = tmp2;

u = sim->u;
v = sim->v;
 
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
void advance_fasterfaster_math_simulation(fasterfaster_math_simulation* sim,
                 unsigned int steps,
                 unsigned int pit, double dt){
    for(unsigned int i=0;i<steps;i++){
            double *restrict b = sim->b;
    const size_t d = sim->d;
    const double sqds = sq(sim->size / (d - 1));
    const double ds = sim->size / (d - 1);
    const double rho = sim->rho;
    const double nu = sim->nu;
    const double dtds = dt*(d - 1)/(sim->size);
    const double dtsqds = dt*(d - 1)*(d - 1)/(sim->size*sim->size);
    const double dt2rhods = dt*(d - 1)/(2*sim->size*rho);
    
    // START BUILD_UP_B

    const double multiplier = (d - 1) / (2.0*sim->size);
    const double sqmultiplier = sq(multiplier);
    const double rdt = 1.0 / dt;

     double *restrict u = sim->u;
     double *restrict v = sim->v;
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
        b[i*d+j] = rho *
            ( rdt *


                ( (u_right - u_left + v_below - v_above) * multiplier )
                - ( (u_right - u_left) * (u_right - u_left)
                     + 2 * ( (u_below - u_above)*(v_right - v_left) )
                      + (v_below - v_above)*(v_below - v_above) )

                      
                * sqmultiplier );


        b[i*d+j] = b[i*d+j] * sqds;
    }
    }
    // END BUILD_UP_B

    // START PRESSURE_POISSON
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

                p[i*d+j] = (pn_right + pn_left + pn_below + pn_above - b[i*d + j]) * 0.25;
            }
        }

        for(size_t i=0;i<d;i++) p[d*i     + d-1] = p[d*i + d-2];
        for(size_t i=0;i<d;i++) p[d*i     + 0  ] = p[d*i + 1];
        for(size_t j=0;j<d;j++) p[d*0     + j  ] = p[d*1 + j];
        for(size_t j=0;j<d;j++) p[d*(d-1) + j  ] = 0;
    }
    // END PRESSURE_POISSON

    //Swap u and un
    double* tmp = sim->u;
    sim->u = sim->un;
    sim->un = tmp;
    //Swap v and vn
    double* tmp2 = sim->v;
    sim->v = sim->vn;
    sim->vn = tmp2;

u = sim->u;
v = sim->v;
 
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
    }
} // Flops ?


void write_fasterfaster_math_simulation(fasterfaster_math_simulation* sim,
                   FILE* fp){
    write_preallocated_simulation(sim, fp);
}

void destroy_fasterfaster_math_simulation(fasterfaster_math_simulation* sim){
    destroy_preallocated_simulation(sim);
}
