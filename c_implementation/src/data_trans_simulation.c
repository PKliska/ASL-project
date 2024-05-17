#include <stdlib.h>
#include <stdio.h>
#include "preallocated_simulation.h"
#include "data_trans_simulation.h"
#include <immintrin.h>
#include <stdint.h>
#include <assert.h>


static const struct simulation_vtable_ DATA_TRANS_SIMULATION_VTABLE[] = {{
    .advance=(void (*)(struct simulation *, unsigned int, unsigned int, double))advance_data_trans_simulation,
    .write=(void (*)(const struct simulation *, FILE *))write_data_trans_simulation,
    .destroy=(void (*)(struct simulation *))destroy_data_trans_simulation
}};


data_trans_simulation* new_data_trans_simulation(
    size_t dimension, double size, double rho, double nu){
    data_trans_simulation* sim = new_preallocated_simulation(dimension,
                                                              size,
                                                              rho, nu);
    sim->base.vtable_ = DATA_TRANS_SIMULATION_VTABLE;
    return sim;
}
static double sq(const double x){
    return x*x;
}

static void build_up_b(const data_trans_simulation* sim,
               double dt){
    double *restrict b = sim->b;
    const size_t d = sim->d;
    const double rho = sim->rho;
    const double multiplier = (d - 1) / (2.0*sim->size);
    const double rdt = 1.0 / dt;
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
        b[i*d+j] = rho *
            (rdt *
                ((u_right - u_left + v_below - v_above) * multiplier)
             - (sq(u_right - u_left)
                + 2*((u_below - u_above)*(v_right - v_left))
                + sq(v_below - v_above))*sq(multiplier));

    }
    }
}
 // ? flops

double* transform_arr(double* arr, int d) {
    int chunkSize = d / 4;
    double* res = _mm_malloc(sizeof(double) * d, 32);
    for (int i = 0; i < chunkSize; i++) {
        for (int j = 0; j < 4; j++) {
            res[i * 4 + j] = arr[j * chunkSize + i];
        }
    }
    return res;
}

static void pressure_poisson(data_trans_simulation* sim,
                 unsigned int pit){
    double *restrict b = sim->b;
    const size_t d = sim->d;
    const double sqds = sq(sim->size / (d - 1));

    __m256d pn_vec_above, pn_vec_below, b_vec, pn_trans_vec_left, pn_trans_vec_right, p_vec_l, p_vec_r, p_vec;
    __m256d constant_quarter = _mm256_set1_pd(0.25);
    __m256d constant_msqds = _mm256_set1_pd(-sqds);

    // ? flops
    for(unsigned int q=0;q<pit;q++){
        //Swap p and pn
        double *tmp = sim->p;
        sim->p = sim->pn;
        sim->pn = tmp;

        double *restrict p = sim->p;
        double *restrict pn = sim->pn;
        
        for(size_t i=1;i<d-1;i++){

            for(size_t j=1;j<4;j++){
            const double pn_left =  pn[d*i     + j-1];
            const double pn_right = pn[d*i     + j+1];
            const double pn_below = pn[d*(i+1) + j  ];
            const double pn_above = pn[d*(i-1) + j  ];
            p[i*d+j] = (pn_right + pn_left + pn_below + pn_above 
                        - b[i*d + j]*sqds) * 0.25;
            }
            //double* aligned_pni = _mm_malloc(sizeof(double) * d, 32);
            //memcpy(aligned_pni, pn[d*i], d * sizeof(double));
            double* pn_trans = transform_arr(&pn[d*i], d);
            for(size_t j=4;j<d-7;j+=4){
                //const double pn_left =  pn[d*i     + j-1];
                //const double pn_right = pn[d*i     + j+1];
                //const double pn_below = pn[d*(i+1) + j  ];
                //const double pn_above = pn[d*(i-1) + j  ];
                //p[i*d+j] = ( ( pn_right + pn_left) 
                //              + (pn_below + pn_above - b[i*d + j]*sqds ) ) * 0.25;
                

                pn_trans_vec_left = _mm256_load_pd(pn_trans + j - 4);
                pn_trans_vec_right = _mm256_load_pd(pn_trans + j + 4);
                p_vec_l = _mm256_add_pd(pn_trans_vec_left, pn_trans_vec_right);
           
                pn_vec_above = _mm256_load_pd(&pn[d*(i+1) +j]);
                pn_vec_below = _mm256_load_pd(&pn[d*(i-1)+j]);
                b_vec = _mm256_load_pd(&b[d*i +j]);
               
                p_vec_r = _mm256_fmadd_pd(b_vec, constant_msqds, pn_vec_above);
                p_vec_r = _mm256_add_pd(p_vec_r, pn_vec_below);
                
                p_vec = _mm256_add_pd(p_vec_l, p_vec_r);
                p_vec_r = _mm256_mul_pd(p_vec_r, constant_quarter);

                //uintptr_t address = (uintptr_t)&p[d*i];
                //assert((address % 32) == 0);
     

                _mm256_store_pd(&p[d*i+j], p_vec_r);

                for (int j1 = 0; j1 < 4; j1++) {
                const double pn_left =  pn[d*i     + j-1+j1];
                const double pn_right = pn[d*i     + j+1+j1];
                //const double pn_below = pn[d*(i+1) + j  ];
                //const double pn_above = pn[d*(i-1) + j  ];
                p[i*d+j+j1] += 0.25* ( pn_right + pn_left);
                }
            }

 

            for(size_t j=d-4;j<d-1;j++){
            const double pn_left =  pn[d*i     + j-1];
            const double pn_right = pn[d*i     + j+1];
            const double pn_below = pn[d*(i+1) + j  ];
            const double pn_above = pn[d*(i-1) + j  ];
            p[i*d+j] = (pn_right + pn_left + pn_below + pn_above 
                        - b[i*d + j]*sqds) * 0.25;
            }
        }

        





        for(size_t i=0;i<d;i++) p[d*i     + d-1] = p[d*i + d-2];
        for(size_t i=0;i<d;i++) p[d*i     + 0  ] = p[d*i + 1];
        for(size_t j=0;j<d;j++) p[d*0     + j  ] = p[d*1 + j];
        for(size_t j=0;j<d;j++) p[d*(d-1) + j  ] = 0;
    }
}
// 22 * (d-2)(d-2)*pit + 2

/* Advance the simulation sim by one step of size dt using pit iterations
 * for the calculation of pressure. Use the array b for the calculation
 * of the intermediate matrix b.
*/
static void step_data_trans_simulation(
    data_trans_simulation* sim, unsigned int pit, double dt){
    
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
void advance_data_trans_simulation(data_trans_simulation* sim,
                 unsigned int steps,
                 unsigned int pit, double dt){
    for(unsigned int i=0;i<steps;i++){
        step_data_trans_simulation(sim, pit, dt);
    }
} // Flops ?


void write_data_trans_simulation(data_trans_simulation* sim,
                   FILE* fp){
    write_preallocated_simulation(sim, fp);
}

void destroy_data_trans_simulation(data_trans_simulation* sim){
    destroy_preallocated_simulation(sim);
}
