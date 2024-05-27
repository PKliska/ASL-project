#include <stdlib.h>
#include <stdio.h>
#include "preallocated_simulation.h"
#include "vanilla_avx_simulation.h"
#include <immintrin.h>
#include <stdint.h>
#include <assert.h>


static const struct simulation_vtable_ VANILLA_AVX_SIMULATION_VTABLE[] = {{
    .advance=(void (*)(struct simulation *, unsigned int, unsigned int, double))advance_vanilla_avx_simulation,
    .write=(void (*)(const struct simulation *, FILE *))write_vanilla_avx_simulation,
    .destroy=(void (*)(struct simulation *))destroy_vanilla_avx_simulation
}};


vanilla_avx_simulation* new_vanilla_avx_simulation(
    size_t dimension, double size, double rho, double nu){
    vanilla_avx_simulation* sim = new_preallocated_simulation(dimension,
                                                              size,
                                                              rho, nu);
    sim->base.vtable_ = VANILLA_AVX_SIMULATION_VTABLE;
    return sim;
}
static double sq(const double x){
    return x*x;
}

static void build_up_b(const vanilla_avx_simulation* sim,
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

static void pressure_poisson(vanilla_avx_simulation* sim,
                 unsigned int pit){
    double *restrict b = sim->b;
    const size_t d = sim->d;
    const double sqds = sq(sim->size / (d - 1));

    __m256d pn_vec_above, pn_vec_below, b_vec, pn_trans_vec_left, pn_trans_vec_right, p_vec_l, p_vec_r, p_vec, pn_vec_left, pn_vec_right;
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



            for(size_t j=1;j<d-4;j+=4){
                //const double pn_left =  pn[d*i     + j-1];
                //const double pn_right = pn[d*i     + j+1];
                //const double pn_below = pn[d*(i+1) + j  ];
                //const double pn_above = pn[d*(i-1) + j  ];
                //p[i*d+j] = ( ( pn_right + pn_left) 
                //              + (pn_below + pn_above - b[i*d + j]*sqds ) ) * 0.25;
                
                pn_vec_left = _mm256_loadu_pd(&pn[d*i + j - 1]);
                pn_vec_right = _mm256_loadu_pd(&pn[d*i + j + 1]);

                pn_vec_above = _mm256_loadu_pd(&pn[d*(i+1) +j]);
                pn_vec_below = _mm256_loadu_pd(&pn[d*(i-1)+j]);
                b_vec = _mm256_loadu_pd(&b[d*i +j]);
               
                p_vec = _mm256_fmadd_pd(b_vec, constant_msqds, pn_vec_above);
                p_vec = _mm256_add_pd(p_vec, pn_vec_below);
                
                p_vec = _mm256_add_pd(p_vec, pn_vec_left);
                p_vec = _mm256_add_pd(p_vec, pn_vec_right);
                p_vec = _mm256_mul_pd(p_vec, constant_quarter);
     

                _mm256_storeu_pd(&p[d*i+j], p_vec);

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
static void step_vanilla_avx_simulation(
    vanilla_avx_simulation* sim, unsigned int pit, double dt){
    
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
   __m256d un_here, un_top, un_bot, un_l, un_r, vn_here, vn_top, vn_bot, vn_l, vn_r, p_l, p_r, p_bot, p_top; 
    __m256d un_here_minus_left, un_here_minus_above, vn_here_mult_un_here_minus_above, un_here_fma_vn_here_mult, un_here_fma_mdtds_un_here;
    __m256d p_r_minus_p_l, dt2rhods_mult_prpl, un_top_fma_m4un_here, un_add_all, u_res;
    __m256d vn_here_minus_left, vn_here_minus_above, vn_here_mult_vn_here_minus_above, un_here_fma_vn_here_mult_v, vn_here_fma_mdtds_un_here;
    __m256d p_bot_minus_p_top, dt2rhods_mult_pbpt, vn_top_fma_m4vn_here, vn_add_all, v_res;
    __m256d constant_mdtds = _mm256_set1_pd(-dtds);
    __m256d constant_dt2rhods = _mm256_set1_pd(dt2rhods);
    __m256d constant_nu_dtsqds = _mm256_set1_pd(nu * dtsqds);
    __m256d constant_m4 = _mm256_set1_pd(-4.0);

    for(size_t i=1;i<d-1;i++){
        for(size_t j=1;j<d-4;j+=4){
            //const double un_here  = un[d*i     + j  ];
            //const double un_left  = un[d*i     + j-1];
            //const double un_right = un[d*i     + j+1];
            //const double un_below = un[d*(i+1) + j  ];
            //const double un_above = un[d*(i-1) + j  ];
            //const double vn_here  = vn[d*i     + j  ];
            //const double vn_left  = vn[d*i     + j-1];
            //const double vn_right = vn[d*i     + j+1];
            //const double vn_below = vn[d*(i+1) + j  ];
            //const double vn_above = vn[d*(i-1) + j  ];
            //const double p_left   =  p[d*i     + j-1];
            //const double p_right  =  p[d*i     + j+1];
            //const double p_below  =  p[d*(i+1) + j  ];
            //const double p_above  =  p[d*(i-1) + j  ];

            un_here = _mm256_loadu_pd(&un[d*i     + j   ]);
            un_l    = _mm256_loadu_pd(&un[d*i     + j-1 ]);
            un_r    = _mm256_loadu_pd(&un[d*i     + j+1]);
            un_bot  = _mm256_loadu_pd(&un[d*(i+1) + j  ]);
            un_top  = _mm256_loadu_pd(&un[d*(i-1) + j  ]);

           vn_here = _mm256_loadu_pd(&vn[d*i     + j   ]);
            vn_l    = _mm256_loadu_pd(&vn[d*i     + j-1 ]);
            vn_r    = _mm256_loadu_pd(&vn[d*i     + j+1]);
            vn_bot  = _mm256_loadu_pd(&vn[d*(i+1) + j  ]);
            vn_top  = _mm256_loadu_pd(&vn[d*(i-1) + j  ]);

           p_l    = _mm256_loadu_pd(&p[d*i     + j-1 ]);
            p_r    = _mm256_loadu_pd(&p[d*i     + j+1]);
            p_bot  = _mm256_loadu_pd(&p[d*(i+1) + j  ]);
            p_top  = _mm256_loadu_pd(&p[d*(i-1) + j  ]);


            un_here_minus_left = _mm256_sub_pd(un_here, un_l);
            un_here_minus_above = _mm256_sub_pd(un_here, un_top);
            vn_here_mult_un_here_minus_above = _mm256_mul_pd(vn_here, un_here_minus_above);
            un_here_fma_vn_here_mult = _mm256_fmadd_pd(un_here, un_here_minus_left, vn_here_mult_un_here_minus_above);
            un_here_fma_mdtds_un_here =  _mm256_fmadd_pd(constant_mdtds,  un_here_fma_vn_here_mult, un_here);
            //
            p_r_minus_p_l = _mm256_sub_pd(p_r, p_l);
            dt2rhods_mult_prpl = _mm256_mul_pd(constant_dt2rhods, p_r_minus_p_l);
//
           //
            un_top_fma_m4un_here = _mm256_fmadd_pd(un_here, constant_m4, un_top);
            un_add_all = _mm256_add_pd(un_r, un_l);
            un_add_all = _mm256_add_pd(un_add_all, un_bot);
            un_add_all = _mm256_add_pd(un_add_all, un_top_fma_m4un_here);
            
            u_res =_mm256_sub_pd(un_here_fma_mdtds_un_here, dt2rhods_mult_prpl);
            u_res = _mm256_fmadd_pd(constant_nu_dtsqds, un_add_all, u_res);

            _mm256_storeu_pd(&u[d*i+j], u_res);
//
            ////u[d*i+j] = un_here_fma_mdtds_un_here -
            ////        dt2rhods_mult_prpl +
            ////        nu * dtsqds *(un_add_all);
//
            ////u[d*i+j] = un_here
            ////            -dtds * (un_here * (un_here - un_left) +
            ////                    vn_here * (un_here - un_above)) -
            ////            dt2rhods * (p_right - p_left) +
            ////            nu * dtsqds *
            ////        (un_right + un_left + un_below + un_above - 4*un_here);
//
//
            vn_here_minus_left = _mm256_sub_pd(vn_here, vn_l);
            vn_here_minus_above = _mm256_sub_pd(vn_here, vn_top);
            vn_here_mult_vn_here_minus_above = _mm256_mul_pd(vn_here, vn_here_minus_above);
            un_here_fma_vn_here_mult_v = _mm256_fmadd_pd(un_here, vn_here_minus_left, vn_here_mult_vn_here_minus_above);
            vn_here_fma_mdtds_un_here =  _mm256_fmadd_pd(constant_mdtds,  un_here_fma_vn_here_mult_v, vn_here);
            
            p_bot_minus_p_top = _mm256_sub_pd(p_bot, p_top);
            dt2rhods_mult_pbpt = _mm256_mul_pd(constant_dt2rhods, p_bot_minus_p_top);
//
           //
            vn_top_fma_m4vn_here = _mm256_fmadd_pd(vn_here, constant_m4, vn_top);
            vn_add_all = _mm256_add_pd(vn_r, vn_l);
            vn_add_all = _mm256_add_pd(vn_add_all, vn_bot);
            vn_add_all = _mm256_add_pd(vn_add_all, vn_top_fma_m4vn_here);
//
            //
            v_res =_mm256_sub_pd(vn_here_fma_mdtds_un_here, dt2rhods_mult_pbpt);
            v_res = _mm256_fmadd_pd(constant_nu_dtsqds, vn_add_all, v_res);
//
            _mm256_storeu_pd(&v[d*i+j], v_res);
//
            //v[d*i+j] = vn_here_fma_mdtds_un_here
            //          -dt2rhods_mult_pbpt
            //          +nu * dtsqds *(vn_add_all);

            //v[d*i+j] = vn_here
            //            -dtds * (un_here * (vn_here - vn_left) +
            //                    vn_here * (vn_here - vn_above))
            //            -dt2rhods * (p_below - p_above)
            //            +nu * dtsqds *
            //            (vn_right + vn_left + vn_below + vn_above
            //                - 4*vn_here);
            
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
void advance_vanilla_avx_simulation(vanilla_avx_simulation* sim,
                 unsigned int steps,
                 unsigned int pit, double dt){
    for(unsigned int i=0;i<steps;i++){
        step_vanilla_avx_simulation(sim, pit, dt);
    }
} // Flops ?


void write_vanilla_avx_simulation(vanilla_avx_simulation* sim,
                   FILE* fp){
    write_preallocated_simulation(sim, fp);
}

void destroy_vanilla_avx_simulation(vanilla_avx_simulation* sim){
    destroy_preallocated_simulation(sim);
}
