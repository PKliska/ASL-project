#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <immintrin.h>
#include <stdalign.h>
#include "utils.h"
#include "faster_math_simulation.h"
#include "aligned_avx_simulation.h"


static const struct simulation_vtable_ ALIGNED_AVX_SIMULATION_VTABLE[] = {{
    .advance=(void (*)(struct simulation *, unsigned int, unsigned int, double))advance_aligned_AVX_simulation,
    .write=(void (*)(const struct simulation *, FILE *))write_aligned_AVX_simulation,
    .destroy=(void (*)(struct simulation *))destroy_aligned_AVX_simulation
}};


aligned_AVX_simulation* new_aligned_AVX_simulation(
    size_t dimension, double size, double rho, double nu){
    aligned_AVX_simulation* sim = malloc(sizeof(*sim));
    init_aligned_AVX_simulation(sim, dimension, size, rho, nu);
    sim->base.base.vtable_ = ALIGNED_AVX_SIMULATION_VTABLE;
    return sim;
}
void init_aligned_AVX_simulation(aligned_AVX_simulation* sim,
    size_t dimension, double size, double rho, double nu){
    assert(dimension % 4 == 0);
    assert(dimension >= 12);
    init_faster_math_simulation(&sim->base, dimension, size, rho, nu);
    free(sim->base.p);
    free(sim->base.pn);
    free(sim->base.b);
    sim->base.p = aligned_zero_array(dimension*dimension);
    sim->base.pn = aligned_zero_array(dimension*dimension);
    sim->base.b = aligned_zero_array(dimension*dimension);
    sim->tmp = aligned_zero_array(dimension*dimension);
}
static double sq(const double x){
    return x*x;
}

static void detransform_data(double *restrict dst, double *restrict src,
                             size_t d){
    const size_t offset = d/4; // d must be divisible by 4
    for(unsigned r=0;r<d;r++){
        double* s_row = src + r*d;
        double* d_row = dst + r*d;
        for(unsigned t=0;t<d/4;t++){
            d_row[t+0*offset] = s_row[4*t+0];
            d_row[t+1*offset] = s_row[4*t+1];
            d_row[t+2*offset] = s_row[4*t+2];
            d_row[t+3*offset] = s_row[4*t+3];
        }
    }
}
static void transform_data(double *restrict dst, double *restrict src,
                             size_t d){
    const size_t offset = d/4; // d must be divisible by 4
    for(unsigned r=0;r<d;r++){
        double* s_row = src + r*d;
        double* d_row = dst + r*d;
        for(unsigned t=0;t<d/4;t++){
            d_row[4*t+0] = s_row[t+0*offset];
            d_row[4*t+1] = s_row[t+1*offset];
            d_row[4*t+2] = s_row[t+2*offset];
            d_row[4*t+3] = s_row[t+3*offset];
        }
    }
}

static void pressure_poisson(aligned_AVX_simulation* sim,
                 unsigned int pit){
    faster_math_simulation* base_sim = &sim->base;
    const size_t d = base_sim->d;
    
    // restore transformed pressure from tmp
    SWAP(double*, sim->tmp, base_sim->p);

    transform_data(sim->tmp, base_sim->b, d);
    SWAP(double*, sim->tmp, base_sim->b);
    double *restrict b = base_sim->b;

    for(unsigned int q=0;q<pit;q++){
        SWAP(double*, base_sim->p, base_sim->pn);
        double *restrict p = base_sim->p;
        double *restrict pn = base_sim->pn;

        for(size_t i=1;i<d-1;i++){
            {
                //unaligned
                const __m256d pn_left  = _mm256_loadu_pd(pn + d*i     + d-5);
                const __m256d pn_right = _mm256_load_pd (pn + d*i     + 4  );
                const __m256d pn_above = _mm256_load_pd (pn + d*(i-1) + 0  );
                const __m256d pn_below = _mm256_load_pd (pn + d*(i+1) + 0  );
                const __m256d b_value  = _mm256_load_pd (b  + d*i     + 0  );
                const __m256d res = (pn_right + pn_left + pn_below + pn_above - b_value) * 0.25;
                _mm256_store_pd(p + d*i, res);
            }
            for(size_t j=1;j<d/4-1;j++){ // d must be divisible by 4
                const __m256d pn_left  = _mm256_load_pd(pn + d*i     + 4*(j-1));
                const __m256d pn_right = _mm256_load_pd(pn + d*i     + 4*(j+1));
                const __m256d pn_below = _mm256_load_pd(pn + d*(i+1) + 4*j    );
                const __m256d pn_above = _mm256_load_pd(pn + d*(i-1) + 4*j    );
                const __m256d b_value  = _mm256_load_pd(b  + d*i     + 4*j    );
                const __m256d res = (pn_right + pn_left + pn_below + pn_above - b_value) * 0.25;
                _mm256_store_pd(p + d*i + 4*j, res);
            }
            {
                //unaligned
                const __m256d pn_right = _mm256_loadu_pd(pn + d*i     + 1  );
                const __m256d pn_left  = _mm256_load_pd (pn + d*i     + d-8);
                const __m256d pn_above = _mm256_load_pd (pn + d*(i-1) + d-4);
                const __m256d pn_below = _mm256_load_pd (pn + d*(i+1) + d-4);
                const __m256d b_value  = _mm256_load_pd (b  + d*i     + d-4);
                const __m256d res = (pn_right + pn_left + pn_below + pn_above - b_value) * 0.25;
                _mm256_store_pd(p + d*i + d-4, res);
            }
            p[d*i + 0  ] = p[d*i + 4];
            p[d*i + d-1] = p[d*i + d-5];
        }
        memcpy(p, p+d, d*sizeof(double));
        memset(p+d*(d-1), 0, d*sizeof(double));
    }
    detransform_data(sim->tmp, base_sim->p, d);
    //store transformed pressure in tmp, and make normal pressure available in p
    SWAP(double*, sim->tmp, base_sim->p);
}
// 22 * (d-2)(d-2)*pit + 2

/* Advance the simulation sim by one step of size dt using pit iterations
 * for the calculation of pressure. Use the array b for the calculation
 * of the intermediate matrix b.
*/
static void step_aligned_AVX_simulation(
    aligned_AVX_simulation* sim, unsigned int pit, double dt){
    struct preallocated_simulation* base_sim = &sim->base;
    faster_math_build_up_b(base_sim, dt);
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
    SWAP(double*, base_sim->u, base_sim->un);
    //Swap v and vn
    SWAP(double*, base_sim->v, base_sim->vn);

    double *restrict u = base_sim->u;
    double *restrict v = base_sim->v;
    double *restrict p = base_sim->p;
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
void advance_aligned_AVX_simulation(aligned_AVX_simulation* sim,
                 unsigned int steps,
                 unsigned int pit, double dt){
    for(unsigned int i=0;i<steps;i++){
        step_aligned_AVX_simulation(sim, pit, dt);
    }
} // Flops ?


void write_aligned_AVX_simulation(aligned_AVX_simulation* sim,
                   FILE* fp){
    write_faster_math_simulation(&sim->base, fp);
}

void destroy_aligned_AVX_simulation(aligned_AVX_simulation* sim){
    deinit_aligned_AVX_simulation(sim);
    free(sim);
}
void deinit_aligned_AVX_simulation(aligned_AVX_simulation* sim){
    deinit_faster_math_simulation(&sim->base);
    free(sim->tmp);
}
