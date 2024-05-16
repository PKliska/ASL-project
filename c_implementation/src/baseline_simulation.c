#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "utils.h"
#include "baseline_simulation.h"

static const struct simulation_vtable_ BASELINE_SIMULATION_VTABLE_[] = {{
    .advance=(void (*)(struct simulation *, unsigned int, unsigned int, double))advance_baseline_simulation,
    .write=(void (*)(const struct simulation *, FILE *))write_baseline_simulation,
    .destroy=(void (*)(struct simulation *))destroy_baseline_simulation
}};

struct baseline_simulation* new_baseline_simulation(
    size_t dimension, double size, double rho, double nu){
    struct baseline_simulation *sim = malloc(
	sizeof(struct baseline_simulation)
    );
    sim->base.vtable_ = BASELINE_SIMULATION_VTABLE_;
    sim->d = dimension;
    sim->rho = rho;
    sim->nu = nu;
    sim->size = size;
    sim->u = zero_array(dimension*dimension);
    sim->v = zero_array(dimension*dimension);
    sim->p = zero_array(dimension*dimension);
    return sim;
}
static double sq(const double x){
    return x*x;
}

static void build_up_b(const struct baseline_simulation* sim,
		       double* b, double dt){
    const size_t d = sim->d;
    const double rho = sim->rho;
    const double ds = sim->size / (d - 1);
    // const double ds = 0.025;
    const double* u = sim->u;
    const double* v = sim->v;

    for(size_t i=1; i < d-1; i++){
		for(size_t j=1; j < d-1; j++){
			const double u_left =  u[d*i     + j-1];
			const double u_right = u[d*i     + j+1];
			const double u_below = u[d*(i+1) + j  ];
			const double u_above = u[d*(i-1) + j  ];
			const double v_left =  v[d*i     + j-1];
			const double v_right = v[d*i     + j+1];
			const double v_below = v[d*(i+1) + j  ];
			const double v_above = v[d*(i-1) + j  ];

			b[i*d+j] = rho * (
				1/dt
				* ((u_right - u_left) / (2 * ds) + (v_below - v_above) / (2 * ds))
				- sq((u_right - u_left)  / (2 * ds))
				- 2*((u_below - u_above) / (2 * ds) *(v_right - v_left) / (2 * ds))
				- sq((v_below - v_above) / (2 * ds))
			);
			// FLOPS: 12 mul, 7 div, 1 add, 9 sub    (d-2)*(d-2)
		}
    }
}

static void pressure_poisson(struct baseline_simulation* sim,
			     unsigned int pit, double* b){
    const size_t d = sim->d;
    const double ds = sim->size / (d - 1);
    // const double ds = 0.025;
    const size_t bytes = d*d*sizeof(double);
    double* p = sim->p;
    double* pn = malloc(bytes);

    for(unsigned int q=0; q < pit; q++){
		memcpy(pn, p, bytes);

		for(size_t i=1; i < d-1; i++){
			for(size_t j=1; j < d-1; j++){
				const double pn_left =  pn[d*i     + j-1];
				const double pn_right = pn[d*i     + j+1];
				const double pn_below = pn[d*(i+1) + j  ];
				const double pn_above = pn[d*(i-1) + j  ];

				p[i*d+j] = (
					((pn_right + pn_left) * sq(ds)
								+(pn_below + pn_above) * sq(ds)) /
						(2 * (sq(ds) + sq(ds))) -
						sq(ds) * sq(ds) / (2 * (sq(ds) + sq(ds))) *
						b[i*d+j]
				);
			}
		} // FLOPS: 14 mul, 2 div, 5 add, 1 sub   (d-2)*(d-2)*pit

		// FLOPS: not counted as just data movement
		for(size_t i=0;i<d;i++) p[d*i      + d-1] = p[d*i + d-2];
		for(size_t i=0;i<d;i++) p[d*i      + 0  ] = p[d*i + 1];
		for(size_t j=0;j<d;j++) p[d*0      + j  ] = p[d*1 + j];
		for(size_t j=0;j<d;j++) p[d*(d-1)  + j  ] = 0;
	}
    free(pn);
}

/* Advance the simulation sim by one step of size dt using pit iterations
 * for the calculation of pressure. Use the array b for the calculation
 * of the intermediate matrix b.
*/
static void step_baseline_simulation(struct baseline_simulation* sim,
				     double *b, unsigned int pit,
				     double dt){
    const size_t d = sim->d;
    const double ds = sim->size / (d - 1); // FLOPS: 1 div
    // const double ds = 0.025; // FLOPS: 1 div
    const double rho = sim->rho;
    const double nu = sim->nu;
    double* u = sim->u;
    double* v = sim->v;
    double* p = sim->p;
    double* un = copy_array(sim->u, d*d);
    double* vn = copy_array(sim->v, d*d);

    build_up_b(sim, b, dt);        // FLOPS: 12 mul, 7 div, 1 add, 9 sub    (d-2)*(d-2)
    pressure_poisson(sim, pit, b); // FLOPS: 14 mul, 2 div, 5 add, 1 sub   (d-2)*(d-2)*pit

	for(size_t i=1; i < d-1; i++){
		for(size_t j=1; j < d-1; j++){
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
					nu * ( dt / sq(ds) *
				(un_right - 2 * un_here + un_left) +
				dt / sq(ds) *
				(un_below - 2 * un_here + un_above) ) );
			// FLOPS: 14 mul, 5 div, 8 sub, 4 add  (d-1)*(d-1)

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
			// FLOPS: 14 mul, 5 div, 8 sub, 4 add  (d-2)*(d-2)
		}
    }
	// FLOPS(both loops): 28 mul, 10 div, 16 sub, 8 add  (d-2)*(d-2)

	// FLOPS: not counted as just data movement
	for(size_t i=0;i<d;i++){
	    u[d*i + 0  ] = 0;
	    u[d*i + d-1] = 0;
	    v[d*i + 0  ] = 0;
	    v[d*i + d-1] = 0;
	}

	// FLOPS: not counted as just data movement
	for(size_t j=0;j<d;j++){
	    u[d*0     + j] = 0;
	    u[d*(d-1) + j] = 1; // cavity lid
	    v[d*0     + j] = 0;
	    v[d*(d-1) + j] = 0;
	}

    free(un);
    free(vn);
}

/* Advance the simulation sim by steps steps of size dt, using pit
 * iterations for the calculation of pressure.
*/
void advance_baseline_simulation(struct baseline_simulation* sim,
				 unsigned int steps,
				 unsigned int pit, double dt){
    const size_t d = sim->d;
    double* b = zero_array(d * d);
    for(unsigned int i=0; i < steps; i++){
		// FLOPS(one iteration): (40 + 14*pit) muls, (17 + 2*pit) divs, (9 + 5*pit) add, (25 + 1*pit) sub    (d-2)*(d-2)
	    step_baseline_simulation(sim, b, pit, dt);
    }
	// FLOPS(total): (40 + 14*pit) muls, (17 + 2*pit) divs, (9 + 5*pit) add, (25 + 1*pit) sub    (d-2)*(d-2)*steps
    free(b);
}


void write_baseline_simulation(struct baseline_simulation* sim, FILE* fp){
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

void destroy_baseline_simulation(struct baseline_simulation* sim){
    free(sim->u);
    free(sim->v);
    free(sim->p);
    free(sim);
}
