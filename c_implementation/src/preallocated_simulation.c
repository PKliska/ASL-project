#include <stdlib.h>
#include <stdio.h>
#include "utils.h"
#include "preallocated_simulation.h"

static const struct simulation_vtable_ PREALLOCATED_SIMULATION_VTABLE_[] = {{
    .advance=advance_preallocated_simulation,
    .write=write_preallocated_simulation,
    .destroy=destroy_preallocated_simulation
}};


struct preallocated_simulation* new_preallocated_simulation(
    size_t nx, size_t ny, double rho, double nu){
    struct preallocated_simulation *sim = malloc(
	sizeof(struct preallocated_simulation)
    );
    sim->base.vtable_ = PREALLOCATED_SIMULATION_VTABLE_;
    sim->nx = nx;
    sim->ny = ny;
    sim->rho = rho;
    sim->nu = nu;
    sim->u = zero_array(nx*ny);
    sim->un = zero_array(nx*ny);
    sim->v = zero_array(nx*ny);
    sim->vn = zero_array(nx*ny);
    sim->p = zero_array(nx*ny);
    sim->pn = zero_array(nx*ny);
    sim->b = zero_array(nx*ny);
    return sim;
}
static double sq(const double x){
    return x*x;
}

static void build_up_b(const struct preallocated_simulation* sim,
		       double dt){
    double *restrict b = sim->b;
    const size_t nx = sim->nx;
    const size_t ny = sim->ny;
    const double rho = sim->rho;
    const double dx = 2.0 / (nx - 1);
    const double dy = 2.0 / (ny - 1);
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
	    b[i*ny+j] = rho * (1 / dt *
		      ((u_right - u_left) / (2 * dx)
		      +(v_below - v_above) / (2 * dy))
		      -sq((u_right - u_left) / (2 * dx))
		      -2*((u_below - u_above) / (2 * dy)
			 *(v_right - v_left) / (2 * dx))
		      -sq((v_below - v_above) / (2 * dy)));
	}
    }
}

static void pressure_poisson(struct preallocated_simulation* sim,
			     unsigned int pit){
    double *restrict b = sim->b;
    const size_t nx = sim->nx;
    const size_t ny = sim->ny;
    const double dx = 2.0 / (nx - 1);
    const double dy = 2.0 / (ny - 1);
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
		p[i*ny+j] = ((pn_right + pn_left) * sq(dy)
	                    +(pn_below + pn_above) * sq(dx)) /
			    (2 * (sq(dx) + sq(dy))) -
			     sq(dx) * sq(dy) / (2 * (sq(dx) + sq(dy))) *
			     b[i*ny+j];
	    }
	}
	for(size_t i=0;i<nx;i++) p[ny*i      + ny-1] = p[ny*i + ny-2];
	for(size_t i=0;i<nx;i++) p[ny*i      + 0   ] = p[ny*i + 1];
	for(size_t j=0;j<ny;j++) p[ny*0      + j   ] = p[ny*1 + j];
	for(size_t j=0;j<ny;j++) p[ny*(nx-1) + j   ] = 0;
    }
}

/* Advance the simulation sim by one step of size dt using pit iterations
 * for the calculation of pressure. Use the array b for the calculation
 * of the intermediate matrix b.
*/
static void step_preallocated_simulation(
    struct preallocated_simulation* sim, unsigned int pit, double dt){
    
    const size_t nx = sim->nx;
    const size_t ny = sim->ny;
    const double dx = 2.0 / (nx - 1);
    const double dy = 2.0 / (ny - 1);
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
		         un_here * dt / dx *
		        (un_here - un_left) -
		         vn_here * dt / dy *
		        (un_here - un_above) -
		         dt / (2 * rho * dx) * (p_right - p_left) +
		         nu * (dt / sq(dx) *
			(un_right - 2 * un_here + un_left) +
			 dt / sq(dy) *
			(un_below - 2 * un_here + un_above)));
	    v[ny*i+j] = (vn_here -
			 un_here * dt / dx *
			(vn_here - vn_left) -
			 vn_here * dt / dy *
			(vn_here - vn_above) -
			 dt / (2 * rho * dy) * (p_below - p_above) +
		         nu * (dt / sq(dx) *
			(vn_right - 2 * vn_here + vn_left) +
			 dt / sq(dy) *
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
}
/* Advance the simulation sim by steps steps of size dt, using pit
 * iterations for the calculation of pressure.
*/
void advance_preallocated_simulation(struct preallocated_simulation* sim,
				 unsigned int steps,
				 unsigned int pit, double dt){
    for(unsigned int i=0;i<steps;i++){
	step_preallocated_simulation(sim, pit, dt);
    }
}


void write_preallocated_simulation(struct preallocated_simulation* sim,
			       FILE* fp){
    const size_t COLUMNS = sim->ny;
    fprintf(fp, "%zu,%zu,%lf,%lf,", sim->nx, sim->ny, sim->rho, sim->nu);

    for(size_t i=4;i<COLUMNS;i++) fputc(',', fp);
    fputc('\n', fp);

    write_matrix(sim->u, sim->nx, sim->ny, fp);

    for(size_t i=0;i<COLUMNS;i++) fputc(',', fp);
    fputc('\n', fp);

    write_matrix(sim->v, sim->nx, sim->ny, fp);

    for(size_t i=0;i<COLUMNS;i++) fputc(',', fp);
    fputc('\n', fp);

    write_matrix(sim->p, sim->nx, sim->ny, fp);
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
