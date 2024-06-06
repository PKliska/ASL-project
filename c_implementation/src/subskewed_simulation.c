#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include "subskewed_simulation.h"
#include "faster_math_simulation.h"
#include "utils.h"
#include "trapeze_macros.h"

#ifndef SSKEWING_BLOCK_SIZE_X
#error "SSKEWING_BLOCK_SIZE_X was not defined, setting to default (36)"
#endif
#ifndef SSKEWING_BLOCK_SIZE_Y
#error "SSKEWING_BLOCK_SIZE_Y was not defined, setting to default (36)"
#endif
#ifndef SSKEWING_TIMESTEPS
#error "SSKEWING_TIMESTEPS was not defined, setting to default (10)"
#endif
#ifndef SSKEWING_SUBBLOCK_SIZE_X
#error "SSKEWING_SUBBLOCK_SIZE_X was not defined"
#endif
#ifndef SSKEWING_SUBBLOCK_SIZE_Y
#error "SSKEWING_SUBBLOCK_SIZE_Y was not defined"
#endif
#if SSKEWING_BLOCK_SIZE_X % SSKEWING_SUBBLOCK_SIZE_X != 0\
|| SSKEWING_BLOCK_SIZE_Y % SSKEWING_SUBBLOCK_SIZE_Y != 0
#error "SSKEWING_SUBBLOCK_SIZE must divide the SSKEWING_BLOCK_SIZE"
#endif
#if SSKEWING_SUBBLOCK_SIZE_X <= SSKEWING_TIMESTEPS \
 || SSKEWING_SUBBLOCK_SIZE_Y <= SSKEWING_TIMESTEPS
#error "SSKEWING_SUBBLOCK_SIZE must be greater than SSKEWING_TIMESTEPS"
#endif



static const struct simulation_vtable_ subskewed_SIMULATION_VTABLE[] = {{
    .advance=(void (*)(struct simulation *, unsigned int, unsigned int, double))advance_subskewed_simulation,
    .write=(void (*)(const struct simulation *, FILE *))write_subskewed_simulation,
    .destroy=(void (*)(struct simulation *))destroy_subskewed_simulation
}};


subskewed_simulation* new_subskewed_simulation(
    size_t dimension, double size, double rho, double nu){
    subskewed_simulation* sim = malloc(sizeof(*sim));
    init_subskewed_simulation(sim, dimension, size, rho, nu);
    return sim;
}

void init_subskewed_simulation(subskewed_simulation* sim,
    size_t dimension, double size, double rho, double nu){
    init_faster_math_simulation(sim, dimension, size, rho, nu);
    assert(dimension % SSKEWING_SUBBLOCK_SIZE_X == 0);
    assert(dimension % SSKEWING_SUBBLOCK_SIZE_Y == 0);
    sim->base.vtable_ = subskewed_SIMULATION_VTABLE;
}


static double sq(const double x){
    return x*x;
}

__attribute__((always_inline))
static inline void subskewed_trapeze(double *b, double *p0, double* p1, size_t d,
                                     unsigned ts,
                                     int x0, int dx0, int x1, int dx1,
                                     int y0, int dy0, int y1, int dy1){
    if(dx0 == 0 && dy0 == 0){ // top left
        DO_TRAPEZE_TOP_LEFT(b, p0, p1, d,
                            1, ts+1,
                            x0+1, 0, x1, dx1,
                            y0+1, 0, y1, dy1);
    }else if(dx0 == 0 && dy1 == 0){ // top right
        DO_TRAPEZE_TOP_RIGHT(b, p0, p1, d,
                             1, ts+1,
                             x0+1, 0, x1, dx1,
                             y0, dy0, y1-1, 0);
    }else if(dx0 == 0){ // top
        DO_TRAPEZE_TOP(b, p0, p1, d,
                       1, ts+1,
                       x0+1, 0, x1, dx1,
                       y0, dy0, y1, dy1);
    }else if(dx1 == 0 && dy0 == 0){ // bottom left
        DO_TRAPEZE_BOTTOM_LEFT(b, p0, p1, d,
                       1, ts+1,
                       x0, dx0, x1-1, 0,
                       y0+1, 0, y1, dy1);

    }else if(dx1 == 0 && dy1 == 0){ //bottom right
        DO_TRAPEZE_BOTTOM_RIGHT(b, p0, p1, d,
                       1, ts+1,
                       x0, dx0, x1-1, 0,
                       y0, dy0, y1-1, 0);
    }else if(dx1 == 0){ // bottom
        DO_TRAPEZE_BOTTOM(b, p0, p1, d,
                       1, ts+1,
                       x0, dx0, x1-1, 0,
                       y0, dy0, y1, dy1);

    }else if(dy0 == 0){ // left
        DO_TRAPEZE_LEFT(b, p0, p1, d,
                       1, ts+1,
                       x0, dx0, x1, dx1,
                       y0+1, 0, y1, dy1);
    }else if(dy1 == 0){ // right
        DO_TRAPEZE_RIGHT(b, p0, p1, d,
                       1, ts+1,
                       x0, dx0, x1, dx1,
                       y0, dy0, y1-1, 0);
    }else{ // mid
        DO_TRAPEZE_MID(b, p0, p1, d,
                       1, ts+1,
                       x0, dx0, x1, dx1,
                       y0, dy0, y1, dy1);
    }
}

__attribute__((always_inline))
static inline void subskewed_layer1(double* b, double* p0, double* p1, size_t d,
                                      unsigned ts,
                                      int bx0, int dx0, int bx1, int dx1,
                                      int by0, int dy0, int by1, int dy1){
    int x0,y0,x1,y1;
    //top left
    x0 = bx0*SSKEWING_SUBBLOCK_SIZE_X; x1 = (bx0+1)*SSKEWING_SUBBLOCK_SIZE_X;
    y0 = by0*SSKEWING_SUBBLOCK_SIZE_Y; y1 = (by0+1)*SSKEWING_SUBBLOCK_SIZE_Y;
    subskewed_trapeze(b, p0, p1, d,
                      ts,
                      x0, dx0, x1, -1,
                      y0, dy0, y1, -1);
    //top
    for(int by = by0+1; by < by1-1; by++){
        int y0 = by*SSKEWING_SUBBLOCK_SIZE_Y, y1 = (by+1)*SSKEWING_SUBBLOCK_SIZE_Y;
        subskewed_trapeze(b, p0, p1, d,
                          ts,
                          x0, dx0, x1, -1,
                          y0, -1, y1, -1);
    }
    //top right
    y0 = (by1-1)*SSKEWING_SUBBLOCK_SIZE_Y, y1 = by1*SSKEWING_SUBBLOCK_SIZE_Y;
    subskewed_trapeze(b, p0, p1, d,
                      ts,
                      x0, dx0, x1, -1,
                      y0, -1, y1, dy1);
    for(int bx = bx0+1; bx < bx1-1; bx++){
        //left
        int x0 = bx*SSKEWING_SUBBLOCK_SIZE_X, x1 = (bx+1)*SSKEWING_SUBBLOCK_SIZE_X;
        int y0 = by0*SSKEWING_SUBBLOCK_SIZE_Y, y1 = (by0+1)*SSKEWING_SUBBLOCK_SIZE_Y;
        subskewed_trapeze(b, p0, p1, d,
                          ts,
                          x0, -1, x1, -1,
                          y0, dy0, y1, -1);
        for(int by = by0+1; by < by1-1; by++){
            //mid
            int y0 = by*SSKEWING_SUBBLOCK_SIZE_Y, y1 = (by+1)*SSKEWING_SUBBLOCK_SIZE_Y;
            subskewed_trapeze(b, p0, p1, d,
                              ts,
                              x0, -1, x1, -1,
                              y0, -1, y1, -1);
        }
        //right
        y0 = (by1-1)*SSKEWING_SUBBLOCK_SIZE_Y, y1 = by1*SSKEWING_SUBBLOCK_SIZE_Y;
        subskewed_trapeze(b, p0, p1, d,
                          ts,
                          x0, -1, x1, -1,
                          y0, -1, y1, dy1);
    }
    //bottom left
    x0 = (bx1-1)*SSKEWING_SUBBLOCK_SIZE_X; x1 = bx1*SSKEWING_SUBBLOCK_SIZE_X;
    y0 = by0*SSKEWING_SUBBLOCK_SIZE_Y; y1 = (by0+1)*SSKEWING_SUBBLOCK_SIZE_Y;
    subskewed_trapeze(b, p0, p1, d,
                      ts,
                      x0, -1, x1, dx1,
                      y0, dy0, y1, -1);
    for(int by = by0+1; by < by1-1; by++){
        //bottom
        int y0 = by*SSKEWING_SUBBLOCK_SIZE_Y, y1 = (by+1)*SSKEWING_SUBBLOCK_SIZE_Y;
        subskewed_trapeze(b, p0, p1, d,
                          ts,
                          x0, -1, x1, dx1,
                          y0, -1, y1, -1);
    }
    //bottom right
    y0 = (by1-1)*SSKEWING_SUBBLOCK_SIZE_Y, y1 = by1*SSKEWING_SUBBLOCK_SIZE_Y;
    subskewed_trapeze(b, p0, p1, d,
                      ts,
                      x0, -1, x1, dx1,
                      y0, -1, y1, dy1);
}
__attribute__((always_inline))
static inline void subskewed_layer2(double* b, double* p0, double* p1, size_t d,
                                      unsigned ts){
    const int bbx0 = 0, bbx1 = d/SSKEWING_BLOCK_SIZE_X; 
    const int bby0 = 0, bby1 = d/SSKEWING_BLOCK_SIZE_Y; 
    const int SBPB_X = SSKEWING_BLOCK_SIZE_X/SSKEWING_SUBBLOCK_SIZE_X;
    const int SBPB_Y = SSKEWING_BLOCK_SIZE_Y/SSKEWING_SUBBLOCK_SIZE_Y;

    int bx0,by0,bx1,by1;
    //top left
    {
    int bx0 = bbx0*SBPB_X, bx1 = (bbx0+1)*SBPB_X;
    int by0 = bby0*SBPB_Y, by1 = (bby0+1)*SBPB_Y;
    subskewed_layer1(b, p0, p1, d,
                      ts,
                      bx0, 0, bx1, -1,
                      by0, 0, by1, -1);
    //top
    for(int bby = bby0+1; bby < bby1-1; bby++){
        int by0 = bby*SBPB_Y, by1 = (bby+1)*SBPB_Y;
        subskewed_layer1(b, p0, p1, d,
                         ts,
                         bx0, 0, bx1, -1,
                         by0, -1, by1, -1);
    }
    //top right
    by0 = (bby1-1)*SBPB_Y, by1 = bby1*SBPB_Y;
    subskewed_layer1(b, p0, p1, d,
                      ts,
                      bx0, 0, bx1, -1,
                      by0, -1, by1, 0);
    }
    for(int bbx = bbx0+1; bbx < bbx1-1; bbx++){
        //left
        int bx0 = bbx*SBPB_X, bx1 = (bbx+1)*SBPB_X;
        int by0 = bby0*SBPB_Y, by1 = (bby0+1)*SBPB_Y;
        subskewed_layer1(b, p0, p1, d,
                          ts,
                          bx0, -1, bx1, -1,
                          by0, 0, by1, -1);
        for(int bby = bby0+1; bby < bby1-1; bby++){
            //mid
            int by0 = bby*SBPB_Y, by1 = (bby+1)*SBPB_Y;
            subskewed_layer1(b, p0, p1, d,
                              ts,
                              bx0, -1, bx1, -1,
                              by0, -1, by1, -1);
        }
        //right
        by0 = (bby1-1)*SBPB_Y, by1 = bby1*SBPB_Y;
        subskewed_layer1(b, p0, p1, d,
                          ts,
                          bx0, -1, bx1, -1,
                          by0, -1, by1, 0);
    }
    {
    //bottom left
    int bx0 = (bbx1-1)*SBPB_X, bx1 = bbx1*SBPB_X;
    int by0 = bby0*SBPB_Y, by1 = (bby0+1)*SBPB_Y;
    subskewed_layer1(b, p0, p1, d,
                      ts,
                      bx0, -1, bx1, 0,
                      by0, 0, by1, -1);
    for(int bby = bby0+1; bby < bby1-1; bby++){
        //bottom
        int by0 = bby*SBPB_Y, by1 = (bby+1)*SBPB_Y;
        subskewed_layer1(b, p0, p1, d,
                          ts,
                          bx0, -1, bx1, 0,
                          by0, -1, by1, -1);
    }
    //bottom right
    by0 = (bby1-1)*SBPB_Y, by1 = bby1*SBPB_Y;
    subskewed_layer1(b, p0, p1, d,
                      ts,
                      bx0, -1, bx1, 0,
                      by0, -1, by1, 0);
    }
}

static void subskewed_pressure_poisson(subskewed_simulation* sim,
                 unsigned int pit){
    const size_t d = sim->d;
    double *restrict b = sim->b;
    double *restrict p0 = sim->p, *restrict p1 = sim->pn;
    
    const int time_blocks = (pit+SSKEWING_TIMESTEPS-1) / SSKEWING_TIMESTEPS;
    const int last_time = (pit-1) % SSKEWING_TIMESTEPS + 1;
    int bt;
    for(bt=0;bt + 1 < time_blocks - 1; bt += 2){
        subskewed_layer2(b, p0, p1, d, SSKEWING_TIMESTEPS);
        if(SSKEWING_TIMESTEPS % 2 == 1) SWAP(double*, p0, p1);
        subskewed_layer2(b, p0, p1, d, SSKEWING_TIMESTEPS);
        if(SSKEWING_TIMESTEPS % 2 == 1) SWAP(double*, p0, p1);
    }
    if(bt < time_blocks - 1){
        subskewed_layer2(b, p0, p1, d, SSKEWING_TIMESTEPS);
        if(SSKEWING_TIMESTEPS % 2 == 1) SWAP(double*, p0, p1);
    }
    subskewed_layer2(b, p0, p1, d, last_time);
    if(SSKEWING_TIMESTEPS % 2 == 1) SWAP(double*, p0, p1);
    
    sim->p = p0; sim->pn = p1;

    sim->p[0] = sim->p[1];
    sim->p[d-1] = sim->p[d-2];
    sim->p[(d-1)*d+0] = sim->p[(d-1)*d+1];
    sim->p[(d-1)*d+d-1] = sim->p[(d-1)*d+d-2];
}
// 22 * (d-2)(d-2)*pit + 2

/* Advance the simulation sim by one step of size dt using pit iterations
 * for the calculation of pressure. Use the array b for the calculation
 * of the intermediate matrix b.
*/
static void step_subskewed_simulation(
    subskewed_simulation* sim, unsigned int pit, double dt){

    faster_math_build_up_b(sim, dt);
    subskewed_pressure_poisson(sim, pit);
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
void advance_subskewed_simulation(subskewed_simulation* sim,
                 unsigned int steps,
                 unsigned int pit, double dt){
    for(unsigned int i=0;i<steps;i++){
        step_subskewed_simulation(sim, pit, dt);
    }
} // Flops ?


void write_subskewed_simulation(subskewed_simulation* sim,
                   FILE* fp){
    write_preallocated_simulation(sim, fp);
}

void destroy_subskewed_simulation(subskewed_simulation* sim){
    deinit_subskewed_simulation(sim);
    free(sim);
}

void deinit_subskewed_simulation(subskewed_simulation *sim){
    deinit_preallocated_simulation(sim);
}
