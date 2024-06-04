from dataclasses import dataclass
import argparse
from pathlib import Path


@dataclass
class Trapeze:
    t0: int
    t1: int
    x0: int
    dx0: int
    x1: int
    dx1: int
    y0: int
    dy0: int
    y1: int
    dy1: int

    def __contains__(self, tup):
        t, x, y = tup
        if not self.t0 <= t < self.t1:
            return False
        x0 = self.x0 + (t-self.t0)*self.dx0
        x1 = self.x1 + (t-self.t0)*self.dx1
        y0 = self.y0 + (t-self.t0)*self.dy0
        y1 = self.y1 + (t-self.t0)*self.dy1
        return x0 <= x < x1 and y0 <= y < y1

    def __str__(self):
        pmz = {-1: 'm', 0: 'z', 1: 'p'}
        return (f"t{self.t1 - self.t0}"
                f"w{self.y1 - self.y0}"
                f"h{self.x1 - self.x0}"
                f"d{pmz[self.dx0]}{pmz[self.dx1]}{pmz[self.dy0]}{pmz[self.dy1]}")


def iterate_trapeze(trap: Trapeze):
    for t in range(trap.t0, trap.t1):
        x0 = trap.x0 + trap.dx0 * (t - trap.t0)
        x1 = trap.x1 + trap.dx1 * (t - trap.t0)
        y0 = trap.y0 + trap.dy0 * (t - trap.t0)
        y1 = trap.y1 + trap.dy1 * (t - trap.t0)
        for x in range(x0, x1):
            for y in range(y0, y1):
                yield (t, x, y)


def safe_int(x):
    if x < 0:
        return f"m{abs(x)}"
    else:
        return x


def compute_s(computed: set, code: list, trap: Trapeze, t, x, y):
    res = f"s_{t}_{safe_int(x)}_{safe_int(y)}"
    if res in computed:
        return res
    first = compute_p(computed, code, trap, t, x, y)
    second = compute_p(computed, code, trap, t, x+1, y+1)
    computed.add(res)
    code.append(f"\tdouble {res} = {first} + {second};\n")
    return res


def compute_p(computed: set, code: list, trap: Trapeze, t, x, y):
    res = f"p_{t}_{safe_int(x)}_{safe_int(y)}"
    if res in computed:
        return res
    if (t, x, y) not in trap:
        computed.add(res)
        code.append(f"\tdouble {res} = p{t % 2}[(x0 + ({x}))*d + y0 + ({y})];\n")
        return res
    if x == trap.x0 and trap.dx0 == 0:
        other = compute_p(computed, code, trap, t, x+1, y)
        computed.add(res)
        code.append(f"\tdouble {res} = {other};\n")
        return res
    if x == trap.x1 - 1 and trap.dx1 == 0:
        other = compute_p(computed, code, trap, t, x-1, y)
        computed.add(res)
        code.append(f"\tdouble {res} = 0;\n")
        return res
    if y == trap.y0 and trap.dy0 == 0:
        other = compute_p(computed, code, trap, t, x, y+1)
        computed.add(res)
        code.append(f"\tdouble {res} = {other};\n")
        return res
    if y == trap.y1 - 1 and trap.dy1 == 0:
        other = compute_p(computed, code, trap, t, x, y-1)
        computed.add(res)
        code.append(f"\tdouble {res} = {other};\n")
        return res
    left = compute_s(computed, code, trap, t-1, x, y-1)
    top = compute_s(computed, code, trap, t-1, x-1, y)
    computed.add(res)
    code.append(f"\tdouble {res} = ({left} + {top} - b[(x0 + ({x}))*d + y0 + ({y})]) * 0.25;\n")
    return res


def generate_offset_trapeze(trap: Trapeze):
    computed = set()
    code = []
    last_assignment = {}
    name = f"generated_offset_trapeze_{trap}"
    signature = ( "__attribute__((always_inline))\n"
                 f"inline void {name}"
                  "(size_t x0, size_t y0, size_t d, "
                  "double *restrict b, "
                  "double *restrict p0, double *restrict p1)")
    code.append(signature)
    code.append("{\n")
    for (t, x, y) in iterate_trapeze(trap):
        res = compute_p(computed, code, trap, t, x, y)
        last_assignment[(t % 2, x, y)] = res
    for (t, x, y), var in last_assignment.items():
        code.append(f"\tp{t}[(x0 + ({x}))*d + y0 + {y}] = {var};\n")
    code.append("}\n")
    return (name, ''.join(code))


def generate_pressure_poisson(block_size, time_step):
    funcs = []
    code = []
    name = f"pressure_poisson_b{block_size}t{time_step}"
    signature = f"void {name}(double *restrict b, double *restrict p0, double *restrict p1, size_t d)"
    code.append(signature)
    code.append("{\n")
    """ Top left """
    f_name, f_code = generate_offset_trapeze(Trapeze(1, time_step+1,
                                                     0, 0, block_size, -1,
                                                     0, 0, block_size, -1))
    funcs.append(f_code)
    code.append(f"\t{f_name}(0, 0, d, b, p0, p1);\n")
    code.append(f"\tfor(size_t j=1;j<d/{block_size}-1;j++)")
    code.append("{\n")
    """ Top """
    f_name, f_code = generate_offset_trapeze(Trapeze(1, time_step+1,
                                                     0, 0, block_size, -1,
                                                     0, -1, block_size, -1))
    funcs.append(f_code)
    code.append(f"\t\t{f_name}(0, j*{block_size}, d, b, p0, p1);\n")
    code.append("\t}\n")
    """ Top right """
    f_name, f_code = generate_offset_trapeze(Trapeze(1, time_step+1,
                                                     0, 0, block_size, -1,
                                                     0, -1, block_size, 0))
    funcs.append(f_code)
    code.append(f"\t{f_name}(0, d-{block_size}, d, b, p0, p1);\n")
    """ Left """
    code.append(f"\tfor(size_t i=1;i<d/{block_size}-1;i++)")
    code.append("{\n")
    f_name, f_code = generate_offset_trapeze(Trapeze(1, time_step+1,
                                                     0, -1, block_size, -1,
                                                     0, 0, block_size, -1))
    funcs.append(f_code)
    code.append(f"\t\t{f_name}(i*{block_size}, 0, d, b, p0, p1);\n")
    code.append(f"\t\tfor(size_t j=1;j<d/{block_size}-1;j++)")
    code.append("{\n")
    """ Mid """
    f_name, f_code = generate_offset_trapeze(Trapeze(1, time_step+1,
                                                     0, -1, block_size, -1,
                                                     0, -1, block_size, -1))
    funcs.append(f_code)
    code.append(f"\t\t\t{f_name}(i*{block_size}, j*{block_size}, d, b, p0, p1);\n")
    code.append("\t\t}\n")
    """ Right """
    f_name, f_code = generate_offset_trapeze(Trapeze(1, time_step+1,
                                                     0, -1, block_size, -1,
                                                     0, -1, block_size, 0))
    funcs.append(f_code)
    code.append(f"\t\t{f_name}(i*{block_size}, d-{block_size}, d, b, p0, p1);\n")
    code.append("\t}\n")
    """ Bottom left """
    f_name, f_code = generate_offset_trapeze(Trapeze(1, time_step+1,
                                                     0, -1, block_size, 0,
                                                     0, 0, block_size, -1))
    funcs.append(f_code)
    code.append(f"\t{f_name}(d-{block_size}, 0, d, b, p0, p1);\n")
    """ Bottom """
    code.append(f"\tfor(size_t j=1;j<d/{block_size}-1;j++)")
    code.append("{\n")
    f_name, f_code = generate_offset_trapeze(Trapeze(1, time_step+1,
                                                     0, -1, block_size, 0,
                                                     0, -1, block_size, -1))
    funcs.append(f_code)
    code.append(f"\t\t{f_name}(d-{block_size}, j*{block_size}, d, b, p0, p1);\n")
    code.append("\t}\n")
    """ Bottom right """
    f_name, f_code = generate_offset_trapeze(Trapeze(1, time_step+1,
                                                     0, -1, block_size, 0,
                                                     0, -1, block_size, 0))
    funcs.append(f_code)
    code.append(f"\t{f_name}(d-{block_size}, d-{block_size}, d, b, p0, p1);\n")

    code.append("}\n")
    return (name, ''.join(funcs + code))


parser = argparse.ArgumentParser()
parser.add_argument("output_path", type=Path, default='.')
parser.add_argument("block_size", type=int, default=16)
parser.add_argument("time_step", type=int, default=10)
args = parser.parse_args()

source_path = args.output_path / 'src' / 'generated_simulation.c'
header_path = args.output_path / 'include' / 'generated_simulation.h'
source_path.parent.mkdir(exist_ok=True)
header_path.parent.mkdir(exist_ok=True)

header_content = f"""
#ifndef GENERATED_H
#define GENERATED_H
#include <stdio.h>
#include "faster_math_simulation.h"
#define GENERATED_BLOCK_SIZE {args.block_size}
#define GENERATED_TIME_STEP {args.time_step}
""" + """

typedef faster_math_simulation generated_simulation;

generated_simulation* new_generated_simulation(
    size_t dimension, double size, double rho, double nu);

void init_generated_simulation(generated_simulation* sim,
    size_t dimension, double size, double rho, double nu);

void advance_generated_simulation(generated_simulation* sim,
                                 unsigned int steps, unsigned int pit,
                                 double dt);

void write_generated_simulation(generated_simulation* sim,
                               FILE* fp);

void destroy_generated_simulation(generated_simulation* sim);
void deinit_generated_simulation(generated_simulation* sim);

#endif
"""
source_content = """
#include <stdlib.h>
#include <assert.h>
#include "utils.h"
#include "generated_simulation.h"

static const struct simulation_vtable_ generated_SIMULATION_VTABLE[] = {{
    .advance=(void (*)(struct simulation *, unsigned int, unsigned int, double))advance_generated_simulation,
    .write=(void (*)(const struct simulation *, FILE *))write_generated_simulation,
    .destroy=(void (*)(struct simulation *))destroy_generated_simulation
}};


generated_simulation* new_generated_simulation(
    size_t dimension, double size, double rho, double nu){
    generated_simulation* sim = malloc(sizeof(*sim));
    init_generated_simulation(sim, dimension, size, rho, nu);
    return sim;
}

void init_generated_simulation(generated_simulation* sim,
    size_t dimension, double size, double rho, double nu){
    init_faster_math_simulation(sim, dimension, size, rho, nu);
    assert(dimension % GENERATED_BLOCK_SIZE == 0);
    sim->base.vtable_ = generated_SIMULATION_VTABLE;
}


static double sq(const double x){
    return x*x;
}

<DEFINE FUNCTIONS>

static void generated_pressure_poisson(generated_simulation* sim,
                 unsigned int pit){
    const size_t d = sim->d;
    double *restrict b = sim->b;
    double *restrict p0 = sim->p;
    double *restrict p1 = sim->pn;
    for(int i=0;i < pit / GENERATED_TIME_STEP;i++){
        <CALL FUNCTION>
        if(GENERATED_TIME_STEP % 2 == 1){
            SWAP(double*, p0, p1);
        }
    }
}
// 22 * (d-2)(d-2)*pit + 2

/* Advance the simulation sim by one step of size dt using pit iterations
 * for the calculation of pressure. Use the array b for the calculation
 * of the intermediate matrix b.
*/
static void step_generated_simulation(
    generated_simulation* sim, unsigned int pit, double dt){
    faster_math_build_up_b(sim, dt);
    generated_pressure_poisson(sim, pit);
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
}

/* Advance the simulation sim by steps steps of size dt, using pit
 * iterations for the calculation of pressure.
*/
void advance_generated_simulation(generated_simulation* sim,
                 unsigned int steps,
                 unsigned int pit, double dt){
    assert(pit % GENERATED_TIME_STEP == 0);
    for(unsigned int i=0;i<steps;i++){
        step_generated_simulation(sim, pit, dt);
    }
}


void write_generated_simulation(generated_simulation* sim,
                   FILE* fp){
    write_preallocated_simulation(sim, fp);
}

void destroy_generated_simulation(generated_simulation* sim){
    deinit_generated_simulation(sim);
    free(sim);
}

void deinit_generated_simulation(generated_simulation *sim){
    deinit_preallocated_simulation(sim);
}
"""

name, generated_fns = generate_pressure_poisson(args.block_size, args.time_step)
source_content = source_content.replace("<CALL FUNCTION>", f"{name}(b, p0, p1, d);", 1)
source_content = source_content.replace("<DEFINE FUNCTIONS>", generated_fns, 1)
with source_path.open('w') as source, header_path.open('w') as header:
    source.write(source_content)
    header.write(header_content)
