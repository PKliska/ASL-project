#ifndef SIMULATION_H
#define SIMULATION_H
#include <stdio.h>

struct simulation{
    const struct simulation_vtable_* vtable_;
};

struct simulation_vtable_{
    void (*advance)(struct simulation* sim, unsigned steps, unsigned pit,
                    double dt);
    void (*write)(const struct simulation* sim, FILE* fp);
    void (*destroy)(struct simulation* sim);
};

void advance_simulation(struct simulation* sim, unsigned steps,
                        unsigned pit, double dt);


void write_simulation(const struct simulation* sim, FILE* fp);

void destroy_simulation(struct simulation* sim);

#endif
