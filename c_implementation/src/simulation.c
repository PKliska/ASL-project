#include "simulation.h"

void advance_simulation(struct simulation* sim, unsigned steps,
                        unsigned pit, double dt){
    sim->vtable_->advance(sim, steps, pit, dt);
}


void write_simulation(const struct simulation* sim, FILE* fp){
    sim->vtable_->write(sim, fp);
}

void destroy_simulation(struct simulation* sim){
    sim->vtable_->destroy(sim);
}
