#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "pySpin.h"
#include "random.h"

int run_kmc_spin(int nsteps,struct SimData *SD){
    int step = 0;
    float t = SD->time;
    float dt = 0;
    for(step=0;step<nsteps;step++){
        //run advance simulation
        //dt = -log(prob)/total_rate
        t += dt;
    }
    SD->time = t;
    return 0;
}

int setup_spin_system(struct SimData *SD){
    return 0;
}

int cleanup_spin_system(struct SimData *SD){
    return 0;
}
