#ifndef __PYSWAPMC_H__
#define __PYSWAPMC_H__

#define MAXZEROONEMODEL 10
#define MAXPLUSMINUSMODEL 20

struct SimData{
    //lattice stuff
    int nsites;
    int nneighbors_per_site;
    int *neighbors;
    //simulation stuff
    int model_number;
    int current_step;
    int n_possible_events;
    int seed;
    double temp;
    double betaexp;
    float time;
    int *configuration;
    // rate stuff
    float *events;        
    float *event_refs;
    float *event_rates;
    float *event_ref_rates;
    float *cumulative_rates;
    // storage stuff
    int *event_storage;
    float *persistence;
};

#endif
