#ifndef __PYSWAPMC_H__
#define __PYSWAPMC_H__

#define MAXZEROONEMODEL 9 // There are 10 of these
#define MAXPLUSMINUSMODEL 19 // There are 10 more of these

struct SimData{
    //lattice stuff
    int nsites;
    int nneighbors_per_site;
    int nneighbors_update_per_site;
    int *neighbors;
    int *neighbors_update;
    //simulation stuff
    int model_number;
    int current_step;
    int n_possible_events;
    int seed;
    double temp;
    double betaexp;
    float time;
    int *configuration;
    int *initial_configuration;
    int *dual_configuration;
    // rate stuff
    int *events;        
    int *event_refs;
    float *event_rates;
    float *event_ref_rates;
    float *cumulative_rates;
    // calculation stuff
    int *persistence_array;
};

#endif
