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
    long current_step;
    int n_possible_events;
    int seed;
    double temp;
    double betaexp;
    double total_energy;
    double time;
    int *configuration;
    int *prev_configuration;
    int *dual_configuration;
    // rate stuff
    double total_rate;
    int n_event_types;
    int *events;        
    int *event_types;
    int *events_by_type;
    int *events_per_type;
    int *event_refs;
    double *event_rates;
//    float *event_ref_rates;
    double *cumulative_rates;
    // calculation stuff
    int *persistence_array;
};

#endif
