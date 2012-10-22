#ifndef __PYSWAPMC_H__
#define __PYSWAPMC_H__

struct SimData{
    //lattice stuff
    int nsites;
    int nneighbors_per_site;
    int *neighbors;
    //simulation stuff
    double temp;
    double betaexp;
    float time;
    int *configuration;
    // rate stuff
    float *events;        
    float *event_refs;
    float *event_rates;
    float *event_ref_rates;
    // storage stuff
};

#endif
