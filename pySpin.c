#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "pySpin.h"
#include "random.h"


int switch_state( int state, int model_number ){
    if(model_number<MAXZEROONEMODEL){ // models with 0's and 1's
        return 1 - state;
    }
    else if(model_number<MAXPLUSMINUSMODEL){ // models with -1's and 1's
        return -1*state;
    }
    else{
        printf("pySpin.c: When updating event, model number out of acceptable range\n");
        return -1;
    }
}

int update_events_i( int event_i, struct SimData *SD){
    int i,j,state_i;
    int n_possible_events = 0;
    float event_rate = 0.0;
    int model_number = SD->model_number;
    int nneighbors_per_site = SD->nneighbors_per_site;

    //first do it for this site
    int site_idx = SD->event_refs[event_i]; // need this below
    i = site_idx;
    state_i = SD->configuration[i];
    event_rate = get_event_rate(SD);
    SD->events[i] = switch_state( state_i, model_number );
    SD->event_rates[i] = event_rate;

    // now do it for neighbors
    for(j=0;j<nneighbors_per_site;j++){
        i = SD->neighbors[nneighbors_per_site*site_idx+j];
        state_i = SD->configuration[i];
        event_rate = get_event_rate(SD);
        SD->events[i] = switch_state( state_i, model_number );
        SD->event_rates[i] = event_rate;
    }
    // now update possible events
    for(i=0;i<SD->nsites;i++){
        SD->event_refs[i] = -1;
        SD->event_ref_rates[i] = 0.0;
        event_rate = SD->event_rates[i];
        if(event_rate>0){
            SD->event_refs[n_possible_events] = i;
            SD->event_ref_rates[n_possible_events] = event_rate;
            n_possible_events++;
        }
    }
    SD->n_possible_events = n_possible_events;
    return n_possible_events;
}

int update_configuration( int event_i, struct SimData *SD){
    int change_idx = SD->event_refs[event_i];
    int result = SD->events[change_idx];
    SD->configuration[change_idx] = result;
    //event_storage[SD->current_step] = change_idx;
    return 0;
}

double sum_rates( struct SimData *SD ){
    int i;
    double total_rate = 0.0;
    double cume_rate;
    if(SD->n_possible_events>0){
        SD->cumulative_rates[0] = SD->event_rates[0];
        cume_rate = SD->event_rates[0];
    }

    for(i=1;i<SD->n_possible_events;i++){
        total_rate = total_rate + SD->event_rates[i];
        cume_rate = cume_rate + SD->event_rates[i];
        SD->cumulative_rates[i] = cume_rate;
    }
    return total_rate;
} 

//binary search through cumulative event probabilities
int b_find_event( float searchval, struct SimData *SD){
    int final_event = SD->n_possible_events - 1;
    int idx0 = 0;
    int idx1 = final_event;
    if ( SD->n_possible_events < 1 ){
        return -1;
    }
    if ( SD->cumulative_rates[0] > searchval ){
        return 0;
    }
    if ( SD->cumulative_rates[final_event] < searchval ) {
        return final_event;
    }
    while ( (idx1-idx0) > 1 ){
        int halfidx = (idx1+idx0)/2;
        float half_prob = SD->cumulative_rates[halfidx];
        if ( half_prob > searchval ){
            idx1 = halfidx;
        }
        else if ( half_prob < searchval ){
            idx0 = halfidx;
        }
        else{
            return halfidx;
        }
    }
    return idx1;
}

int run_kmc_spin(int nsteps,struct SimData *SD){
    int step = 0;
    int return_val = 0;
    int n_possible_events = SD->n_possible_events;
    float t = SD->time;
    float dt = 0;
    for(step=0;step<nsteps;step++){
        if(n_possible_events<1){
            return_val = -1;
            break;
        }
        float prob = get_frandom();
        float total_rate = sum_rates(SD);
        int event_i = b_find_event( prob*total_rate, SD);
        update_configuration( event_i, SD );
        n_possible_events = update_events_i( event_i, SD );
        dt = -log(prob)/total_rate;
        t += dt;
    }
    SD->time = t;
    return return_val;
}

int setup_spin_system(struct SimData *SD){
    //all_events(SD);
    set_seed(SD->seed);
    return 0;
}

int cleanup_spin_system(struct SimData *SD){
    return 0;
}
