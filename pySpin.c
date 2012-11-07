#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "pySpin.h"
#include "random.h"

int get_event_type(int site_idx, struct SimData *SD){
    int j;
    int event_type = -1;
    int nneighbors_per_site = SD->nneighbors_per_site;
    //float event_rate = 0.0;
    if(SD->model_number < 2){ // FA or East
        float constraint = 0.0;
        int state_i = SD->configuration[site_idx];
        //note this depends on the lattice being properly declared such that 
        //  all neighbors are those that actually affect the given spin
        for(j=0;j<nneighbors_per_site;j++){
            constraint+=SD->configuration[SD->neighbors[nneighbors_per_site*site_idx+j]];
        }
        event_type = 2*constraint+state_i;
//        event_rate = (1-state_i)*constraint*SD->betaexp + state_i*constraint;
//        printf("%i %i %f %f %f\n",site_idx,state_i,constraint,SD->betaexp,event_rate);
    }
    else if(SD->model_number==10){ // Plaquette
        int excitations_created = 0;
        excitations_created += ( 1-2*SD->dual_configuration[site_idx] ); // changes 0 -> 1 and 1-> -1
        // note, on the next line it is just the first nneighbors_per_site out of the neighbors_update list which have this spin in their plaquette. e.g. 3 more spins out of eight neighbors for square plaquette
        for(j=0;j<nneighbors_per_site;j++){
            int neighbor_site = SD->neighbors_update[SD->nneighbors_update_per_site*site_idx+j];
            excitations_created += ( 1 - 2*SD->dual_configuration[neighbor_site] );
        }
        //metropolis rate
        //event_rate = (excitations_created>0)*pow(SD->betaexp,excitations_created) + (excitations_created<=0); // if n plaquettes are excited, rate is exp(-n * beta). otherwise rate 1 (decreases energy or keeps same)

//WARNING, next is for square plaquette only. need to generalize

        event_type = excitations_created/2+2;
    }
    else{
       printf("pySpin.c: Model number not yet supported by get_event_type\n");
    }
    return event_type;
}

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


int all_events( struct SimData *SD){
    int i,j;
    int n_possible_events = 0;
    for(i=0;i<SD->nsites;i++){
        SD->events[i] = 0;
        SD->event_types[i] = -1;
        SD->event_refs[i] = -1;
        for(j=0;j<SD->n_event_types;j++){
            SD->events_by_type[j*SD->nsites+i]=0;
        }
    }
    for(i=0;i<SD->n_event_types;i++){
//        SD->event_rates[i] = 0.0;
        SD->events_per_type[i] = 0.0;
        SD->cumulative_rates[i] = 0;
    }

    for(i=0;i<SD->nsites;i++){
//        float event_rate = get_event_rate(i,SD);
        int event_type = get_event_type(i,SD);
        SD->events[i] = switch_state(SD->configuration[i],SD->model_number);
        SD->event_types[i] = event_type;
        if(event_type>0){
            int nevents_type_i = SD->events_per_type[event_type];
            SD->events_by_type[event_type*SD->nsites+nevents_type_i] = i;
            SD->event_refs[i] = nevents_type_i;
            SD->events_per_type[event_type]++;
            n_possible_events++;
        }
    }
    SD->n_possible_events = n_possible_events;
    for(i=0;i<SD->n_event_types;i++){
        printf("type %i) rate: %f nevents: %i\n",i,SD->event_rates[i],SD->events_per_type[i]);
    }
    return n_possible_events;
}

int change_event_type( int i, int old_event_type, int new_event_type, struct SimData *SD){
//first remove it from the old list (note, this hopefully should work even if last event in list)
    int old_final_event_ref = SD->events_per_type[old_event_type];
    int old_event_ref = SD->events_per_type[old_event_type];
    // replace this event with event from end of list (may be same if this was last event)
    SD->events_by_type[old_event_ref] = SD->events_by_type[old_final_event_ref];
    // now negate final event
    SD->events_by_type[old_event_type*SD->nsites+old_final_event_ref] = -1;
    // and decrement number of events
    SD->events_per_type[old_event_type] = old_final_event_ref - 1;

// insert this site at end of new_event_type list
    int new_event_ref = SD->events_per_type[new_event_type]+1;
    SD->events_per_type[new_event_type] = new_event_ref;
    SD->events_by_type[new_event_type*SD->nsites+new_event_ref] = i;
    SD->event_types[i] = new_event_type;
    SD->event_refs[i] = new_event_ref;

}

int update_events_i( int move_site, struct SimData *SD){
    int i,j,state_i;
    //int n_possible_events = 0;
    //float event_rate = 0.0;
    int model_number = SD->model_number;
    int nneighbors_update_per_site = SD->nneighbors_update_per_site;

    //first do it for this site
//    int site_idx = SD->events_by_type[SD->event_types[event_i]*SD->nsites+SD->event_refs[event_i]]; // need this below
    i = move_site;
    state_i = SD->configuration[i];
    int old_event_type = SD->event_types[i];
    int new_event_type = get_event_type(i,SD);
    SD->events[i] = switch_state( state_i, model_number );
    //now change where it is in the list of events for that type
    change_event_type( i, old_event_type, new_event_type, SD );

    // now do it for neighbors this site affects
    for(j=0;j<nneighbors_update_per_site;j++){
        i = SD->neighbors_update[nneighbors_update_per_site*move_site+j];
        state_i = SD->configuration[i];
        old_event_type = SD->event_types[i];
        new_event_type = get_event_type(i,SD);
        SD->events[i] = switch_state( state_i, model_number );
        change_event_type( i, old_event_type, new_event_type, SD );
    }
    return 0;
}

int update_configuration( int event_i, struct SimData *SD){
    int j,k;
    int change_idx = SD->event_refs[event_i];
    int result = SD->events[change_idx];
    // Debug code:
    //printf("Changing spin %i to %i, which effects %i %i %i with rate %f\n",change_idx,result, SD->neighbors_update[SD->nneighbors_update_per_site*change_idx+0], SD->neighbors_update[SD->nneighbors_update_per_site*change_idx+1], SD->neighbors_update[SD->nneighbors_update_per_site*change_idx+2],SD->event_rates[change_idx]);
    SD->configuration[change_idx] = result;
    if(SD->model_number<MAXZEROONEMODEL) SD->persistence_array[change_idx] = 0;

    //also update dual representation
    if(SD->model_number==10){
        //change this site
        SD->dual_configuration[change_idx] = 1 - SD->dual_configuration[change_idx];
        //now find subset of affected neigbhors which have this spin in their plaquette
        for(j=0;j<SD->nneighbors_update_per_site;j++){
            int affected_neighbor_idx = SD->neighbors_update[SD->nneighbors_update_per_site*change_idx+j];
            for(k=0;k<SD->nneighbors_per_site;k++){
                int neighbors_neighbor_idx = SD->neighbors[SD->nneighbors_per_site*affected_neighbor_idx+k];
                if(neighbors_neighbor_idx==change_idx){
                    SD->dual_configuration[affected_neighbor_idx] = 
                        1 - SD->dual_configuration[affected_neighbor_idx];
                    break;
                }
            }
        }
    }
    //event_storage[SD->current_step] = change_idx;
    return 0;
}

double sum_rates( struct SimData *SD ){
    int i;
    double total_rate = 0.0;
    int n_possible_events = 0;
    if(SD->n_event_types>0){
        total_rate = SD->events_per_type[0]*SD->event_rates[0];
        SD->cumulative_rates[0] = total_rate;
        n_possible_events = SD->events_per_type[0];
    }

    for(i=1;i<SD->n_event_types;i++){
        float rate_type_i = SD->events_per_type[i]*SD->event_rates[i];
        total_rate = total_rate + rate_type_i;
        SD->cumulative_rates[i] = SD->cumulative_rates[i-1] + rate_type_i;
        n_possible_events += SD->events_per_type[i];
    }
    SD->total_rate = total_rate;
    SD->n_possible_events = n_possible_events;
    return total_rate;
} 

//binary search through cumulative event probabilities
int b_find_event( float searchval, struct SimData *SD){
    int final_event = SD->n_event_types - 1;
    int idx0 = 0;
    int idx1 = final_event;
    if ( SD->n_event_types < 1 ){
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

int run_kmc_spin_steps(int nsteps,struct SimData *SD){
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

        //int j;
        //printf("prob*total_rate: %f\n",prob*total_rate);
        //for(j=0;j<SD->n_possible_events;j++){ 
        //    printf("%f ",SD->cumulative_rates[j]);
        //}
        int event_i = b_find_event( prob*total_rate, SD);
        //printf("\nevent_i: %i\n",event_i);
        update_configuration( event_i, SD );
        n_possible_events = update_events_i( event_i, SD );
        dt = -log(prob)/total_rate;
        t += dt;
    }
    SD->time = t;
    return return_val;
}
int copy_configuration_prev(struct SimData *SD){
    int i;
    for(i=0;i<SD->nsites;i++){ 
        SD->prev_configuration[i] = SD->configuration[i];
    }
    return 0;
}

int run_kmc_spin(float stop_time,struct SimData *SD){
    int step = 0;
    int n_possible_events = SD->n_possible_events;
    int return_val = 0;
    float t = SD->time;
    float dt = 0;
    while(t<stop_time){
        if(n_possible_events<1){
            return_val = -1;
            break;
        }

        float prob = get_frandom();
        float total_rate = sum_rates(SD);
        dt = -log(prob)/total_rate;
        t += dt;
        if(t >= stop_time ){
            copy_configuration_prev(SD);
        }

        int event_type_i = b_find_event( prob*total_rate, SD);
        int rand_event = get_irandom( 0, SD->events_per_type[event_type_i] );
        int move_site = SD->events_by_type[event_type_i*SD->nsites+rand_event];
        update_configuration( move_site, SD );
        update_events_i( move_site, SD );
        step++;
    }
    SD->time = t;
    return return_val;
}

int setup_spin_system(struct SimData *SD){
    set_seed(SD->seed);
    all_events(SD);
    sum_rates(SD);
    return 0;
}

int cleanup_spin_system(struct SimData *SD){
    return 0;
}

