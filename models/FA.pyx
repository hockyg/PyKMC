# Fredrickson Andersen (FA) Model 
# Original Source:
#     G. H. Fredrickson and H. C. Andersen
#     Kinetic Ising Model of the Glass Transition. PRL 53, 1244 (1984).
import numpy as np
cimport numpy as np
changes_per_step = 1

cdef ConstraintI( int site_idx, np.ndarray[np.int_t,ndim=1] configuration, int nsites, np.ndarray[np.int_t,ndim=2] neighbors, int nneighbors_per_site ):
    """ This function gives the FA model 'constraint' 

    Functions adjacent to activations are activated.
    Functions bounded by n active sites are n-fold activated.
    """
    cdef float constraint = 0.0
    cdef int j
    for j in range(nneighbors_per_site):
        constraint = constraint+configuration[neighbors[site_idx,j]]
    return constraint

def AllEvents( float betaexp, np.ndarray[np.int_t,ndim=2] events, np.ndarray[np.float_t,ndim=1] event_rates, np.ndarray[np.int_t,ndim=1] configuration, int nsites, np.ndarray[np.int_t,ndim=2] neighbors, int nneighbors_per_site ):
    """ This function generates the full list of allowed FA moves and their respective rates """
    cdef int i, state_i, n_possible_events
    cdef int j = 0
    cdef float constraint
    cdef float event_rate = 0.0
    n_possible_events = 0

    # first clear out rates, and events, just in case
    for j in range(n_possible_events):
        events[n_possible_events,0] = 0
        events[n_possible_events,1] = 0
        event_rates[j] = 0.0

    for i in range(nsites):
        state_i = configuration[i]
        constraint = ConstraintI( i, configuration, nsites, neighbors, nneighbors_per_site )
        # next line means that if up, flip down and if down, flip up
        # next line is equivalent to
        # rate_{0->1} = e^(-beta) * constraint
        # rate_{1->0} = constraint
        event_rate = (1-state_i)*constraint*betaexp + state_i*constraint
        if( event_rate > 0 ):
            events[n_possible_events,0] = i
            events[n_possible_events,1] = (1-state_i) 
            event_rates[n_possible_events] = event_rate
            n_possible_events += 1
    return n_possible_events

def UpdateConfiguration( np.ndarray[np.int_t,ndim=1] configuration, np.ndarray[np.int_t,ndim=2] events, int event_i ):
    cdef int change_idx = events[event_i,0]
    cdef int result = events[event_i,1]
    configuration[change_idx] = result
