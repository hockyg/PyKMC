import numpy as np
cimport numpy as np
from math import ceil

def BSearchProb( float prob, int nevents, 
                  np.ndarray[np.float_t,ndim=1] cumulative_event_probs ):
    cdef int final_event = nevents-1
    cdef int idx0 = 0 
    cdef int idx1 = final_event
    cdef int halfidx
    cdef float half_prob
    if cumulative_event_probs[0] > prob:
        return 0
    if cumulative_event_probs[final_event] < prob:
        return final_event

    while (idx1-idx0)>1:
        halfidx = (idx1+idx0)/2
        half_prob = cumulative_event_probs[halfidx]
        if half_prob > prob:
            idx1 = halfidx
        elif half_prob < prob:
            idx0 = halfidx
        else:
            return halfidx 

    return idx1
