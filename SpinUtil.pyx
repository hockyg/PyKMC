import numpy as np
cimport numpy as np
from math import ceil

def persistence( int nsites, int ndownspins_start, np.ndarray[np.int_t,ndim=1] persistence_array ):
    cdef int i
    cdef downpersist = 0
    cdef totalpersist = 0
    for i in range(nsites):
        if persistence_array[i]<0:
            downpersist = downpersist+1
            totalpersist = totalpersist+1
        elif persistence_array[i]>0:
            totalpersist = totalpersist+1
    return float(totalpersist)/float(nsites), float(downpersist)/float(ndownspins_start)
    

def textline_box( lineoftext ):
    characters = len( lineoftext )
    header = "+" + "-"*characters + "+"
    center = "|"+lineoftext.strip()+"|"
    footer = "+" + "-"*characters + "+"
    return "\n".join( (header,center,footer) )

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
    if nevents < 1:
        return -1

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
