import numpy as np
cimport numpy as np
from math import ceil
import sys

def tee_logfile( output_prefix ):
    from SpinObj import tee
    logfile = output_prefix+'.log'
    fh = open(logfile,'w')
    stdout_sv = sys.stdout
    stderr_sv = sys.stderr
    sys.stdout = tee( stdout_sv, fh )
    sys.stderr = tee( stderr_sv, fh )

def c_to_T_ideal( int nsites, np.ndarray[np.int32_t,ndim=1] configuration ):
    """converts the concentration of up spins for a lattice gas of zeros and ones to an effective temperature"""
    cdef int i
    cdef float c = 0.0
    for i in range(nsites):
        c = c + configuration[i]
    c = c/nsites
    # c = 1/(1+exp(1/T))
    return 1./np.log(1/c-1)

def persistence( int nsites, int ndownspins_start, np.ndarray[np.int32_t,ndim=1] persistence_array ):
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
