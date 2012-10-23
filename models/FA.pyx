# Fredrickson Andersen (FA) Model 
# Original Source:
#     G. H. Fredrickson and H. C. Andersen
#     Kinetic Ising Model of the Glass Transition. PRL 53, 1244 (1984).
import numpy as np
cimport numpy as np
import ctypes as ct
changes_per_step = 1
model_name = "FA"

def InitializeArrays( int nsites, int max_steps ):
    cdef np.ndarray events = np.zeros(nsites,dtype=ct.c_int)
    cdef np.ndarray event_rates = np.zeros(nsites,dtype=ct.c_float)
    cdef np.ndarray event_refs = -1*np.ones(nsites,dtype=ct.c_int)
    cdef np.ndarray event_ref_rates = np.zeros(nsites,dtype=ct.c_float)
    cdef np.ndarray cumulative_rates = np.zeros(nsites,dtype=ct.c_float)
    cdef np.ndarray persistence_array= np.ones(nsites,dtype=ct.c_int)
    return {"events": events,
            "event_rates": event_rates, 
            "event_refs": event_refs,
            "event_ref_rates": event_ref_rates,
            "cumulative_rates": cumulative_rates,
            "persistence_array": persistence_array,
           }

def RandomConfiguration( int nsites, float temperature ):
    """ Generates a random configuration commensurate with the temperature
       
        The average excitation value of a site in a non-interacting lattice gas with H = \sum_i n_i
        is \langle n_i \rangle = 1/(\exp(\beta)+1)

    """
    cdef int i
    cdef np.ndarray configuration = np.zeros(nsites,dtype=ct.c_int)
    cdef np.ndarray rand_values = np.random.random(size=nsites)
    cdef double beta = 1/temperature
    cdef double avg_site_value = 1./(np.exp(beta)+1.0)

    for i in range(nsites):
        if rand_values[i] < avg_site_value:
            configuration[i] = 1
    return configuration
