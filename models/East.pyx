# East model
# Original Source:
#
import numpy as np
cimport numpy as np
import ctypes as ct
changes_per_step = 1
model_name = "East"

class CubeClass(object):
    def __init__(self, int side_length ):
        self.side_length = side_length
        self.nsites = side_length*side_length*side_length
        self.n_event_types = 8

    def EventRates(self,double temp,int dynamics_number):
        if ( dynamics_number != 0 ): raise AssertionError("East model only supports metropolis dynamics")
        cdef int i,j
        cdef np.ndarray event_rates = np.zeros(self.n_event_types,dtype=ct.c_double)
        cdef np.ndarray constraints = np.array( np.arange(0.,4.,1), dtype=ct.c_double) # 0,1,2,3 # number of affecting neighbors
        cdef np.ndarray spin_states = np.array( np.arange(0.,2.,1), dtype=ct.c_double) # 0 and 1
        cdef float excitations_created
        for i in range(4):
            constraint = constraints[i]
            for j in range(2):
                spin_state = spin_states[j]
                event_rates[2*constraint+spin_state] = (1-spin_state)*constraint*np.exp(-1/temp) + spin_state*constraint
        return event_rates

    def Neighbors(self):
        """ Calculates all neighbors for all sites """
        cdef int nsites = self.nsites
        cdef int nneighbors_per_site = 3
        cdef int site_idx, j
        cdef np.ndarray[np.int32_t,ndim=2] neighbors = np.zeros((nsites,nneighbors_per_site),dtype=ct.c_int)
        cdef np.ndarray[np.int32_t,ndim=2] neighbors_update = np.zeros((nsites,nneighbors_per_site),dtype=ct.c_int)
        for site_idx in range(nsites):
            neighbors_i, neighbors_i_update = self.NeighborsI(site_idx)
            for j in range(nneighbors_per_site):
                neighbors[site_idx,j] = neighbors_i[j]
            for j in range(nneighbors_per_site):
                neighbors_update[site_idx,j] = neighbors_i_update[j]
        return nneighbors_per_site, nneighbors_per_site, neighbors, neighbors_update
    
    def NeighborsI( self, int site_idx ):
        """ Returns the neighbors of lattice site site_idx """
        # first find x and y position
        cdef int index2 = site_idx%self.side_length
        cdef int index1 = ((site_idx-index2)/self.side_length)%self.side_length
        cdef int index0 = (site_idx-index1*self.side_length-index2)/self.side_length/self.side_length
        
        cdef int e_idx = (index2+1)%self.side_length #east
        cdef int w_idx = (index2-1)%self.side_length #west
        cdef int f_idx = (index1+1)%self.side_length #front
        cdef int b_idx = (index1-1)%self.side_length #back
        cdef int n_idx = (index0-1)%(self.side_length) #north
        cdef int s_idx = (index0+1)%(self.side_length) #south

        return [ self.row_col_to_idx(s_idx,index1,index2), self.row_col_to_idx(index0,b_idx,index2), self.row_col_to_idx(index0,index1,w_idx) ], [ self.row_col_to_idx(n_idx,index1,index2), self.row_col_to_idx(index0,f_idx,index2), self.row_col_to_idx(index0,index1,e_idx) ]

    def row_col_to_idx(self, int index0, int index1, int index2):
        return self.side_length*self.side_length*index0+index1*self.side_length+index2

    def RandomConfiguration(self, double temperature):
        return RandomConfiguration(self.nsites, temperature)

class SquareClass(object):
    def __init__(self, int side_length ):
        self.side_length = side_length
        self.nsites = side_length*side_length
        self.n_event_types = 6

    def EventRates(self,double temp, int dynamics_number):
        if ( dynamics_number != 0 ): raise AssertionError("East model only supports metropolis dynamics")
        cdef int i,j
        cdef np.ndarray event_rates = np.zeros(self.n_event_types,dtype=ct.c_double)
        cdef np.ndarray constraints = np.array( np.arange(0.,3.,1), dtype=ct.c_double) # 0,1,2 # number of affecting neighbors
        cdef np.ndarray spin_states = np.array( np.arange(0.,2.,1), dtype=ct.c_double) # 0 and 1
        cdef float excitations_created
        for i in range(3):
            constraint = constraints[i]
            for j in range(2):
                spin_state = spin_states[j]
                event_rates[2*constraint+spin_state] = (1-spin_state)*constraint*np.exp(-1/temp) + spin_state*constraint
        return event_rates

    def Neighbors(self):
        """ Calculates all neighbors for all sites """
        cdef int nsites = self.nsites
        cdef int nneighbors_per_site = 2
        cdef int site_idx, j
        cdef np.ndarray[np.int32_t,ndim=2] neighbors = np.zeros((nsites,nneighbors_per_site),dtype=ct.c_int)
        cdef np.ndarray[np.int32_t,ndim=2] neighbors_update = np.zeros((nsites,nneighbors_per_site),dtype=ct.c_int)
        for site_idx in range(nsites):
            neighbors_i, neighbors_i_update = self.NeighborsI(site_idx)
            for j in range(nneighbors_per_site):
                neighbors[site_idx,j] = neighbors_i[j]
            for j in range(nneighbors_per_site):
                neighbors_update[site_idx,j] = neighbors_i_update[j]
        return nneighbors_per_site, nneighbors_per_site, neighbors, neighbors_update
   
    def NeighborsI( self, int site_idx ):
        """ Returns the neighbors of lattice site site_idx """
        # first find x and y position
        cdef int col_num = site_idx%self.side_length
        cdef int row_num = (site_idx-col_num)/self.side_length
        
        cdef int f_column = (col_num+1)%self.side_length
        cdef int b_column = (col_num-1)%self.side_length
        cdef int u_row = (row_num-1)%(self.side_length)
        cdef int d_row = (row_num+1)%(self.side_length)

        return [ self.row_col_to_idx(row_num,b_column), self.row_col_to_idx(d_row,col_num) ], [ self.row_col_to_idx(row_num,f_column), self.row_col_to_idx(u_row,col_num) ]

    def row_col_to_idx(self, int row, int col):
        return row*self.side_length+col

    def RandomConfiguration(self, double temperature):
        return RandomConfiguration(self.nsites, temperature)

class LinearClass(object):
    def __init__(self, int side_length ):
        self.nsites = side_length
        self.side_length = side_length
        self.n_event_types = 4

    def EventRates(self,double temp,int dynamics_number):
        if ( dynamics_number != 0 ): raise AssertionError("East model only supports metropolis dynamics")
        cdef int i,j
        cdef np.ndarray event_rates = np.zeros(self.n_event_types,dtype=ct.c_double)
        cdef np.ndarray constraints = np.array( np.arange(0.,2.,1), dtype=ct.c_double) # 0 and 1
        cdef np.ndarray spin_states = np.array( np.arange(0.,2.,1), dtype=ct.c_double) # 0 and 1
        cdef float excitations_created
        for i in range(2):
            constraint = constraints[i]
            for j in range(2):
                spin_state = spin_states[j]
                event_rates[2*constraint+spin_state] = (1-spin_state)*constraint*np.exp(-1/temp) + spin_state*constraint
        return event_rates

    def Neighbors(self):
        """ Calculates all neighbors for all sites """
        cdef int nsites = self.nsites
        cdef int nneighbors_per_site = 1
        cdef int site_idx, j
        cdef np.ndarray[np.int32_t,ndim=2] neighbors = np.zeros((nsites,nneighbors_per_site),dtype=ct.c_int)
        cdef np.ndarray[np.int32_t,ndim=2] neighbors_update = np.zeros((nsites,nneighbors_per_site),dtype=ct.c_int)
        for site_idx in range(nsites):
            neighbors_i, neighbors_i_update = self.NeighborsI(site_idx)
            for j in range(nneighbors_per_site):
                neighbors[site_idx,j] = neighbors_i[j]
            for j in range(nneighbors_per_site):
                neighbors_update[site_idx,j] = neighbors_i_update[j]
        return nneighbors_per_site, nneighbors_per_site, neighbors, neighbors_update
 
    def NeighborsI( self, int site_idx):
        """ Returns the neighbors of lattice site site_idx """
        # forward neighbor
        cdef int neighborf = site_idx + 1
        # note, in the next line, would want floor( neighborf/self.side_length ) but since these are ints, uses auto flooring property of python.
        # seems to work and be over 10x faster, but potential problem area later
        neighborf = neighborf - self.side_length*( neighborf/ self.side_length )
        # backward neighbor
        cdef int neighborb = site_idx - 1
        neighborb = neighborb - self.side_length*( neighborb/ self.side_length )
    
        return [neighborb], [neighborf] # second is update

    def RandomConfiguration(self, double temperature):
        return RandomConfiguration(self.nsites, temperature)

LatticeRegistry = {"linear":LinearClass, "square":SquareClass, "cube":CubeClass }

def InitializeArrays( int nsites, int n_event_types, int nactive):
    cdef np.ndarray events = np.zeros(nsites,dtype=ct.c_int)
    cdef np.ndarray event_types = np.zeros(nsites,dtype=ct.c_int)
    cdef np.ndarray events_by_type = np.zeros((n_event_types,nsites),dtype=ct.c_int)
    cdef np.ndarray events_per_type = np.zeros(n_event_types,dtype=ct.c_int)
    cdef np.ndarray event_refs = -1*np.ones(nsites,dtype=ct.c_int)
    cdef np.ndarray event_rates = np.zeros(n_event_types,dtype=ct.c_double)
#    cdef np.ndarray event_ref_rates = np.zeros(nsites,dtype=ct.c_double)
    cdef np.ndarray cumulative_rates = np.zeros(n_event_types,dtype=ct.c_double)
    cdef np.ndarray persistence_array= np.ones(nsites,dtype=ct.c_int)
    # note: this will set activelist as 0-(nactive-1) which is only correct for nactive=nsites
    #         hence this will be written over in the initialization step
    cdef np.ndarray activelist = np.arange(nactive,dtype=ct.c_int)
    # all sites will be set to active by default
    cdef np.ndarray isactivelist = np.ones(nsites,dtype=ct.c_int)
    return {"events": events,
            "event_types": event_types,
            "events_by_type": events_by_type,
            "events_per_type": events_per_type,
            "event_refs": event_refs,
            "event_rates": event_rates, 
            "cumulative_rates": cumulative_rates,
            "persistence_array": persistence_array,
            "activelist": activelist,
            "isactivelist": isactivelist,
           }

def RandomConfiguration( int nsites, double temperature ):
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
