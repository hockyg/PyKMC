# Fredrickson Andersen (FA) Model 
# Original Source:
#     G. H. Fredrickson and H. C. Andersen
#     Kinetic Ising Model of the Glass Transition. PRL 53, 1244 (1984).
import numpy as np
cimport numpy as np
import ctypes as ct
changes_per_step = 1
model_name = "FA"

class CubeClass(object):
    def __init__(self, int side_length ):
        self.side_length = side_length
        self.nsites = side_length*side_length*side_length

    def Neighbors(self):
        """ Calculates all neighbors for all sites """
        cdef int nsites = self.nsites
        cdef int nneighbors_per_site = 6
        cdef int site_idx, j
        cdef np.ndarray[np.int_t,ndim=2] neighbors = np.zeros((nsites,nneighbors_per_site),dtype=ct.c_int)
        for site_idx in range(nsites):
            neighbors_i = self.NeighborsI(site_idx)
            for j in range(nneighbors_per_site):
                neighbors[site_idx,j] = neighbors_i[j]
        return nneighbors_per_site, neighbors
    
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
        cdef int n_idx = (index0+1)%(self.side_length) #north
        cdef int s_idx = (index0-1)%(self.side_length) #south

        return self.row_col_to_idx(n_idx,index1,index2), self.row_col_to_idx(s_idx,index1,index2), self.row_col_to_idx(index0,f_idx,index2), self.row_col_to_idx(index0,b_idx,index2), self.row_col_to_idx(index0,index1,e_idx), self.row_col_to_idx(index0,index1,w_idx)

    def row_col_to_idx(self, int index0, int index1, int index2):
        return self.side_length*self.side_length*index0+index1*self.side_length+index2

class SquareClass(object):
    def __init__(self, int side_length ):
        self.side_length = side_length
        self.nsites = side_length*side_length

    def Neighbors(self):
        """ Calculates all neighbors for all sites """
        cdef int nsites = self.nsites
        cdef int nneighbors_per_site = 4
        cdef int site_idx, j
        cdef np.ndarray[np.int_t,ndim=2] neighbors = np.zeros((nsites,nneighbors_per_site),dtype=ct.c_int)
        for site_idx in range(nsites):
            neighbors_i = self.NeighborsI(site_idx)
            for j in range(nneighbors_per_site):
                neighbors[site_idx,j] = neighbors_i[j]
        return nneighbors_per_site, neighbors
    
    def NeighborsI( self, int site_idx ):
        """ Returns the neighbors of lattice site site_idx """
        # first find x and y position
        cdef int col_num = site_idx%self.side_length
        cdef int row_num = (site_idx-col_num)/self.side_length
        
        cdef int f_column = (col_num+1)%self.side_length
        cdef int b_column = (col_num-1)%self.side_length
        cdef int u_row = (row_num+1)%(self.side_length)
        cdef int d_row = (row_num-1)%(self.side_length)

        return self.row_col_to_idx(row_num,f_column), self.row_col_to_idx(row_num,b_column), self.row_col_to_idx(u_row,col_num), self.row_col_to_idx(d_row,col_num)

    def row_col_to_idx(self, int row, int col):
        return row*self.side_length+col

class LinearClass(object):
    def __init__(self, int side_length ):
        self.nsites = side_length

    def Neighbors(self):
        """ Calculates all neighbors for all sites """
        # note that there are exactly 2*nsites neighbors
        cdef int nsites = self.nsites
        cdef int nneighbors_per_site = 2
        cdef int site_idx, j
        cdef np.ndarray[np.int_t,ndim=2] neighbors = np.zeros((nsites,nneighbors_per_site),dtype=ct.c_int)
        for site_idx in range(nsites):
            neighbors_i = self.NeighborsI(site_idx, nsites)
            for j in range(nneighbors_per_site):
                neighbors[site_idx,j] = neighbors_i[j]
        return nneighbors_per_site, neighbors
    
    def NeighborsI( self, int site_idx, int nsites ):
        """ Returns the neighbors of lattice site site_idx """
        # forward neighbor
        cdef int neighborf = site_idx + 1
        # note, in the next line, would want floor( neighborf/nsites ) but since these are ints, uses auto flooring property of python.
        # seems to work and be over 10x faster, but potential problem area later
        neighborf = neighborf - nsites*( neighborf/ nsites )
        # backward neighbor
        cdef int neighborb = site_idx - 1
        neighborb = neighborb - nsites*( neighborb/ nsites )
    
        return neighborb, neighborf

#    def PeriodicDistance( self, site_idx1, site_idx2, nsites ):
#        dx = 1.*( site_idx2 - site_idx1 )
#        return dx - nsites*np.floor( dx/nsites + 0.5 )

LatticeRegistry = {"linear":LinearClass, "square":SquareClass, "cube":CubeClass }

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
