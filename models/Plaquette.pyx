# Plaquette models
# Square Plaquette model:
#    See e.g. Robert L. Jack, Ludovic Berthier, Juan P. Garrahan.
#             "Static and dynamic length scales in a simple glassy plaquette model". 
#             Phys Rev E 72, 016103 (2005).
#    Hamiltonian will actually be H = -1/2 \sum_{ij}( sig_{i,j} sig_{i,j+1} sig_{i+1,j} sig_{i+1,j+1} - 1 )
#       Such that p_{ij} = sig_{i,j} sig_{i,j+1} sig_{i+1,j} sig_{i+1,j+1}
#       With a dual hamiltonian of H = \sum_{ij} ( 1 - p_{ij} )/2 = \sum_{ij} n_{ij}
#       And p_ij = -1 give \n_ij = 1
#       And p_ij = +1 give \n_ij = 0
import numpy as np
cimport numpy as np
import ctypes as ct
changes_per_step = 1
model_name = "Plaquette"

cdef gen_pascal_parity( np.ndarray[np.int32_t,ndim=2] pascal_parity, int linear_size ):
    cdef int i,j
    for i in range(1,linear_size):
        for j in range(0,linear_size-2*i):
            pascal_parity[i,i+j] = (pascal_parity[i-1,i+j]+pascal_parity[i,i+j-1])%2
            pascal_parity[i+j,i] = (pascal_parity[i+j-1,i]+pascal_parity[i+j,i-1])%2

def test_triangle_dual( configuration, dual_configuration, int linear_size ):
    cdef int i,j
#    cdef np.ndarray dual_configuration = np.zeros(configuration.shape,dtype=ct.c_int)
    cdef int sumvar = 0
    for i in range(linear_size):
        for j in range(linear_size):
            spinij = configuration[i,j]
            spinip1j = configuration[(i+1)%linear_size,j]
            spinip1jp1 = configuration[(i+1)%linear_size,(j+1)%linear_size]
            sumvar = sumvar + abs(dual_configuration[i,j] - ( 1-spinij*spinip1j*spinip1jp1 )/2)
    return sumvar

cdef pascal_parity_f(int row, int col):
    """get the parity of a pascal triangle entry using lucas's theorem.
       a corrolary of this theorem is that (m,n) is divisble by 2 
         iff 1 digit in base 2 of n is greater than that digit of m"""
    #note, row should always be larger than col
    cdef int digit_m, digit_n
    cdef int tmp_m = row
    cdef int tmp_n = col

    # next while loop generates reverse binary representation of number
    while( tmp_m >1 and tmp_n > 1):
        digit_m = tmp_m%2
        digit_n = tmp_n%2
        if (digit_n > digit_m ):
            return 0
        tmp_m = tmp_m//2
        tmp_n = tmp_n//2

    if ( tmp_n%2 > tmp_m%2 ):
        return 0

    return 1

def pascal_parity_f_slow(row, col):
    """get the parity of a pascal triangle entry using lucas's theorem.
       a corrolary of this theorem is that (m,n) is divisble by 2 
         iff 1 digit in base 2 of n is greater than that digit of m"""
    #note, row should always be larger than col
    row2 = np.binary_repr(row)
    col2 = np.binary_repr(col,len(row2))
    for i in range(len(row2)):
        if col2[i]>row2[i]:
            return 0
    return 1

class TriangleClass(object):
    def __init__(self, int side_length ):
        self.side_length = side_length
        self.nsites = side_length*side_length
        self.n_event_types = 4

    def EventRates(self,double temp,int dynamics_number):
        cdef int i
        cdef np.ndarray event_rates = np.zeros(self.n_event_types,dtype=ct.c_double)
        # can create or distroy one excitation with a single spin flip
        cdef np.ndarray excitations_created_array = np.array( np.arange(-3.,4.,2), dtype=ct.c_double)
        cdef float excitations_created
        for i in range(self.n_event_types):
            excitations_created = excitations_created_array[i];
            # Glauber dynamics
            if dynamics_number == 1:
                event_rates[i] = 1./(1+np.exp(excitations_created/temp))
            # Metropolis
            else:
                event_rates[i] = (excitations_created>0)*np.exp(-excitations_created/temp) + (excitations_created<=0)
        return event_rates

    def Neighbors(self):
        """ Calculates all neighbors for all sites """
        cdef int nsites = self.nsites
        cdef int nneighbors_per_site = 2
        cdef int nneighbors_update_per_site = 2
        cdef int site_idx, j
        cdef np.ndarray[np.int32_t,ndim=2] neighbors = np.zeros((nsites,nneighbors_per_site),dtype=ct.c_int)
        cdef np.ndarray[np.int32_t,ndim=2] neighbors_update = np.zeros((nsites,nneighbors_update_per_site),dtype=ct.c_int)
        for site_idx in range(nsites):
            neighbors_i, neighbors_i_update = self.NeighborsI(site_idx)
            for j in range(nneighbors_per_site):
                neighbors[site_idx,j] = neighbors_i[j]
            for j in range(nneighbors_update_per_site):
                neighbors_update[site_idx,j] = neighbors_i_update[j]
        return nneighbors_per_site, nneighbors_update_per_site, neighbors, neighbors_update
   
    def NeighborsI( self, int site_idx ):
        """ Returns the neighbors of lattice site site_idx """
        # first find x and y position
        cdef int col_num = site_idx%self.side_length
        cdef int row_num = (site_idx-col_num)/self.side_length
        
        cdef int f_column = (col_num+1)%self.side_length
        cdef int b_column = (col_num-1)%self.side_length
        cdef int u_row = (row_num-1)%(self.side_length)
        cdef int d_row = (row_num+1)%(self.side_length)

        # first neighbors, then neighbors to update
        # first 3 of neighbor update have this site in their plaquette
        return [ 
                 self.row_col_to_idx(d_row,col_num),
                 self.row_col_to_idx(d_row,f_column) ], \
               [ 
                 self.row_col_to_idx(u_row,col_num), 
                 self.row_col_to_idx(u_row,b_column), ],
#                 self.row_col_to_idx(d_row,col_num), 
#                 self.row_col_to_idx(d_row,f_column),  ],

    def row_col_to_idx(self, int row, int col):
        return row*self.side_length+col

    def introduce_defect( self, row, col, configuration ):
        cdef int linear_size = self.side_length
        cdef int i,j
        # note, one 2x speedup could come from realizing that pascal parity is 2 fold symmetric
        #    and filling in configuration from two sides
        for i in range(linear_size):
            for j in range(i+1):
                # i.e., if i Choose j is odd, flip spin
                if( pascal_parity_f(i,j) > 0 ):
                    configuration[(row-i)%linear_size,(col-j)%linear_size]*=-1

    def RandomConfiguration( self, double temperature ):
        """ Generates a random configuration commensurate with the temperature
           
        """
        import sys
        cdef int i,j,k,l
        cdef int d_idx, r_idx, dr_value, pij
        cdef int side_length = self.side_length
        cdef int half_side_length = side_length/2

        cdef float size_log2 = np.log(self.side_length)/np.log(2)
        if  (size_log2 - np.floor(size_log2))>0:
            print "Side length must be a power of 2 to generate a random triangular plaquette configuration"
            sys.exit(1)

        cdef np.ndarray[np.int32_t,ndim=2] dual_configuration = RandomConfigurationIdeal( self.nsites, temperature ).reshape((side_length,side_length))
        # note, start all spin up
        cdef np.ndarray[np.int32_t,ndim=2] configuration = np.ones((side_length,side_length),dtype=ct.c_int)

        for i in range(side_length):
            for j in range(side_length):
                if dual_configuration[i,j]>0: self.introduce_defect( i, j, configuration )

        #tmp_dual = test_dual(configuration,side_length)
        #print ( dual_configuration-tmp_dual ).sum()
        return configuration.flatten(), dual_configuration.flatten()


class SquareClass(object):
    def __init__(self, int side_length ):
        self.side_length = side_length
        self.nsites = side_length*side_length
        self.n_event_types = 5

    def EventRates(self,double temp, int dynamics_number):
        cdef int i
        cdef np.ndarray event_rates = np.zeros(self.n_event_types,dtype=ct.c_double)
        cdef np.ndarray excitations_created_array = np.array( np.arange(-4.,5.,2), dtype=ct.c_double)
        cdef float excitations_created
        for i in range(self.n_event_types):
            excitations_created = excitations_created_array[i];
            # Glauber dynamics
            if dynamics_number == 1:
                event_rates[i] = 1./(1+np.exp(excitations_created/temp))
            # Metropolis
            else:
                event_rates[i] = (excitations_created>0)*np.exp(-excitations_created/temp) + (excitations_created<=0)
        return event_rates

    def Neighbors(self):
        """ Calculates all neighbors for all sites """
        cdef int nsites = self.nsites
        cdef int nneighbors_per_site = 3
        cdef int nneighbors_update_per_site = 3
        cdef int site_idx, j
        cdef np.ndarray[np.int32_t,ndim=2] neighbors = np.zeros((nsites,nneighbors_per_site),dtype=ct.c_int)
        cdef np.ndarray[np.int32_t,ndim=2] neighbors_update = np.zeros((nsites,nneighbors_update_per_site),dtype=ct.c_int)
        for site_idx in range(nsites):
            neighbors_i, neighbors_i_update = self.NeighborsI(site_idx)
            for j in range(nneighbors_per_site):
                neighbors[site_idx,j] = neighbors_i[j]
            for j in range(nneighbors_update_per_site):
                neighbors_update[site_idx,j] = neighbors_i_update[j]
        return nneighbors_per_site, nneighbors_update_per_site, neighbors, neighbors_update
   
    def NeighborsI( self, int site_idx ):
        """ Returns the neighbors of lattice site site_idx """
        # first find x and y position
        cdef int col_num = site_idx%self.side_length
        cdef int row_num = (site_idx-col_num)/self.side_length
        
        cdef int f_column = (col_num+1)%self.side_length
        cdef int b_column = (col_num-1)%self.side_length
        cdef int u_row = (row_num-1)%(self.side_length)
        cdef int d_row = (row_num+1)%(self.side_length)

        # first neighbors, then neighbors to update
        # first 3 of neighbor update have this site in their plaquette
        return [ 
                 self.row_col_to_idx(row_num,f_column),
                 self.row_col_to_idx(d_row,col_num),
                 self.row_col_to_idx(d_row,f_column) ], \
               [ 
                 self.row_col_to_idx(row_num,b_column), 
                 self.row_col_to_idx(u_row,col_num), 
                 self.row_col_to_idx(u_row,b_column), ]
#                 self.row_col_to_idx(row_num,f_column), 
#                 self.row_col_to_idx(u_row,f_column), 
#                 self.row_col_to_idx(d_row,col_num), 
#                 self.row_col_to_idx(d_row,b_column), 
#                 self.row_col_to_idx(d_row,f_column),  ],

    def row_col_to_idx(self, int row, int col):
        return row*self.side_length+col

    def RandomConfiguration( self, double temperature ):
        """ Generates a random configuration (almost) commensurate with the temperature
           
        """
        import sys
        cdef int i,j,k,l
        cdef int d_idx, r_idx, dr_value, pij
        cdef int side_length = self.side_length
        cdef int half_side_length = side_length/2
        if not (side_length%4==0):
            print "Side length must be divisible by 4 to generate a random square plaquette configuration"
            sys.exit(1)
        dual_sub_configuration_flat = RandomConfigurationIdeal( self.nsites/4, temperature )

        cdef np.ndarray[np.int32_t,ndim=2] dual_configuration = np.zeros((side_length,side_length),dtype=ct.c_int)
        cdef np.ndarray[np.int32_t,ndim=2] dual_sub_configuration = dual_sub_configuration_flat.reshape((half_side_length,half_side_length))
 
        for k in range(2):
            for l in range(2):
                for i in range(half_side_length):
                    for j in range(half_side_length):
                        dual_configuration[k*half_side_length+i,l*half_side_length+j] = dual_sub_configuration[i,j]
                        
        cdef np.ndarray[np.int32_t,ndim=2] configuration = np.zeros((side_length,side_length),dtype=ct.c_int)
        cdef np.ndarray[np.int32_t,ndim=1] rand_spins

        # first fill in top row
        rand_spins = np.array( 2*np.random.randint(2,size=side_length)-1, dtype=ct.c_int )
        for i in range(side_length):
            configuration[0,i] = rand_spins[i]
        # then fill in first column
        rand_spins = np.array( 2*np.random.randint(2,size=side_length)-1, dtype=ct.c_int )
        for i in range(side_length):
            configuration[i,0] = rand_spins[i]

        # now fill in rest
        for i in range(side_length-1):
            for j in range(side_length-1):
                # note, plaquettes are defined as down and right, and so we will look down and right and fill in the down, right square
                d_idx = (i+1)%side_length
                r_idx = (j+1)%side_length
                # dual_configuration[i,j] = n_ij = ( 1 - p_{ij} )/2 
                # p_{ij} = 1- 2*n_ij
                pij = 1 - 2*dual_configuration[i,j] 
                dr_value = pij/configuration[i,j]/configuration[d_idx,j]/configuration[i,r_idx]
                configuration[d_idx,r_idx] = dr_value

        return configuration.flatten(), dual_configuration.flatten()

           
LatticeRegistry = {"square":SquareClass, 
                   "triangle":TriangleClass}

def PlaquetteEnergy(np.ndarray[np.int32_t,ndim=1] configuration,np.ndarray[np.int32_t,ndim=2] neighbors, int nsites, int nneighbors_per_site):
    cdef int i,j
    cdef int spin_prod
    cdef double site_e, total_e
    cdef int neighbor_idx
    total_e = 0.0
    for i in range(nsites):
        spin_prod = configuration[i]
        for j in range(nneighbors_per_site):
            neighbor_idx = neighbors[i,j] 
            spin_prod = spin_prod*configuration[neighbor_idx]
        site_e = -(spin_prod-1)/2
        total_e = total_e + site_e
    return total_e
 
def InitializeArrays( int nsites, int n_event_types ):
    cdef np.ndarray events = np.zeros(nsites,dtype=ct.c_int)
    cdef np.ndarray event_types = np.zeros(nsites,dtype=ct.c_int)
    cdef np.ndarray events_by_type = np.zeros((n_event_types,nsites),dtype=ct.c_int)
    cdef np.ndarray events_per_type = np.zeros(n_event_types,dtype=ct.c_int)
    cdef np.ndarray event_refs = -1*np.ones(nsites,dtype=ct.c_int)
    cdef np.ndarray event_rates = np.zeros(n_event_types,dtype=ct.c_double)
#    cdef np.ndarray event_ref_rates = np.zeros(nsites,dtype=ct.c_double)
    cdef np.ndarray cumulative_rates = np.zeros(n_event_types,dtype=ct.c_double)
    cdef np.ndarray persistence_array= np.ones(nsites,dtype=ct.c_int)
    return {"events": events,
            "event_types": event_types,
            "events_by_type": events_by_type,
            "events_per_type": events_per_type,
            "event_refs": event_refs,
            "event_rates": event_rates, 
            "cumulative_rates": cumulative_rates,
            "persistence_array": persistence_array,
           }

def RandomConfigurationIdeal( int nsites, double temperature ):
    """ Generates a random configuration commensurate with the temperature for an ideal lattice gas
       
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
