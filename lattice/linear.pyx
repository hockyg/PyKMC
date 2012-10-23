""" One-dimensional lattice

    Code for a one-dimensional lattice where neighbors are just
        adjacent sites.

"""

lattice_name = "linear"

import numpy as np
cimport numpy as np
from math import floor
import ctypes as ct

def Neighbors(nsites): 
    """ Calculates all neighbors for all sites """
    # note that there are exactly 2*nsites neighbors
    cdef int nneighbors_per_site = 2
    cdef int site_idx, j
    cdef np.ndarray[np.int_t,ndim=2] neighbors = np.zeros((nsites,nneighbors_per_site),dtype=ct.c_int)
    for site_idx in range(nsites):
        neighbors_i = NeighborsI(site_idx, nsites)
        for j in range(nneighbors_per_site):
            neighbors[site_idx,j] = neighbors_i[j]
    return nneighbors_per_site, neighbors

cdef NeighborsI( int site_idx, int nsites ):
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

def PeriodicDistance( site_idx1, site_idx2, nsites ):
    dx = 1.*( site_idx2 - site_idx1 )
    return dx - nsites*floor( dx/nsites + 0.5 )
