import numpy as np
cimport numpy as np

from util import BSearchProb

LatticeRegistry = {}
ModelRegistry = {}

# try to make this more general. perhaps move later?
from lattice import linear
from models import FA

LatticeRegistry[linear.lattice_name] = linear
ModelRegistry[FA.model_name] = FA

class Simulation(object):
    def __init__(self,lattice,model,nsites,max_steps,configuration=None):
        self.lattice = lattice
        self.model = model
        self.nsites = nsites
        self.max_steps = max_steps
        self.initial_configuration = configuration

    def setup_lattice(self):
        nneighbors_per_site, neighbors = self.lattice.Neighbors(self.nsites)
        return nneighbors_per_site, neighbors

    def initialize(self, temperature):
        if self.initial_configuration is None:
            self.initial_configuration = self.model.RandomConfiguration( self.nsites, temperature )
        return self.model.InitializeArrays( self.nsites )

cdef SumRates( np.ndarray[np.float_t,ndim=1] rates, np.ndarray[np.float_t,ndim=1] cumulative_rates,int n_possible_rates ):
    cdef i
    cdef double total_rate = 0.0
    for i in range(n_possible_rates):
        total_rate = total_rate + rates[i]
        if i>0:
            cumulative_rates[i] = cumulative_rates[i-1]+rates[i]
        else:
            cumulative_rates[i] = rates[i]
    return total_rate

def driveKCM( simulation, int info_steps, float temperature, int seed = 0 ):
    cdef int nstages
    cdef int n_possible_events
    cdef int step_i, stage_index
    cdef float t = 0.0
    cdef float prob, total_rate

    cdef double beta = 1/temperature
    cdef double betaexp = np.exp(-beta)

    model = simulation.model
    lattice = simulation.lattice 
    cdef int max_steps = simulation.max_steps
    cdef int nsites = simulation.nsites

    np.random.seed( seed )

    cdef np.ndarray random_floats = np.zeros( info_steps, dtype=np.float )

    nneighbors_per_site, neighbors = lattice.Neighbors(nsites)
    events, event_rates, event_refs, event_ref_rates = simulation.initialize(temperature)
    cdef np.ndarray initial_configuration = simulation.initial_configuration
    cdef np.ndarray cumulative_rates = np.zeros(len(event_ref_rates),dtype=np.float)
    # initialize or load configuration
    configuration = initial_configuration.copy()
    n_possible_events = model.AllEvents( betaexp, events, event_rates, event_refs, event_ref_rates, configuration, nsites, neighbors, nneighbors_per_site)

    for step_i in range( max_steps ):
        stage_index = step_i % info_steps
        if( step_i % info_steps == 0):
            random_floats = np.random.random(size=info_steps)
            print "%i / %i: %f"%( step_i, max_steps, t )
            print "\t",configuration
        if n_possible_events < 1:
            print "Terminating. No more moves possible."
            break
        prob = random_floats[ stage_index ]
        total_rate = SumRates( event_ref_rates, cumulative_rates, n_possible_events )
        event_i = BSearchProb( prob*total_rate, n_possible_events, cumulative_rates )
        model.UpdateConfiguration( configuration, events, event_refs, event_i )
        n_possible_events = model.UpdateEventsI( event_refs[event_i], betaexp, events, event_rates, event_refs, event_ref_rates, configuration, nsites, neighbors, nneighbors_per_site)
        dt = np.log(1/prob)/total_rate
        t += dt
