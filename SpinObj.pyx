import os
import sys
import numpy as np
cimport numpy as np
import copy
import time
import datetime
import cPickle as pickle
import ctypes as ct
import gzip

from SpinUtil import *

LatticeRegistry = {}
ModelRegistry = {}

# try to make this more general. perhaps move later?
from lattice import linear
from models import FA

LatticeRegistry[linear.lattice_name] = linear
ModelRegistry[FA.model_name] = FA

model_dict = {"FA":0}

#definitions
c_int = ct.c_int
c_double = ct.c_double
c_float = ct.c_float
c_double_p = ct.c_double
c_void_p = ct.c_void_p

class Timer:
    def __init__(self):
        self.starttime = time.time()
    def gettime(self):
        return time.time() - self.starttime
    def time(self):
        print self.gettime()

# Taken from: http://shallowsky.com/blog/programming/python-tee.html
class tee :
    def __init__(self, _fd1, _fd2) :
        self.fd1 = _fd1
        self.fd2 = _fd2

    def __del__(self) :
        if self.fd1 != sys.stdout and self.fd1 != sys.stderr :
            self.fd1.close()
        if self.fd2 != sys.stdout and self.fd2 != sys.stderr :
            self.fd2.close()

    def write(self, text) :
        self.fd1.write(text)
        self.fd2.write(text)

    def flush(self) :
        self.fd1.flush()
        self.fd2.flush()

class NullDevice():
    def write(self, s):
        pass

class Simulation(object):
    def __init__(self,lattice_name,model_name,nsites,max_steps,configuration=None):
        self.initial_configuration = None
        self.command_line_options = None
        self.final_options = None

        self.lattice_name = lattice_name
        self.model_name = model_name

        self.nsites = nsites

        self.reset_for_continue(reset_time=True)

        self.max_steps = max_steps
        if configuration is not None:
             self.initial_configuration = np.copy(configuration)

    def reset_for_continue(self,reset_time=True):
        self.steps_to_take = []

        if reset_time:
            self.frame_step = 0
            self.frame_time = 0

        self.frame_num = 0

        self.trj_file = None
        self.start_file = None
        self.restart_file = None
        self.restarted_from = None

    def setup_output_files(self):
        if self.final_options.trj_time > 0:
            self.trj_file_name = self.final_options.output_prefix+'.spintrj.gz'
        if self.final_options.restart_time > 0:
            self.restart_file_name = self.final_options.output_prefix+'.spinrestart.gz'

    def setup_steps_to_take(self):
        pass

    def save_state(self,filename,compresslevel=3):
        if not os.path.splitext(filename)[-1]==".gz":
            filename=filename+".gz"
        pickle.dump(self,gzip.open(filename,'wb',compresslevel),protocol=-1)

    def load_state(self,filename):
        tmp=pickle.load(gzip.open(filename,'rb'))
        for item in tmp.__dict__:
            setattr(self,item,getattr(tmp,item))
        setattr(self,"restarted_from",os.path.abspath(os.path.expanduser(filename)))

    def start_new_from_file(self,filename,reset_time=True):
        self.load_state(filename)
        self.reset_for_continue(reset_time=True)

#    def open_final_files(self,mode="w",compresslevel=3):
#    def open_start_files(self,mode="w",compresslevel=3):
    def print_state(self):
        for key in sorted(self.__dict__):
            if hasattr(self, key):
                print key, getattr(self, key)

class SpinSys(object):
    def __init__(self):
        import socket
        self.SD = SimData()
        self.SD_p = ct.pointer(self.SD)
        self.creation_date = datetime.datetime.now()
        self.created_on = socket.gethostname()
        self.creation_host_system_info = os.uname()

    def print_state(self):
        for key in sorted(SimDataFields):
            if hasattr(self, key):
                print key, getattr(self, key)

    def __getstate__(self):
        state = { }
        #for key in classvars:
        for key in SimDataFields:
            if hasattr(self, key):
                state[key] = getattr(self, key)

        return state

    def __setstate__(self, state):
        self.__init__()   # allocate arrays and stuff

        #for key in classvars:
        for key in SimDataFields:
            if not state.has_key(key):
                continue
            setattr(self, key, state[key])

        #from exceptions list
# e.g.       setattr(self,"coords_verlet",self.coords.copy())

    def __getattr__(self, attrname):
        if attrname in SimDataFields:
           if attrname in self.exception_list: return None

           return getattr(self.SD, attrname)
        else:
            raise AttributeError("No Such Attribute: %s"%attrname)

    def __setattr__(self, name, value):
        if type(value)==np.ndarray:
            if name in SimDataFields:
                setattr(self.SD,name,value.ctypes.data)
            self.__dict__[name] = value
        elif name in SimDataFields:
            setattr(self.SD, name, value)
        else:
            self.__dict__[name] = value

    def copy(self):
        return copy.copy(self)

class SimData(ct.Structure):
    _fields_=[
        #lattice stuff
                ("nsites",c_int),
                ("nneighbors_per_site",c_int), 
                ("neighbors",c_void_p),
        # simulation stuff
                ("model_number",c_int),
                ("current_step",c_int),
                ("n_possible_events",c_int),
                ("seed",c_int),
                ("temp",c_double),
                ("betaexp",c_double), 
                ("time",c_float), 
                ("configuration",c_void_p),
        # rate stuff
                ("events",c_void_p),
                ("event_refs",c_void_p),
                ("event_rates",c_void_p),
                ("event_ref_rates",c_void_p),
                ("cumulative_rates",c_void_p),
        # storage stuff
                ("event_storage",c_void_p),
                ("persistence",c_void_p),
    ]

SimData_p = ct.POINTER(SimData)
SimDataFields = {}
for name, t in SimData._fields_:
    SimDataFields[name] = True

########################################
#Old version


#cdef SumRates( np.ndarray[np.float_t,ndim=1] rates, np.ndarray[np.float_t,ndim=1] cumulative_rates,int n_possible_rates ):
#    cdef i
#    cdef double total_rate = 0.0
#    for i in range(n_possible_rates):
#        total_rate = total_rate + rates[i]
#        if i>0:
#            cumulative_rates[i] = cumulative_rates[i-1]+rates[i]
#        else:
#            cumulative_rates[i] = rates[i]
#    return total_rate
#
#def driveKCM( simulation, int info_steps, float temperature, int seed = 0 ):
#    cdef int nstages
#    cdef int n_possible_events
#    cdef int step_i, stage_index
#    cdef float t = 0.0
#    cdef float prob, total_rate
#
#    cdef double beta = 1/temperature
#    cdef double betaexp = np.exp(-beta)
#
#    model = simulation.model
#    lattice = simulation.lattice 
#    cdef int max_steps = simulation.max_steps
#    cdef int nsites = simulation.nsites
#
#    np.random.seed( seed )
#
#    cdef np.ndarray random_floats = np.zeros( info_steps, dtype=np.float )
#
#    nneighbors_per_site, neighbors = lattice.Neighbors(nsites)
#    events, event_rates, event_refs, event_ref_rates, event_storage = simulation.initialize(temperature)
#    cdef np.ndarray initial_configuration = simulation.initial_configuration
#    cdef np.ndarray cumulative_rates = np.zeros(len(event_ref_rates),dtype=np.float)
#    cdef np.ndarray times = np.zeros(max_steps+1,dtype=np.float)
#
#    # initialize or load configuration
#    configuration = initial_configuration.copy()
#    n_possible_events = model.AllEvents( betaexp, events, event_rates, event_refs, event_ref_rates, configuration, nsites, neighbors, nneighbors_per_site)
#
#    for step_i in range( max_steps ):
#        stage_index = step_i % info_steps
#        if( step_i % info_steps == 0):
#            random_floats = np.random.random(size=info_steps)
#            print "%i / %i: %f"%( step_i, max_steps, t )
#            #print "\t",configuration
#        if n_possible_events < 1:
#            print "Terminating. No more moves possible."
#            break
#        prob = random_floats[ stage_index ]
#        total_rate = SumRates( event_ref_rates, cumulative_rates, n_possible_events )
#        event_i = BSearchProb( prob*total_rate, n_possible_events, cumulative_rates )
#        model.UpdateConfiguration( configuration, events, event_refs, event_i, event_storage, step_i )
#        n_possible_events = model.UpdateEventsI( event_refs[event_i], betaexp, events, event_rates, event_refs, event_ref_rates, configuration, nsites, neighbors, nneighbors_per_site)
#        dt = -np.log(prob)/total_rate
#        t += dt
#        times[step_i+1] = t 
#
#    return times, initial_configuration, event_storage
