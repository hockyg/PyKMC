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

#LatticeRegistry = {}
#LatticeRegistry[linear.lattice_name] = linear
ModelRegistry = {}

# try to make this more general. perhaps move later?
from models import FA, East, Plaquette

ModelRegistry[FA.model_name] = FA
ModelRegistry[East.model_name] = East
ModelRegistry[Plaquette.model_name] = Plaquette

model_dict = {"FA":0, "East":1, "Plaquette":10}
dynamics_dict = {"Metropolis":0, "Glauber":1}
has_dual = ["Plaquette"]
frozen_geometries = {None:0,"CAVITY":1,"WALL":2,"SANDWICH":3,"RANDOM":4}

#definitions
c_int = ct.c_int
c_long = ct.c_long
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
    def __init__(self):
        self.set_host_info()

    def set_host_info(self):
        import socket
        self.creation_date = datetime.datetime.now()
        self.created_on = socket.gethostname()
        self.creation_host_system_info = os.uname()

    # count initial non-excited spins, and also set persistence array
    def set_initial_nonexcited(self):
        self.system.persistence_array[:] = 1
        if self.model_name in ["FA","East"]:
            self.initial_nonexcited = (self.initial_configuration<0.1).sum() # zeros
            self.system.persistence_array[self.configuration == 0] = -1 # down spins are -1, up spins are 1
        elif self.model_name in ["Plaquette"]:
            self.initial_nonexcited = (self.initial_configuration>0.1).sum() # ones
            self.system.persistence_array[self.dual_configuration == 0] = -1 # down spins are -1, up spins are 1
        else:
            print "In SpinObj.pyx, have not defined which spin values are non-excited"
            sys.exit(2)

    def initialize_new(self,lattice_name,model_name,dynamics_type,linear_size,temperature,max_time,seed=0,activelist=None,input_cfgs=None):
        self.initial_configuration = None
        self.dual_configuration = None
        self.command_line_options = None
        self.final_options = None
        self.linear_size = linear_size

        self.lattice_name = lattice_name
        self.model_name = model_name


        self.max_time = max_time

        # set up lattice
        model = ModelRegistry[self.model_name]

        lattice = model.LatticeRegistry[lattice_name](linear_size)
        nneighbors_per_site,nneighbors_update_per_site, neighbors, neighbors_update = lattice.Neighbors()
        self.nsites = nsites = lattice.nsites
        if activelist is not None:
            self.nactive = len(activelist)
        else:
            self.nactive = self.nsites

        self.n_event_types = n_event_types = lattice.n_event_types
        event_rates = lattice.EventRates(temperature,dynamics_dict[dynamics_type])

        if input_cfgs is None:
            if self.model_name in has_dual:
                self.configuration, self.dual_configuration = lattice.RandomConfiguration( temperature )
            else:
                self.configuration = lattice.RandomConfiguration( temperature )
                self.dual_configuration = self.configuration # shouldn't take up any space, just pointer to same array
        else:
            if self.model_name in has_dual:
                self.configuration, self.dual_configuration = input_cfgs
            else:
                self.configuration = input_cfgs[0]
                self.dual_configuration = self.configuration # shouldn't take up any space, just pointer to same array
            

        self.initial_configuration = self.configuration.copy()
        self.prev_configuration = self.configuration.copy()

        # set up system object
        self.system = SpinSys()
        self.system.nsites = nsites
        self.system.nactive = self.nactive
        self.system.nneighbors_per_site = nneighbors_per_site
        self.system.nneighbors_update_per_site = nneighbors_update_per_site
        self.system.neighbors = neighbors
        self.system.neighbors_update = neighbors_update
        self.system.model_number = model_dict[self.model_name]
        self.system.current_step = 0
        self.system.n_possible_events = 0
        self.system.n_event_types = self.n_event_types
        self.system.seed = self.seed = seed
        self.system.temp = temperature
        self.system.betaexp = np.exp(-1./temperature)
        self.system.total_energy = 0.0
        self.system.time = 0.0

        self.system.configuration = self.configuration
        self.system.prev_configuration = self.prev_configuration
        self.system.dual_configuration = self.dual_configuration

        sys_arrays = ModelRegistry[self.model_name].InitializeArrays( self.nsites, self.n_event_types, self.nactive )
        sys_arrays["event_rates"] = np.array( event_rates, dtype=c_double )
        for key in sys_arrays.keys():
            setattr( self.system, key, sys_arrays[key] )
        self.set_initial_nonexcited()

        cdef int i
        # overwrite activelist if it is already specified. then use activelist to set isactivelist
        if activelist is not None:
            # first clear isactivelist by making everything inatctive
            for i in range(self.nsites):
                self.system.isactivelist[i] = 0
            for i in range(self.nactive):
                self.system.activelist[i] = activelist[i]
                self.system.isactivelist[self.system.activelist[i]] = 1

    def reset_for_continue(self,reset_time=True):
        self.set_host_info()
        self.stop_times = None
        if reset_time:
            self.system.current_step = 0
            self.system.time = 0

        self.trj_file_name = None
        self.trj_file = None

        #This is a good place to check that configuration, initial configuration, and dual_configuration all work as expected. for example, the next line works to make sure  the configuration and dual configuration match, but not the initial configuration and dual configuration (unless restarting from a state where that was true)
        #print "tpm test 1: (expect 0)",Plaquette.test_triangle_dual( self.configuration.reshape((self.linear_size,-1)), self.dual_configuration.reshape((self.linear_size,-1)), self.linear_size )
        #print "tpm test 2: (expect >0)",Plaquette.test_triangle_dual( self.initial_configuration.reshape((self.linear_size,-1)), self.dual_configuration.reshape((self.linear_size,-1)), self.linear_size )
        self.prev_configuration = self.configuration.copy()
        self.initial_configuration = self.configuration.copy()
        self.set_initial_nonexcited()

    def freeze_random(self, frozen_fraction):
        frac_active = 1 - frozen_fraction
        nactive = int(np.round( frac_active*self.nsites ))
        activelist = np.arange(self.nsites,dtype=c_int)
        np.random.shuffle( activelist )
        activelist = activelist[:nactive]
        activelist.sort()
        self.set_active(activelist)
    
    def set_active(self, activelist):
        nactive = len(activelist)
        isactivelist = np.zeros(self.nsites,dtype=c_int)
        isactivelist[activelist] = 1
        self.nactive = self.system.nactive = nactive
        self.system.activelist = activelist
        self.system.isactivelist = isactivelist

    def write_frame(self):
        pickle.dump(self.system.get_frame_state(only_active=True), self.trj_file )

    def setup_output_files(self,mode="w",compresslevel=3):
        if self.final_options.output_prefix is None: return
        if self.final_options.write_trj > 0:
            self.trj_file_name = self.final_options.output_prefix+'.spintrj.gz'
            self.trj_file = gzip.open(self.trj_file_name,mode+'b',compresslevel)
#        if self.final_options.restart_time > 0:
#            self.restart_file_name = self.final_options.output_prefix+'.spinrestart.gz'

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

    def print_state(self):
        for key in sorted(self.__dict__):
            if hasattr(self, key) and key != "trj_file":
                print key, getattr(self, key)

class SpinSys(object):
    def __init__(self):
        import socket
        self.SD = SimData()
        self.SD_p = ct.pointer(self.SD)
        # Note, these are automatically reset every time a new simulation is started, even in restart mode
        self.creation_date = datetime.datetime.now()
        self.created_on = socket.gethostname()
        self.creation_host_system_info = os.uname()
        self.exception_list = []
        self.stop_time = 0.0
        # exceptions for saving a single frame/configuration
        self.frame_exception_list = ["neighbors", "neighbors_update",                                                    "events", "event_types", "events_by_type",
                                     "events_per_type", "event_refs",
                                     "event_rates", "cumulative_rates",
                                     ]
        self.save_fields = ["creation_date","creation_host_system_info","created_on", "stop_time"]

    def get_frame_state(self,only_active=False):
        state = { }
        for key in SimDataFields.keys()+self.save_fields:
            if hasattr(self, key) and key not in self.frame_exception_list:
                if only_active and getattr(self,"nactive")<getattr(self,"nsites") and key.find("configuration")>-1:
                    tmp_cfg = getattr(self,key)[getattr(self,"activelist")]
                    state[key] = tmp_cfg
                else:
                    state[key] = getattr(self, key)

        return state

    def print_state(self):
        for key in sorted(self.__dict__):
        #for key in sorted(SimDataFields):
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
                ("nactive",c_int),
                ("nneighbors_per_site",c_int), 
                ("nneighbors_update_per_site",c_int), 
                ("activelist",c_void_p),
                ("isactivelist",c_void_p),
                ("neighbors",c_void_p),
                ("neighbors_update",c_void_p),
        # simulation stuff
                ("model_number",c_int),
                ("current_step",c_long),
                ("n_possible_events",c_int),
                ("seed",c_int),
                ("temp",c_double),
                ("betaexp",c_double), 
                ("total_energy",c_double), 
                ("time",c_double), 
                ("configuration",c_void_p),
                ("prev_configuration",c_void_p),
                ("dual_configuration",c_void_p),
        # rate stuff
                ("total_rate",c_double),
                ("n_event_types",c_int),
                ("events",c_void_p),
                ("event_types",c_void_p),
                ("events_by_type",c_void_p),
                ("events_per_type",c_void_p),
                ("event_refs",c_void_p),
                ("event_rates",c_void_p),
#                ("event_ref_rates",c_void_p),
                ("cumulative_rates",c_void_p),
        # calculation stuff
                ("persistence_array",c_void_p),
    ]

SimData_p = ct.POINTER(SimData)
SimDataFields = {}
for name, t in SimData._fields_:
    SimDataFields[name] = True
