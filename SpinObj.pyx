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
has_dual = ["Plaquette"]

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
    def __init__(self):
        pass

    def initialize_new(self,lattice_name,model_name,side_length,temperature,max_time,seed=0):
        self.initial_configuration = None
        self.dual_configuration = None
        self.command_line_options = None
        self.final_options = None
        self.side_length = side_length

        self.lattice_name = lattice_name
        self.model_name = model_name

        #self.reset_for_continue(reset_time=True)

        self.max_time = max_time

        # set up lattice
        model = ModelRegistry[self.model_name]

        lattice = model.LatticeRegistry[lattice_name](side_length)
        nneighbors_per_site,nneighbors_update_per_site, neighbors, neighbors_update = lattice.Neighbors()
        self.nsites = nsites = lattice.nsites

        if self.model_name in has_dual:
            self.configuration, self.dual_configuration = lattice.RandomConfiguration( temperature )
        else:
            self.configuration = lattice.RandomConfiguration( temperature )
            self.dual_configuration = self.configuration # shouldn't take up any space, just pointer to same array

        self.initial_configuration = self.configuration.copy()
        self.prev_configuration = self.configuration.copy()
        if self.model_name in ["FA","East"]:
            self.initial_nonexcited = (self.initial_configuration<0.1).sum() # zeros
        elif self.model_name in ["Plaquette"]:
            self.initial_nonexcited = (self.initial_configuration>0.1).sum() # ones
        else:
            print "In SpinObj.pyx, have not defined which spin values are non-excited"
            sys.exit(2)


        # set up system object
        self.system = SpinSys()
        self.system.nsites = nsites
        self.system.nneighbors_per_site = nneighbors_per_site
        self.system.nneighbors_update_per_site = nneighbors_update_per_site
        self.system.neighbors = neighbors
        self.system.neighbors_update = neighbors_update
        self.system.model_number = model_dict[self.model_name]
        self.system.current_step = 0
        self.system.n_possible_events = 0
        self.system.seed = self.seed = seed
        self.system.temp = temperature
        self.system.betaexp = np.exp(-1./temperature)
        self.system.time = 0.0

        self.system.configuration = self.configuration
        self.system.prev_configuration = self.prev_configuration
        self.system.dual_configuration = self.dual_configuration

        sys_arrays = ModelRegistry[self.model_name].InitializeArrays( self.nsites )
        for key in sys_arrays.keys():
            setattr( self.system, key, sys_arrays[key] )

        self.system.persistence_array[self.configuration == 0] = -1 # down spins are -1, up spins are 1

#    def reset_for_continue(self,reset_time=True):
#        self.steps_to_take = []
#
#        if reset_time:
#            self.frame_step = 0
#            self.frame_time = 0
#
#        self.frame_num = 0
#
#        self.trj_file = None
#        self.start_file = None
#        self.restart_file = None
#        self.restarted_from = None

    def write_frame(self):
        pickle.dump(self.system, self.trj_file )

    def setup_output_files(self,mode="w",compresslevel=3):
        if self.final_options.trj_time > 0:
            self.trj_file_name = self.final_options.output_prefix+'.spintrj.gz'
            self.trj_file = gzip.open(self.trj_file_name,mode+'b',compresslevel)
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
        self.exception_list = []

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
                ("nneighbors_update_per_site",c_int), 
                ("neighbors",c_void_p),
                ("neighbors_update",c_void_p),
        # simulation stuff
                ("model_number",c_int),
                ("current_step",c_int),
                ("n_possible_events",c_int),
                ("seed",c_int),
                ("temp",c_double),
                ("betaexp",c_double), 
                ("time",c_float), 
                ("configuration",c_void_p),
                ("prev_configuration",c_void_p),
                ("dual_configuration",c_void_p),
        # rate stuff
                ("total_rate",c_float),
                ("events",c_void_p),
                ("event_refs",c_void_p),
                ("event_rates",c_void_p),
                ("event_ref_rates",c_void_p),
                ("cumulative_rates",c_void_p),
        # calculation stuff
                ("persistence_array",c_void_p),
    ]

SimData_p = ct.POINTER(SimData)
SimDataFields = {}
for name, t in SimData._fields_:
    SimDataFields[name] = True
