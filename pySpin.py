#!/usr/bin/python
import pyximport; pyximport.install()
import os
import sys
from SpinObj import *
import cPickle as pickle
import numpy as np

libname = "pySpin.so"
verbose_out = None

def simulate(options):
    timer = Timer()
    lattice = LatticeRegistry[options.lattice]
    model = ModelRegistry[options.model]

    if options.output_prefix:
        output_pickle = options.output_prefix+'_results.pickle'
        pickle.dump(results,open(output_pickle,'wb'), protocol=-1)

    print >>verbose_out, textline_box("Initializing simulator")
    # set up c functions
    pysimlib = os.path.join(os.getcwd(),libname)
    if not os.path.exists(pysimlib):
        sourcedir=os.path.abspath(os.path.dirname(sys.argv[0]))
        pysimlib=os.path.join(sourcedir,libname)

    print >>verbose_out, "Using simulation object:",pysimlib
    C = ct.CDLL(pysimlib)
    C.setup_spin_system.restype = c_int
    C.setup_spin_system.argtypes = (SimData_p,)
    C.cleanup_spin_system.restype = c_int
    C.cleanup_spin_system.argtypes = (SimData_p,)
    C.run_kmc_spin.restype = c_int
    C.run_kmc_spin.argtypes=(c_int,SimData_p) # firt arg is number of steps

    if options.seed is not None:
        seed=options.seed
    else:
        import random
        seed=random.SystemRandom().randint(0,10000000)
    np.random.seed(seed)

    if options.input is not None:
        print "Reading input not currently supported"
        sys.exit(1)
    else:
        initial_configuration = model.RandomConfiguration( options.nsites, options.temperature )

    nneighbors_per_site, neighbors = lattice.Neighbors(options.nsites)
 
    simulation = Simulation( lattice.lattice_name, model.model_name, options.nsites,
                             options.max_steps, configuration=initial_configuration )
    simulation.seed = seed

def main():
    from optparse import OptionParser, OptionGroup
    parser = OptionParser()
    parser.add_option('-i','--input',default=None,
                      help="Read in stored configuration (default: use random)")
    parser.add_option('-m', '--model', default="FA", 
                      help="Model to simulate (default: %default)" )
    parser.add_option('-l', '--lattice', default="linear", 
                      help="Model to simulate (default: %default)" )
    parser.add_option('-T', '--temperature', default=1.0, type=float,
                      help="Temperature to use (default: %default)" )
    parser.add_option('-N', '--nsites', default=10, type=int,
                      help="Number of sites to simulate (default: %default)" )
    parser.add_option('-s', '--max_steps', default=100, type=int,
                      help="Maximum number of steps to simulate (default: %default)" )
    parser.add_option('--seed', default=1, type=int,
                      help="Random seed for simulation (default: %default)" )
    parser.add_option("-o","--output_prefix",default=None,help="Set prefix for output files (default:none)")
    parser.add_option("-v","--verbose",help="Print useful execution information",dest="verbose",default=True,action="store_true")
    options, args = parser.parse_args()
 
    # add a more complicated option for this later, once deciding on writing out
    info_steps = options.max_steps/10
 
    if options.seed < 1:
        parser.error("Seed must be assigned a positive value")
 
    if not options.model in ModelRegistry:
        print "Model not yet defined. Please select from:"
        print "\t",ModelRegistry.keys()
        parser.print_help()
        sys.exit(1)
 
    if not options.lattice in LatticeRegistry:
        print "Lattice not yet defined. Please select from:"
        print "\t",LatticeRegistry.keys()
        parser.print_help()
        sys.exit(1)

    if options.verbose is False:
        verbose_out = NullDevice()
    else:
        verbose_out = sys.stdout

    if options.output_prefix:
        tee_logfile( options.output_prefix )
 
    simulate(options)

if __name__ == "__main__":
    main()
