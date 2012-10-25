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
    model = ModelRegistry[options.model]

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
    C.get_event_rate.restype = c_float
    C.get_event_rate.argtypes = (c_int, SimData_p,)

    if options.seed is not None:
        seed=options.seed
    else:
        import random
        seed=random.SystemRandom().randint(0,10000000)
    np.random.seed(seed)
    print >>verbose_out,"\tUsing seed: %i"%seed

    if options.input is not None:
        print "Reading input not currently supported"
        sys.exit(1)
    else:
        simulation = Simulation()
        simulation.initialize_new( options.lattice, options.model, options.side_length,
                                   options.temperature, options.max_steps, seed=seed )
        simulation.command_line_options = options
        # set up write, restart, and info times
        options.trj_time = 1
        options.restart_time = 1

        simulation.final_options = options
        if options.output_prefix:
            simulation.setup_output_files()

    try:
        print >>verbose_out, textline_box("Running simulation: (setup time = %f )"%timer.gettime())
        C.setup_spin_system(simulation.system.SD)
        p1,p2 = persistence( simulation.nsites, simulation.initial_nonexcited, simulation.system.persistence_array)
        
        print simulation.system.time, p1, p2, c_to_T_ideal( simulation.nsites, simulation.system.dual_configuration ),model.SquareEnergy( simulation.system.configuration, simulation.system.neighbors, simulation.system.nsites,simulation.system.nneighbors_per_site )
#        print simulation.system.configuration.reshape((simulation.side_length,simulation.side_length))
#        print simulation.system.dual_configuration.reshape((simulation.side_length,simulation.side_length))

        for i in range(simulation.max_steps):
#            print simulation.system.event_ref_rates[:simulation.system.n_possible_events]
            return_val = C.run_kmc_spin(1,simulation.system.SD)
#            np.savetxt(sys.stdout,simulation.system.event_rates,fmt="%3.2f")
            if return_val == -1: 
                print "No more possible moves"
                break
            p1,p2 = persistence( simulation.nsites, simulation.initial_nonexcited, simulation.system.persistence_array)
            if options.output_prefix:
                simulation.write_frame()
            print "Time: %.2e"%simulation.system.time,p1,p2, c_to_T_ideal( simulation.nsites, simulation.system.dual_configuration ), model.SquareEnergy( simulation.system.configuration, simulation.system.neighbors, simulation.system.nsites,simulation.system.nneighbors_per_site )

#            tmp_configuration = simulation.system.configuration.copy()
#            for i in range(simulation.nsites):
#                tmp_configuration[i] = -1*tmp_configuration[i]
#                if np.abs(C.get_event_rate(i,simulation.system.SD)-min(1, np.exp(-1*(model.SquareEnergy( tmp_configuration, simulation.system.neighbors, simulation.system.nsites,simulation.system.nneighbors_per_site ) - model.SquareEnergy( simulation.system.configuration, simulation.system.neighbors, simulation.system.nsites,simulation.system.nneighbors_per_site ))/simulation.system.temp)))>0.001: print i, C.get_event_rate(i,simulation.system.SD), model.SquareEnergy( tmp_configuration, simulation.system.neighbors, simulation.system.nsites,simulation.system.nneighbors_per_site ) - model.SquareEnergy( simulation.system.configuration, simulation.system.neighbors, simulation.system.nsites,simulation.system.nneighbors_per_site )
#                tmp_configuration[i] = -1*tmp_configuration[i]


#            print simulation.system.configuration.reshape((simulation.side_length,simulation.side_length))
#            print simulation.system.dual_configuration.reshape((simulation.side_length,simulation.side_length))

        print >>verbose_out, "Simulation Finished!"
        C.cleanup_spin_system(simulation.system.SD)
    except KeyboardInterrupt:
        print "Terminating simulation..."
        C.cleanup_spin_system(simulation.system.SD)

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
    parser.add_option('-L', '--side_length',dest="side_length",default=10, type=int,
                      help="Side length of linear, square or cubic lattice to simulate. Total number of sites for other lattices, if ever implemented (default: %default)" )
    parser.add_option('-s', '--max_steps', default=100, type=int,
                      help="Maximum number of steps to simulate (default: %default)" )
    parser.add_option('--seed', default=None, type=int,
                      help="Random seed for simulation (default: %default)" )
    parser.add_option("-o","--output_prefix",default=None,help="Set prefix for output files (default:none)")
    parser.add_option("-v","--verbose",help="Print useful execution information",dest="verbose",default=True,action="store_true")
    options, args = parser.parse_args()
 
    # add a more complicated option for this later, once deciding on writing out
    info_steps = options.max_steps/10
 
    if options.seed is not None and options.seed < 1:
        parser.error("Seed must be assigned a positive value")
 
    if not options.model in ModelRegistry:
        print "Model not yet defined. Please select from:"
        print "\t",ModelRegistry.keys()
        parser.print_help()
        sys.exit(1)

    LatticeRegistry = ModelRegistry[options.model].LatticeRegistry
 
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
