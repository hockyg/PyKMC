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
    C.run_kmc_spin.argtypes=(c_double,SimData_p) # firt arg is stop_time
#    C.get_event_type.restype = c_int
#    C.get_event_type.argtypes = (c_int, SimData_p,)

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
        simulation.initialize_new( options.lattice, options.model, options.dynamics_type,
                                   options.side_length, options.temperature, 
                                   options.max_time, seed=seed )
        simulation.command_line_options = options

        simulation.final_options = options
        if options.output_prefix:
            simulation.setup_output_files()

        if options.info_time <= 0 or options.info_time > options.max_time:
            options.info_time = options.max_time

    try:
        print >>verbose_out, textline_box("Running simulation: (setup time = %f )"%timer.gettime())
        C.setup_spin_system(simulation.system.SD)
        average_time_per_step = 1/simulation.system.SD.total_rate

        # now set up desired write out times
        if options.linear_time:
            simulation.stop_times = np.array(np.arange(options.info_time,options.max_time+options.info_time,options.info_time),dtype=c_double)
            if simulation.stop_times[-1] > options.max_time:
                simulation.stop_times[-1] = options.max_time
            simulation.nstages = len(simulation.stop_times)
            
            pass
        else: # log time
            min_time = average_time_per_step
            min_time_log = int(np.ceil(np.log10(min_time)))
            max_time_log = int(np.ceil(np.log10(options.max_time)))
            prelim_stop_times = np.logspace( min_time_log, max_time_log, num=(max_time_log-min_time_log)*options.stops_per_decade+1)
            simulation.stop_times = prelim_stop_times[prelim_stop_times<=options.max_time]
            simulation.nstages = len(simulation.stop_times)
          
        # before starting, write initial frame
        if options.output_prefix and options.write_trj:
            simulation.write_frame()

        #p1,p2 = persistence( simulation.nsites, simulation.initial_nonexcited, simulation.system.persistence_array)
        
        #print "Time: %.2e"%simulation.system.time, p1, p2, c_to_T_ideal( simulation.nsites, simulation.system.dual_configuration ),model.SquareEnergy( simulation.system.configuration, simulation.system.neighbors, simulation.system.nsites,simulation.system.nneighbors_per_site )
#        print "Time: %.2e"%simulation.system.time, c_to_T_ideal( simulation.nsites, simulation.system.dual_configuration ),model.SquareEnergy( simulation.system.configuration, simulation.system.neighbors, simulation.system.nsites,simulation.system.nneighbors_per_site )
#uncomment for newest
        print "Time: %.2e Energy: %f"%( simulation.system.time, simulation.system.total_energy ),c_to_T_ideal( simulation.nsites, simulation.system.dual_configuration )

        sim_timer = Timer()
        prev_time = sim_timer.gettime()
        prev_stop_time = 0
        total_steps=0L
        last_time = simulation.stop_times[-1]
        for frame_idx,stop_time in enumerate(simulation.stop_times):
            #print simulation.system.configuration.reshape((8,-1))
            #print simulation.system.dual_configuration.reshape((8,-1))
            #print simulation.system.event_types.reshape((8,-1))
            #print "Total rate:", simulation.system.total_rate, "Should be:", (simulation.system.events_per_type*simulation.system.event_rates).sum()
            #print simulation.system.events_per_type
            #print simulation.system.event_rates
            return_val = C.run_kmc_spin(stop_time, simulation.system.SD)
            if return_val == -1: 
                print "No more possible moves"
                break
            #p1,p2 = persistence( simulation.nsites, simulation.initial_nonexcited, simulation.system.persistence_array)
            if options.output_prefix and options.write_trj:
                simulation.write_frame()
            #print "Time: %.2e"%simulation.system.time,p1,p2, c_to_T_ideal( simulation.nsites, simulation.system.dual_configuration ), model.SquareEnergy( simulation.system.configuration, simulation.system.neighbors, simulation.system.nsites,simulation.system.nneighbors_per_site )
#            print "Time: %.2e"%simulation.system.time, c_to_T_ideal( simulation.nsites, simulation.system.dual_configuration ), model.SquareEnergy( simulation.system.configuration, simulation.system.neighbors, simulation.system.nsites,simulation.system.nneighbors_per_site )

            elapsed_time = stop_time - prev_stop_time
            avg_dt = elapsed_time/simulation.system.SD.current_step
            elapsed_time = sim_timer.gettime()
            time_remaining = last_time - stop_time
#            The next line is a timing estimate that works, but the replaced code gives a better estimate
#            est_final_sim_time = elapsed_time / ( 1 - (time_remaining)/last_time )
             # then this is remaining time guess: est_final_sim_time - elapsed_time
            stage_elapsed_time = sim_timer.gettime()-prev_time
            wall_time_per_sim_time = stage_elapsed_time / ( stop_time - prev_stop_time )
            wall_time_remaining = time_remaining * wall_time_per_sim_time
            total_steps = total_steps+simulation.system.SD.current_step
            efficiency = simulation.system.SD.current_step/stage_elapsed_time

            #this should be last major thing in loop
            prev_time = sim_timer.gettime()
            prev_stop_time = stop_time

#uncomment for newest
            print "Time: %.2e Energy: %f Dt: %3.2e SimTime: %f (etr: %f) Eff: %3.2e"%( simulation.system.time, simulation.system.total_energy, avg_dt, elapsed_time, wall_time_remaining, efficiency ), c_to_T_ideal( simulation.nsites, simulation.system.dual_configuration )

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
    parser.add_option('--dynamics_type', default="Metropolis",
                      help="Type of dynamics to run (default: %default)")
    parser.add_option('-l', '--lattice', default="linear", 
                      help="Model to simulate (default: %default)" )
    parser.add_option('-T', '--temperature', default=1.0, type=float,
                      help="Temperature to use (default: %default)" )
    parser.add_option('-L', '--side_length',dest="side_length",default=10, type=int,
                      help="Side length of linear, square or cubic lattice to simulate. Total number of sites for other lattices, if ever implemented (default: %default)" )

#    parser.add_option('-s', '--max_steps', default=100, type=int,
#                      help="Maximum number of steps to simulate (default: %default)" )
    parser.add_option('-t', '--max_time', default=1, type=float,
                      help="Maximum time up to which to simulate (default: %default)" )
    parser.add_option('--seed', default=None, type=int,
                      help="Random seed for simulation (default: %default)" )

    write_group=OptionGroup(parser,"Options for writing out","Options for writing out information to the screen and files")
    write_group.add_option("-o","--output_prefix",default=None,help="Set prefix for output files (default:none)")
    write_group.add_option("-v","--verbose",help="Print useful execution information",dest="verbose",default=True,action="store_true")
    write_group.add_option("--linear_time",help="Save information on a linear time scale (default: logarithmic)",default=False,action="store_true")
    write_group.add_option("--info_time",default=-1,type=float,help="Set how often simulation info and configurations are written for linear time (default:none)")
    write_group.add_option("--stops_per_decade",default=10,type=int,help="For logarithmic writing. Set how many times to write per decade of simulation time (default:none)")
    write_group.add_option("--write_trj",default=False,action="store_true",help="Specify whether or not to write a trajectory (default:False)")

    parser.add_option_group(write_group)
    options, args = parser.parse_args()
 
    if options.seed is not None and options.seed < 1:
        parser.error("Seed must be assigned a positive value")

    if not options.model in ModelRegistry:
        print "Model not yet defined. Please select from:"
        print "\t",sorted(ModelRegistry.keys())
        parser.print_help()
        sys.exit(1)

    LatticeRegistry = ModelRegistry[options.model].LatticeRegistry
 
    if not options.lattice in LatticeRegistry:
        print "Lattice not yet defined. Please select from:"
        print "\t",sorted(LatticeRegistry.keys())
        parser.print_help()
        sys.exit(1)

    if not options.dynamics_type in dynamics_dict:
        print "Dynamics type selected not available. Please select from:"
        print "\t",sorted(dynamics_dict.keys())
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
