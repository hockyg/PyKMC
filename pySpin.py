#!/usr/bin/python
import pyximport; pyximport.install()
import os
import sys
from SpinObj import *
import cPickle as pickle
import numpy as np

libname = "pySpin.so"
verbose_out = None

def print_start_options(options,simulation):
    print >>verbose_out,"Lattice parameters:"
    print >>verbose_out,"\tmodel: %s"%(options.model)
    print >>verbose_out,"\tlattice: %s"%(options.lattice)
    print >>verbose_out,"\tlinear_size = %i (nsites = %i)"%(simulation.side_length,simulation.nsites)
    print >>verbose_out

    print >>verbose_out,"Simulation parameters:"
    print >>verbose_out,"\tT = %f (beta = %f)"%(simulation.system.temp, 1./simulation.system.temp)
    print >>verbose_out,"\tmaximum_time = %e"%(simulation.max_time)
    print >>verbose_out,"\tdynamics_type: %s"%options.dynamics_type
    print >>verbose_out,"\tseed: %i"%simulation.seed
    if options.output_prefix:
        print >>verbose_out,"\toutput_prefix: %s"%(options.output_prefix) 
        print >>verbose_out,"\twriting_trajectory:",options.write_trj
    print >>verbose_out

def simulate(options):
    timer = Timer()
    model = ModelRegistry[options.model]

    print >>verbose_out, textline_box("Initializing simulator")
    # set up c functions
    pysimlib = os.path.join(os.getcwd(),libname)
    if not os.path.exists(pysimlib):
        sourcedir=os.path.abspath(os.path.dirname(sys.argv[0]))
        pysimlib=os.path.join(sourcedir,libname)

    print >>verbose_out, "Using simulation object: %s\n"%(pysimlib)
    C = ct.CDLL(pysimlib)
    C.setup_spin_system.restype = c_int
    C.setup_spin_system.argtypes = (SimData_p,)
    C.cleanup_spin_system.restype = c_int
    C.cleanup_spin_system.argtypes = (SimData_p,)
    C.run_kmc_spin.restype = c_int
    C.run_kmc_spin.argtypes=(c_double,SimData_p) # firt arg is stop_time

    if options.seed is not None:
        seed=options.seed
    else:
        import random
        seed=random.SystemRandom().randint(0,10000000)
    np.random.seed(seed)

    if options.restart is not None:
        simulation = load_object(options.restart)
        command_line_options = copy.copy(options)
        options = simulation.final_options
        if simulation.system.time == 0:
            openmode='w' 
        else:
            openmode='a'
        simulation.setup_output_files(mode=openmode)
        print "Doing full restart with:"
        np.random.seed(simulation.seed)
        print_start_options(options,simulation)

    elif options.input is not None:
        simulation = load_object(options.input)
        command_line_options = copy.copy(options)
        options = simulation.final_options
        
        # Note, this will have a broken trj_file attribute
        simulation.print_state()
        #remove this later
        simulation.setup_output_files()
    else:
        simulation = Simulation()
        simulation.initialize_new( options.lattice, options.model, options.dynamics_type,
                                   options.side_length, options.temperature, 
                                   options.max_time, seed=seed )
        simulation.command_line_options = options
        simulation.seed = seed

        simulation.final_options = options
        if options.output_prefix:
            simulation.setup_output_files()

        if options.info_time <= 0 or options.info_time > options.max_time:
            options.info_time = options.max_time
        
        print_start_options(options,simulation)

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
            if len(simulation.stop_times)==0:
                # if simulation time is less than estimated dt (happens for very small lattices)
                #     create an array of size 1 with just max time
                simulation.stop_times = np.array(np.ones(1)*options.max_time,dtype=c_double)
            simulation.nstages = len(simulation.stop_times)
          
        # before starting, write initial frame and start file
        if options.output_prefix and options.write_trj:
            simulation.write_frame()
            simulation.save_state(options.output_prefix+'.start.spinsim.gz')

        E_per_site = simulation.system.total_energy/simulation.system.nsites
        Teff = 1/np.log(1/E_per_site-1)
        print "Time: %.2e Dt: %.2e | Energy: %.2e Teff: %.4e"%( simulation.system.time, average_time_per_step, E_per_site, Teff )

        sim_timer = Timer()
        prev_time = sim_timer.gettime()
        prev_stop_time = 0
        total_steps=0L
        last_time = simulation.stop_times[-1]
        for frame_idx,stop_time in enumerate(simulation.stop_times):
            # last check for cases of restart to make sure time has been set back up properly
            if simulation.system.time > simulation.stop_times[-1]: break

            return_val = C.run_kmc_spin(stop_time, simulation.system.SD)
            if return_val == -1: 
                print "No more possible moves"
                break
            if options.output_prefix and options.write_trj:
                simulation.write_frame()

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
            predicted_total_walltime = wall_time_remaining+elapsed_time
            total_steps = total_steps+simulation.system.SD.current_step
            efficiency = simulation.system.SD.current_step/stage_elapsed_time

            #this should be last major thing in loop, before printing logging info
            prev_time = sim_timer.gettime()
            prev_stop_time = stop_time

            E_per_site = simulation.system.total_energy/simulation.system.nsites
            Teff = 1/np.log(1/E_per_site-1)
            print "Time: %.2e Dt: %3.2e | Energy: %.2e Teff: %.4e\n\tElapsed: %.2e Pred: %.2e Etr: %.2e Eff: %3.2e"%( simulation.system.time, avg_dt, E_per_site, Teff, elapsed_time, predicted_total_walltime, wall_time_remaining, efficiency )

        # after finishing, write final state
        if options.output_prefix and options.write_trj:
            simulation.save_state(options.output_prefix+'.final.spinsim.gz')
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
    parser.add_option('-r','--restart',default=None,
                      help="Read in stored configuration and do full restart (default: use random)")
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

    freezing_group=OptionGroup(parser,"Options for freezing spins","Options for freezing spins")
    freezing_group.add_option('--save_frozen',dest='discard_frozen',default=True,action="store_false",
                              help="Store frozen spins in trajectory. Always saved in 'start' and 'final' files. (default: false)")
    freezing_group.add_option('--center_coord','--ccoord',dest='ccoord',default=None,type='string',
                              help="Set a position to use for a cavity, wall or sandwich geometry")
    freezing_group.add_option('--cradius','--frozen_radius',dest='frozen_radius',default=None,type=float,
                              help="Set half size of frozen area")
    freezing_group.add_option('--frozen_fraction',default=None,type=float,
                              help="Set fraction of frozen particles in random freezing geometry")
    freezing_group.add_option('--frozen_geometry',default=None,type='string',
                              help="Set a frozen geometry to use. Available options: %s"%sorted( frozen_geometries.keys()))
    freezing_group.add_option('--frozen_dimension',default=0,type=int,
                              help="Set axis perpendicular to frozen slab if WALL or SANDWICH geometry (default: %default)")
    parser.add_option_group(freezing_group)

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
