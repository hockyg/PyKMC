#!/usr/bin/python
import pyximport; pyximport.install()
import os
import sys
from SpinObj import *
import cPickle as pickle

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
    options, args = parser.parse_args()
 
    # add a more complicated option for this later, once deciding on writing out
    info_steps = options.max_steps/10
 
    if options.seed < 1:
        parser.error("Seed must be assigned a positive value")
 
    if options.model in ModelRegistry:
        model = ModelRegistry[options.model]
    else:
        print "Model not yet defined. Please select from:"
        print "\t",ModelRegistry.keys()
        parser.print_help()
        sys.exit(1)
 
    if options.lattice in LatticeRegistry:
        lattice = LatticeRegistry[options.lattice]
    else:
        print "Lattice not yet defined. Please select from:"
        print "\t",LatticeRegistry.keys()
        parser.print_help()
        sys.exit(1)
 
    if options.input is not None:
        print "Reading input not currently supported"
        sys.exit(1)
    else:
        initial_configuration = model.RandomConfiguration( options.nsites, options.temperature )

    nneighbors_per_site, neighbors = lattice.Neighbors(options.nsites)
 
    simulation = Simulation( lattice.lattice_name, model.model_name, options.nsites,
                             options.max_steps, configuration=initial_configuration )
    #results = driveKCM( simulation, info_steps, options.temperature, options.seed )
 
    if options.output_prefix:
        output_pickle = options.output_prefix+'_results.pickle'
        pickle.dump(results,open(output_pickle,'wb'), protocol=-1)

if __name__ == "__main__":
    main()
