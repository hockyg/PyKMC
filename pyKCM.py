#!/usr/bin/python
import os
import sys
from driver import *

def main():
   from optparse import OptionParser, OptionGroup
   parser = OptionParser()
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
   parser.add_option('--seed', default=0, type=int,
                     help="Random seed for simulation (default: %default)" )
   options, args = parser.parse_args()

   # add a more complicated option for this later, once deciding on writing out
   info_steps = 10

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

   # later, allow reading in of the configuration here
   simulation = Simulation( lattice, model, options.nsites,
                            options.max_steps, configuration=None )
   # if options.mode == "KCM" 
   driveKCM( simulation, info_steps, options.temperature, options.seed )
   # create simulation object
   # pass simulation object to driver

if __name__ == "__main__":
    main()
