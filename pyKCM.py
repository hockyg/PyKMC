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
   parser.add_option('-N', '--nsites', default=100, type=int,
                     help="Number of sites to simulate (default: %default)" )
   options, args = parser.parse_args()

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

   # create simulation object
   # pass simulation object to driver

if __name__ == "__main__":
    main()
