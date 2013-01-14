#!/usr/bin/python
import numpy as np
import pyximport;pyximport.install(setup_args={'include_dirs': np.get_include()})
import os
import sys
import cPickle as pickle
import gzip
from PyKMC.SpinObj import *
from PyKMC.analyze import util
from optparse import OptionParser,OptionGroup
import time

usage="%prog [options] spintrj "
parser=OptionParser(usage=usage)
parser.add_option('-o',dest="output_prefix",default=None,help="Prefix for plots and pickle of results")
parser.add_option('-m','--maxtouse',default=200,type=int,help="Max frame intervals to use (default %default)")
#parser.add_option('--first_frame',default=None,type=int,help="First frame in trajectory to use (default: 0)")
parser.add_option('--last_frame',default=None,type=int,help="Last frame in trajectory to use (default: all)")

(options, args) = parser.parse_args()

if options.output_prefix is not None:
    sys.stdout=open(options.output_prefix+'_coft.txt', 'w')

maxtouse = options.maxtouse
maxframe=options.last_frame
nsamples = len(args)
coft_list = []
times = []
time_array = None
c0_list = []
defect_list = []
for fileidx, filename in enumerate(args):
    try:
        starttime = time.time()
        times, stop_times, trajectory, dual_trajectory = util.get_spintrj(filename,maxframe=maxframe,plus_dual=True)
#        tmp_times, tmp_stop_times, dual_trajectory = util.get_spintrj(filename,maxframe=maxframe,cfgname="dual_configuration")
        print "# Read time:",time.time()-starttime
    except IOError:
        continue
   
    # get number of defects in each frame and put in defect list 
    # also time average site values
    nsites=len(dual_trajectory[0])
    nframes = len(dual_trajectory)
    defect_array = np.zeros(nframes)
    site_values = np.zeros(nsites)
    for i in range(nframes):
        defect_array[i] = dual_trajectory[i].sum()
        site_values = site_values + trajectory[i]
    defect_list.append(defect_array)
    site_values = site_values/float(nframes)

    frame_intervals = util.logframes(len(times))
    nframes = len(trajectory)
    dt = stop_times[1]-stop_times[0]
    counts = np.zeros(len(frame_intervals))
    coft_array = np.zeros((len(frame_intervals),maxtouse))
    time_array = frame_intervals*dt

    # time average spin values
#    site_values = trajectory.mean(axis=0)
    c0 = (site_values*site_values).mean()
    c0_list.append(c0)
    print "# %s: c0 = %f"%(filename,c0)

    means = np.zeros(len(frame_intervals))
    for idx0,frame_interval in enumerate(frame_intervals):
        starting_frames = util.get_starting_frames( nframes, frame_interval, maxtouse=maxtouse )
        for starting_frame in starting_frames:
            if counts[idx0] >= maxtouse: continue
            cfg0 = trajectory[starting_frame] 
            cfg1 = trajectory[starting_frame+frame_interval] 
            coft = (cfg0*cfg1).mean()
            coft_array[idx0,counts[idx0]] = coft
            counts[idx0]+=1
        means[idx0] = coft_array[idx0,:counts[idx0]].mean()

    coft_list.append(means)

defect_array = np.array(defect_list)
c0_array = np.array(c0_list)
coft_all = np.array(coft_list)
print "#Avg: %f"%(c0_array.mean())

#longvalue = 

#print coft_all.shape
final_avg = coft_all.mean(axis=0)
output = np.append(time_array.reshape(-1,1),final_avg.reshape(-1,1),axis=1)
np.savetxt(sys.stdout,output,fmt="%f")

# now save important information as dict
if options.output_prefix:
    results_dict = { 'time_array': time_array, 'coft': final_avg, 'c0': c0_array.mean(),'defect_array':defect_array,'nsites':nsites }
    pickle.dump(results_dict, gzip.GzipFile( options.output_prefix+'_coft.pickle.gz', 'wb'),protocol=-1)

print "# Total time:",time.time()-starttime

#    if fileidx == 0:
#        coft_array = np.zeros((len(times),nsamples))
#        coft_array[:,0] = np.array(cofts)

#coft_means = coft_array.mean(axis=1)
#coft_vars = coft_array.var(axis=1)
#for i in range(len(times)):
#    print times[i],coft_means[i],coft_vars[i]*configuration.size
