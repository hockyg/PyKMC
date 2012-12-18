#!/usr/bin/python
import pyximport; pyximport.install()
import os
import sys
import cPickle as pickle
import gzip
from SpinObj import *
from PyKMC.analyze import util

maxtouse = 200
maxframe=None
#maxframe=100
nsamples = len(sys.argv[1:])
coft_list = []
times = []
time_array = None
c0_list = []
for fileidx, filename in enumerate(sys.argv[1:]):
    try:
        times, stop_times, trajectory = util.get_spintrj(filename,maxframe=maxframe)
    except IOError:
        continue
    frame_intervals = util.logframes(len(times))
    nframes = trajectory.shape[0]
    dt = stop_times[1]-stop_times[0]
    counts = np.zeros(len(frame_intervals))
    coft_array = np.zeros((len(frame_intervals),maxtouse))
    time_array = frame_intervals*dt

    # time average spin values
    site_values = trajectory.mean(axis=0)
    c0 = (site_values*site_values).mean()
    c0_list.append(c0)
    print "# %s: c0 = %f"%(filename,c0)

    means = np.zeros(len(frame_intervals))
    for idx0,frame_interval in enumerate(frame_intervals):
        starting_frames = util.get_starting_frames( nframes, frame_interval, maxtouse=maxtouse )
        for starting_frame in starting_frames:
            if counts[idx0] >= maxtouse: continue
            cfg0 = trajectory[starting_frame,:] 
            cfg1 = trajectory[starting_frame+frame_interval,:] 
            coft = (cfg0*cfg1).mean()
            coft_array[idx0,counts[idx0]] = coft
            counts[idx0]+=1
        means[idx0] = coft_array[idx0,:counts[idx0]].mean()

    coft_list.append(means)

c0_array = np.array(c0_list)
coft_all = np.array(coft_list)
print "#Avg: %f"%(c0_array.mean())

#longvalue = 

#print coft_all.shape
final_avg = coft_all.mean(axis=0)
output = np.append(time_array.reshape(-1,1),final_avg.reshape(-1,1),axis=1)
np.savetxt(sys.stdout,output,fmt="%f")


#    if fileidx == 0:
#        coft_array = np.zeros((len(times),nsamples))
#        coft_array[:,0] = np.array(cofts)

#coft_means = coft_array.mean(axis=1)
#coft_vars = coft_array.var(axis=1)
#for i in range(len(times)):
#    print times[i],coft_means[i],coft_vars[i]*configuration.size
