#!/usr/bin/python
import pyximport; pyximport.install()
import os
import sys
import cPickle as pickle
import gzip
from SpinObj import *
from PyKMC.analyze import util

maxtouse = 500
maxframe=1000
nsamples = len(sys.argv[1:])
coft_list = []
times = []
time_array = None
for fileidx, filename in enumerate(sys.argv[1:]):
    times, trajectory = util.get_spintrj(filename,maxframe=maxframe)
    frame_intervals = util.logframes(len(times))
    nframes = trajectory.shape[0]
    dt = times[1]-times[0]
    counts = np.zeros(len(frame_intervals))
    coft_array = np.zeros((len(frame_intervals),maxtouse))
    time_array = frame_intervals*dt

    means = np.zeros(len(frame_intervals))
    for idx0,frame_interval in enumerate(frame_intervals):
        cfg0 = trajectory[idx0,:] 
        starting_frames = util.get_starting_frames( nframes, frame_interval, maxtouse=maxtouse )
        for idx1,starting_frame in enumerate(starting_frames):
            if counts[idx0] >= maxtouse: continue
            cfg1 = trajectory[idx1,:] 
            coft = (cfg0*cfg1).mean()
            coft_array[idx0,counts[idx0]] = coft
            counts[idx0]+=1
        means[idx0] = coft_array[idx0,:counts[idx0]].mean()

    coft_list.append(means)

coft_all = np.array(coft_list)
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
