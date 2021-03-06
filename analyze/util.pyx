#!/usr/bin/python
import pyximport; pyximport.install()
import os
import sys
import cPickle as pickle
import gzip
from PyKMC.SpinObj import *

def get_starting_frames( nframes,frame_interval,maxtouse=None):
    if maxtouse is None:
        return range( nframes-frame_interval )
    else:
        #note, next line is deliberately integer math
        if maxtouse>nframes:
           step=1
        else:
           step=int(nframes)/int(maxtouse)
        return range(0,nframes-frame_interval,step)

def logframes(total_steps,steps_per_decade=10):
    linearsteps=np.arange(1,total_steps+1,1.0)
    logtimepoints=np.log10(linearsteps)
    minlogtimepoint=np.min( logtimepoints )
    maxlogtimepoint=np.max( logtimepoints )
    evenlogspace=np.arange(minlogtimepoint,maxlogtimepoint,1.0/steps_per_decade)
    logdist=np.array( np.unique( np.around( np.power(10,evenlogspace) ) )  ,dtype='int')

    return logdist

def get_last_cfg( spintrjfile ):
    fh = gzip.GzipFile(spintrjfile,'r')

    try:
        count = 0
        while True:
            simdata = pickle.load(fh)
            configuration = simdata["configuration"]
            dual_configuration = simdata["dual_configuration"]

            count = count+1
    except EOFError:
        pass
    except IOError:
        raise
    return configuration, dual_configuration

def get_spintrj( spintrjfile, cfgname="configuration",maxframe=None, plus_dual=False):
    fh = gzip.GzipFile(spintrjfile,'r')
    times = []
    stop_times = []
    trajectory = []
    dual_trajectory = []

    try:
        count = 0
        while True:
            if (maxframe is not None and count>maxframe): break
            simdata = pickle.load(fh)
            sim_time = simdata["time"]
            stop_time = simdata["stop_time"]
            configuration = simdata[cfgname]
            times.append(sim_time)
            stop_times.append(stop_time)
            trajectory.append(configuration)
            if plus_dual is True:
                configuration = simdata["dual_configuration"]
                dual_trajectory.append(configuration)

            count = count+1
    except EOFError:
        pass
    except IOError:
        raise
    #trajectory_array = np.array(trajectory,dtype=int)
    times_array = np.array(times)
    stop_times_array = np.array(stop_times)

    if plus_dual is True:
        return times_array, stop_times_array, trajectory, dual_trajectory
    else:
        return times_array, stop_times_array, trajectory
    #return times_array, stop_times_array, trajectory_array

if __name__ == "__main__":
    get_spintrj(sys.argv[1])
