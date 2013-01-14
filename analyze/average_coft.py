#!/usr/bin/python
import os
import sys
import gzip
import cPickle as pickle
import numpy as np

from optparse import OptionParser,OptionGroup

usage="%prog [options] spintrj "
parser=OptionParser(usage=usage)
parser.add_option('-o',dest="output_prefix",default=None,help="Prefix for plots and pickle of results")

(options, args) = parser.parse_args()

if options.output_prefix is not None:
    sys.stdout = open(options.output_prefix+'_coft_avg_std.txt','w')

t0 = None
file0 = None
coft_list = []
c0_list = []
ndefect_avg_list = []

for fileidx, picklefile in enumerate(args):
    result_dict = pickle.load(gzip.open(picklefile,'rb'))
    if fileidx == 0:
        t0 = result_dict['time_array']
        file0 = None
    else:
        if not np.all(t0 == result_dict['time_array']):
            print "# Times in file",file0,"do not match those in",picklefile," - Skipping",picklefile
            continue
    coft = result_dict['coft']
    coft_list.append(coft)
    c0_list.append(result_dict['c0'])
    if "defect_array" in result_dict:
        ndefect_avg_list.append( result_dict["defect_array"].mean() )
        nsites = result_dict["nsites"]

coft_array = np.array(coft_list)
means = coft_array.mean(axis=0).reshape(-1,1)
stds = coft_array.std(axis=0).reshape(-1,1)
times = t0.reshape(-1,1)
c0_array = np.array(c0_list)

result_array = np.concatenate(( times, means, stds ), axis=1)
np.savetxt(sys.stdout, result_array, fmt="%f")

print "# c0 mean",c0_array.mean()
#print c0_array

if options.output_prefix is not None:
    #histogram c0 values
    #binstep=0.0666
    binstep=0.05
    bins = np.arange(0,1+binstep,binstep)
    hist,edges = np.histogram(c0_array,bins=bins)
    #hist,edges = np.histogram(c0_array,bins=20)
    #bin_centers = edges[:-1]+binstep/2
    bin_centers = (edges[:-1]+edges[1:])/2
    bin_centers = bin_centers.reshape(-1,1)
    hist = (hist/float(hist.sum())).reshape(-1,1)
    hist_result = np.concatenate(( bin_centers, hist ),axis=1)
    
    np.savetxt(open(options.output_prefix+'_pofq.txt','w'),hist_result, fmt="%f")

    if len(ndefect_avg_list)>0:
       ndefect_avg_array = np.array(ndefect_avg_list).reshape(-1,1)
       fdefect_avg_array = ndefect_avg_array/nsites
       defect_scatter_result = np.concatenate(( ndefect_avg_array, fdefect_avg_array, c0_array.reshape(-1,1) ),axis=1)
       np.savetxt(open(options.output_prefix+'_defect_v_c0.txt','w'), defect_scatter_result, fmt="%f")
       bins = np.array(np.arange(0,ndefect_avg_array.max()+0.5,1),dtype=int).reshape(-1,1)
       ndefect_v_c0_avg = np.zeros(len(bins)).reshape(-1,1)
       ndefect_v_c0_std = np.zeros(len(bins)).reshape(-1,1)
       ndefect_v_c0_error = np.zeros(len(bins)).reshape(-1,1)
       ndefect_v_c0_nentries = np.zeros(len(bins)).reshape(-1,1)
       ndefect_avg_bins = np.floor(ndefect_avg_array)
       
       for i in bins:
           tmp_array = c0_array[ np.where(ndefect_avg_bins==i)[0] ] 
           nitems = len(tmp_array)
           if nitems>0:
               ndefect_v_c0_avg[i] = tmp_array.mean()
               ndefect_v_c0_std[i] = tmp_array.std()
               ndefect_v_c0_error[i] = ndefect_v_c0_std[i]/np.sqrt(nitems)
               ndefect_v_c0_nentries[i] = nitems/float(len(ndefect_avg_bins))
       defect_bin_result = np.concatenate(( bins, ndefect_v_c0_avg, ndefect_v_c0_std, ndefect_v_c0_error, ndefect_v_c0_nentries ),axis=1)
       np.savetxt(open(options.output_prefix+'_ndefect_bin_v_c0.txt','w'), defect_bin_result, fmt="%f")
       

