#!/usr/bin/python
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import pandas as pd
from scipy import interpolate
import numpy as np
import math
import glob
import re
import sys
import argparse
from collections import defaultdict

numbers = re.compile(r'(\d+)')
def numericalSort(value):           #sorting routine to sort files read in from directory
    parts = numbers.split(value)    #to get correct timeseries
    parts[1::2] = map(int, parts[1::2])
    return parts

#Just a parser takes volume and beta as arguments
parser = argparse.ArgumentParser(description='Wilson Flow Auto Correlation')
parser.add_argument('v', metavar='volume', type=str,nargs=1, help='Enter the volume you would like to find the Auto Correlation time for')
parser.add_argument('b', metavar='beta', type=str,nargs=1, help='Enter the beta you would like to find the Auto Correlation time for')
args    = parser.parse_args()
argdict = vars(args)
v       = argdict['v'].pop()
b       = argdict['b'].pop()
flav    = '12'
c       = 0.30 #common values used in analysis
L       = int(v[:len(v)/2]) #get L from volume i.e. 12 from 1212
data    = []

#rips through files, interpolates data in files, generates timeseries data np array
filelist = sorted(glob.glob(flav+'flav_'+v+'/wflow/dat/wflow_'+v+'_'+b+'_0.0.*'), key=numericalSort)
for filename in filelist:
  frame = pd.read_table(filename, sep=' ')
  f = interpolate.interp1d(frame['t'],frame['E'])
  t = ((c*L)**2)/8                #find the amount of wilson flow
  data.append(float(f(t)))     #interpolate data to the desired wilson flow
npdata  = np.array(data)       #stores data in pandas series
print "Length of timeseries:  " + str(npdata.size)

#calculates the jack knife error estimate over a range of block sizes
bins = [(2*x) for x in range(1,(npdata.size//10))]
jack_err = []
for binsize in (bins):
  binmeans = npdata[:(npdata.size // binsize) * binsize].reshape(-1, binsize).mean(axis=1) 
  total = np.sum(binmeans)
  avg   = np.mean(binmeans)
  n     = binmeans.size - 1
  accum = 0
  for x in np.nditer(binmeans):
    jackavg = (total - x) / n
    diff    = jackavg-avg
    accum   += math.pow(diff,2)
  N = float(binmeans.size)
  jack_err.append(math.sqrt((N-1)/N*accum))

fig1  = plt.figure(1)
stdev = [np.std(npdata)/math.sqrt(npdata.size)]*len(bins)
x     = np.array(bins)
y     = np.array(jack_err)
y2    = np.array(stdev)
y3    = y2/(math.sqrt(npdata.size))
splt1 = fig1.add_subplot(111)
splt1.set_title("Jackknife analysis of volume " + v + " beta " + b)
splt1.set_xlabel('block size')
splt1.set_xscale('log')
splt1.ticklabel_format(axis='y',style='sci',scilimits=(-2,2))
splt1.set_ylabel('error estimate on <t^2*E>')
splt1.plot(x, y, linestyle='None', marker='.', color='red')
splt1.plot(x, y2, color='blue',label='stdev')
splt1.plot(x, y3, color='green',label='sdev_mean')
splt1.legend(loc=2)
fig1.savefig("plots/jack_"+v+"_"+b+".png", format='png')
