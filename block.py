#!/usr/bin/python

import pandas as pd
from scipy import interpolate
import numpy as np
import matplotlib.pyplot as plt
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

#rips through files, interpolates data in files, generates timeseries data structure
filelist = sorted(glob.glob(flav+'flav_'+v+'/wflow/dat/wflow_'+v+'_'+b+'_0.0.*'), key=numericalSort)
for filename in filelist:
  frame = pd.read_table(filename, sep=' ')
  f = interpolate.interp1d(frame['t'],frame['E'])
  t = ((c*L)**2)/8                #find the amount of wilson flow
  data.append(float(f(t)))     #interpolate data to the desired wilson flow
npdata  = np.array(data)       #stores data in pandas series
print "Length of timeseries:  " + str(npdata.size)

#calculates the autocorrelation based on each timeseries for each c
maxbin = int(math.log(npdata.size,2))-1
bins = [math.pow(2,x) for x in range(0,maxbin)]
print bins
for binsize in (bins):
  print binsize
  binmeans = npdata[:(npdata.size // binsize) * binsize].reshape(-1, binsize).mean(axis=1) 
  print type(binmeans)
  jack_blk = 
