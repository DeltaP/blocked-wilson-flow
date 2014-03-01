#!/usr/bin/python

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from collections import defaultdict
from scipy import interpolate
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob
import re
import sys
import argparse

numbers = re.compile(r'(\d+)')
def numericalSort(value):           #sorting routine to sort files read in from directory
    parts = numbers.split(value)    #to get correct timeseries
    parts[1::2] = map(int, parts[1::2])
    return parts

#Just a parser takes volume and beta as arguments
#v     = '1212'
#beta  = np.arange(3.4,8.0,0.2)
#v     = '1616'
#beta  = np.arange(3.0,8.0,0.2)
v     = '1818'
beta  = np.arange(3.0,8.5,0.5)
#v     = '2424'
#beta  = (4.0, 4.25, 4.5, 4.75, 5.0, 5.25, 5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 7.0, 7.5, 8.0)
#v     = '3232'
#beta  = (4.0, 4.5, 4.75, 5.0, 5.25, 5.5, 5.75, 6.0, 6.5, 7.0, 8.0)
#v     = '3636'
#beta  = (4.5, 5.0, 5.5, 6.0, 6.5, 8.0)
flav  = '12'
c     = 0.30                        #chose 0.3 because it has longest autocorr
L     = int(v[:len(v)/2])           #get L from volume i.e. 12 from 1212
data  = defaultdict(list)

#rips through files, interpolates data in files, generates timeseries data structure
for b in beta:
  filelist = sorted(glob.glob(flav+'flav_'+v+'/wflow/dat/wflow_'+v+'_'+str(b)+'_0.0.*'), key=numericalSort)
  for filename in filelist:
    frame = pd.read_table(filename, sep=' ')
    f = interpolate.interp1d(frame['t'],frame['E'])
    t = ((c*L)**2)/8                #find the amount of wilson flow
    data[b].append(float(f(t)))     #interpolate data to the desired wilson flow
frame_data  = pd.Series(data)       #stores data in pandas series

#calculates the autocorrelation based on each timeseries for each c
iacf = []
for i in beta:
  timeseries = frame_data[i]
  timeseries -= np.mean(timeseries) 
  #calculates autocorrelation function
  autocorr_f = np.correlate(timeseries, timeseries, mode='full')
  #sums the forward autocorrelation function until it reaches zero
  iacf.append(np.max(np.cumsum(autocorr_f[(autocorr_f.size/2.0):]/autocorr_f[(autocorr_f.size/2.0)])))

fig1=plt.figure(1)
x = np.array(beta)
y = np.array(iacf)
splt1 = fig1.add_subplot(111)
splt1.set_title("Volume " + v + " IACT")
splt1.set_xlabel('beta')
splt1.set_ylabel('Integrated Auto Correlation Time')
splt1.plot(x, y, linestyle='None', marker='.', color='red', label=v)
fig1.savefig("plots/IACT_"+v+".png", format='png')
