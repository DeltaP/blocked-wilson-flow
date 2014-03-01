#!/usr/bin/python

import sys
import csv
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import math
import scipy.stats
import scipy
import glob
import re
import argparse
import time
import datetime
from collections import defaultdict
from scipy.optimize import curve_fit
from scipy import optimize 

parser = argparse.ArgumentParser(description='Wilson Flow Crossing.')
parser.add_argument('c', metavar='constant', type=float, nargs=1, help='Enter the smearing constant that defines the extent of the smearing you are interested in.')
parser.add_argument('v', metavar='volume', type=str, nargs=1, help='Enter the volume.')
args = parser.parse_args()
argdict = vars(args)
c    = argdict['c'].pop()
vol = argdict['v'].pop()

d_c = {
    0.2  : -0.005,
    0.25 : -0.012,
    0.3  : -0.03,
    0.35 : -0.045,
    0.4  : -0.07,
    }
dc = d_c[c]


bin_df = pd.read_csv("binsize.csv", index_col=[0,1], skipinitialspace=True)

t_start = time.time()
print "Working on %s:" % vol
L = int(vol[:len(vol)/2])
t = ((c*L)**2)/8
data     = defaultdict(list)
tmparry  = []
tmpbetal = []

filelist = glob.glob('12flav_'+vol+'/wflow/dat/*')
for filename in filelist:
  tmparry = re.split('_', filename)
  beta = tmparry[3]
  tmpbetal.append(beta)
  frame = pd.read_table(filename, sep=' ')
  f = scipy.interpolate.interp1d(frame['t'],frame['E'])
  data[beta].append(f(t))

Xval = []
Yval = []
Eval = []
betal=(set(tmpbetal))
for b in sorted(betal):
  npdata = np.array(data[b])
  binsize = int(bin_df.ix[(int(vol),float(b))]['bin'])
  bindata = npdata[:(npdata.size // binsize) * binsize].reshape(-1, binsize).mean(axis=1) 
  E     = np.mean(bindata)
  total = np.sum(bindata)
  n     = bindata.size - 1
  accum = 0
  for x in np.nditer(bindata):
    jackavg = (total - x) / n
    diff    = jackavg-E
    accum   += math.pow(diff,2)
  N     = float(bindata.size)
  E_err = math.sqrt((N-1)/N*accum)
  g2    = (128*math.pi**2*E*t**2)/(3*(3**2-1)*(1+dc))
  g2_err= (128*math.pi**2*E_err*t**2)/(3*(3**2-1)*(1+dc))

  Xval.append(b)
  Yval.append(g2)
  Eval.append(g2_err)

rows = zip(Xval,Yval,Eval)
fname = str(vol) + "_" + str(c) + ".csv"
with open(fname, 'wb') as f:
  writer = csv.writer(f)
  writer.writerows(rows)
