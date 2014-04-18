#!/usr/bin/python

import os.path
import csv
import numpy as np
import pandas as pd
import glob
import re
import argparse
import math
from collections import defaultdict

#sorting routine to sort files read in from the directory to get the correct timeseries
numbers = re.compile(r'(\d+)')
def numericalSort(value):
  parts = numbers.split(value)
  parts[1::2] = map(int, parts[1::2])
  return parts

#checks that the smearings are the same for all configurations, returns a list of smearing values
#right now not 100% accomplishing that goal
def checkT(d):
  #same_smearing = len(set(d.values()))==1
  #print same_smearing
  for key in d:
    uniq_smearing = d[key]
    next
  return(uniq_smearing)

#clears dictionaries
def initDict(T,E,t2E,dt2E):
  T.clear()
  E.clear()
  t2E.clear()
  dt2E.clear()

#average calculation
def Avg(a):                                                     #argument is list of numbers
  npdata = np.array(a)
  mean = np.mean(npdata)
  return mean

#jackknife error calculation
def jackErr(a,b):                                               #arguments are list of nubmers and beta
  npdata  = np.array(a)
  binsize = int(bin_df.ix[(int(vol),float(b))]['bin'])
  bindata = npdata[:(npdata.size // binsize) * binsize].reshape(-1, binsize).mean(axis=1) 
  avg     = np.mean(bindata)
  total   = np.sum(bindata)
  n       = bindata.size - 1
  accum   = 0

  for x in np.nditer(bindata):
    jackavg = (total - x) / n
    diff    = jackavg-avg
    accum   += math.pow(diff,2)

  N       = float(bindata.size)
  err     = math.sqrt((N-1)/N*accum)
  return(err)

#calculates average, error, and prints to file
def Crunch(b):                                                  #argument is beta
  print "Calculating statistics on volume: %s, beta: %s" % (vol,b)
  smear_t = checkT(T)
  avg_E, avg_t2E, avg_dt2E = [], [], []
  err_E, err_t2E, err_dt2E = [], [], []

  for tau in smear_t:
    avg_E.append(Avg(E[tau]))
    avg_t2E.append(Avg(t2E[tau]))
    avg_dt2E.append(Avg(dt2E[tau]))

    err_E.append(jackErr(E[tau],b))
    err_t2E.append(jackErr(t2E[tau],b))
    err_dt2E.append(jackErr(dt2E[tau],b))

  fileout = '12flav_'+vol+'/wflow/avg/avgwflow_'+vol+'_'+b+'_0.0.csv'
  rows = zip(smear_t, avg_E, err_E, avg_t2E, err_t2E, avg_dt2E, err_dt2E)
  with open(fileout, 'wb') as f:
    print "Writing summary file: %s\n" % fileout
    writer = csv.writer(f)
    writer.writerows(rows)

  initDict(T,E,t2E,dt2E)

parser = argparse.ArgumentParser(description='Print averages and errors for each beta of a given volume for every t')
parser.add_argument('v', metavar='volume', type=str, nargs=1, help='Enter the volume.')
args = parser.parse_args()
argdict = vars(args)
vol = argdict['v'].pop()

bin_df = pd.read_csv("binsize.csv", index_col=[0,1], skipinitialspace=True)

T        = defaultdict(list)                                    #initialize dicts
E        = defaultdict(list)
t2E      = defaultdict(list)
dt2E     = defaultdict(list)
old_beta = 0

filelist = sorted(glob.glob('12flav_'+vol+'/wflow/dat/*'), key=numericalSort)
for filename in filelist:
  mont_t  = filename.rsplit('.',1)[1]                           #parse the file name
  tmparry = re.split('_', filename)
  beta    = tmparry[3]
  
  if (beta != old_beta) and (old_beta != 0):                    #if we have a new beta caluclate!
    Crunch(old_beta)

  with open(filename,'rb') as csvfile:                          #read file
    reader = csv.reader(csvfile, delimiter=' ')
    reader.next()
    print "Parsing: %s" % filename
    for line in reader:
      tau = float(line[0])                                      #smearing
      T[mont_t].append(tau)                                     #dict list with smearing
      E[tau].append(float(line[2]))                             #dict list with E
      t2E[tau].append(float(line[3]))                           #dict list with t2E
      dt2E[tau].append(float(line[4]))                          #dict list with dt2E

  old_beta = beta

Crunch(beta)                                                    #have to calculate last beta 
