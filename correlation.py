#!/usr/bin/python

import pandas as pd
from pandas import Series, DataFrame, MultiIndex
import numpy as np
import matplotlib.pyplot as plt
import glob
import re
import sys
import getopt

numbers = re.compile(r'(\d+)')
def numericalSort(value):
    parts = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])
    return parts

def uniq(seq, idfun=None):  
    # order preserving 
    if idfun is None: 
        def idfun(x): return x 
    seen = {} 
    result = [] 
    for item in seq: 
        marker = idfun(item) 
        # in old Python versions: 
        # if seen.has_key(marker) 
        # but in new ones: 
        if marker in seen: continue 
        seen[marker] = 1 
        result.append(item) 
    return result

def main(argv):
  flav = ''
  tag = ''
  try:
    opts, args = getopt.getopt(argv,"hf:t:",["flav=","tag="])
  except getopt.GetoptError:
    print 'correlation.py <flavour> <tag>'
    sys.exit(2)
  for opt, arg in opts:
    if opt == '-h':
       print 'correlation.py <flavour> <tag>'
       sys.exit()
    elif opt in ("-f", "--flavour"):
      flav = arg
    elif opt in ("-t", "--tag"):
      tag = arg

  parse_tag = re.split('_', tag)
  vol = parse_tag[2]
  smear = []
  block = []
  obser = []
  value = []
  
  filelist = sorted(glob.glob(flav+'flav_'+vol+'/'+tag+'.*'), key=numericalSort)
  for filename in filelist:
    f = open(filename, 'r')
    ftext = f.readlines()
    for line in ftext:
      if re.match('^LOOPS.*', line):
        line=line.strip()
        junk, smear_t, loop, junk, blk, junk, val = re.split('\s+', line)
        smear.append(smear_t)
        block.append(blk)
        obser.append(loop)
        value.append(float(val))
  zipped = zip(smear, block, obser)
  ind = MultiIndex.from_tuples(zipped, names=['smearing', 'block_level', 'observable'])
  data = Series(value, index=ind)
  iact = []
  uni_zipped = uniq(zipped)
  for i in uni_zipped:
    timeseries = data.ix[i].values
    mean = np.mean(timeseries)
    timeseries = np.array([x - mean for x in timeseries])
    autocorr_f = np.correlate(timeseries, timeseries, mode='full')
    iact.append(sum(autocorr_f[autocorr_f.size/2:]/autocorr_f[autocorr_f.size/2]))
    #temp = autoc[autoc.size/2:]/autoc[autoc.size/2]
    #plt.plot(temp)
    #plt.show()
  ind = MultiIndex.from_tuples(uni_zipped, names=['smearing', 'block_level', 'observable'])
  result = Series(iact, index=ind)
  print result

if __name__ == "__main__":
   main(sys.argv[1:])
