#!/usr/bin/python

import pandas as pd
from pandas import Series, DataFrame, MultiIndex
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import re
import sys
import getopt
import random

numbers = re.compile(r'(\d+)')
def numericalSort(value):
    parts = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])
    return parts

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
    print 'Reading filename:  '+filename
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
  data = Series(value, index=[smear, block, obser])
  data_list = data.ix[('0.07','4','4')].values
  autoc = np.correlate(data_list, data_list, mode='full')
  iact = sum(autoc[autoc.size/2:]/autoc[autoc.size/2])


if __name__ == "__main__":
   main(sys.argv[1:])
