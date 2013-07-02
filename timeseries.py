#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import glob
import re
import sys
import getopt
import argparse
from collections import defaultdict

numbers = re.compile(r'(\d+)')
def numericalSort(value):
    parts = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])
    return parts

parser = argparse.ArgumentParser(description='Produce timeseries plots of blocked observables.')
parser.add_argument('volume', metavar='volume', type=int, nargs=1, help='Enter the volume you are intersted in.')
parser.add_argument('tag', metavar='tag', type=str, nargs=1, help='Enter the base name of the run you are intersted in.')
parser.add_argument('tsmear', metavar='tsmear', type=float, nargs=1, help='Enter the Wflow smearing time you are interested in.')
args = parser.parse_args()
argdict = vars(args)
volume = argdict['volume'].pop()
tag = argdict['tag'].pop()
tsmear = argdict['tsmear'].pop()

value = defaultdict(list)
search = '12flav_'+str(volume)+'/'+tag
filelist = sorted(glob.glob(search+'.*'), key=numericalSort)

for filename in filelist:
  f = open(filename, 'r')
  ftext = f.readlines()
  for line in ftext:
    if re.match('^LOOPS.*', line):
      line=line.strip()
      junk, smear_t, loop, junk, blk, junk, val = re.split('\s+', line)
      if float(smear_t) == tsmear:
        if loop == '0':
          value[blk].append(float(val))

for key in value:
  plt.title('Blocking step '+key)
  plt.plot(value[key])
  plt.show()
