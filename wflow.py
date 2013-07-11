#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import math
import glob
import re
import sys
import getopt
import argparse
from collections import defaultdict

s  = 2
dc = 0

parser = argparse.ArgumentParser(description='Wilson Flow Matching.')
parser.add_argument('c', metavar='constant', type=float, nargs=1, help='Enter the smearing constant that defines the extent of the smearing you are interested in.')
parser.add_argument('volumes', metavar='volume', type=str, nargs= '+', help='Enter the volume(s) you are intersted in.')
args = parser.parse_args()
argdict = vars(args)
c = argdict['c'].pop()
vol = defaultdict(list)
vol = argdict['volumes']

for volume in vol:
  L = int(volume[:len(volume)/2])
  t = (c*L)**2/8
  betal    = defaultdict(list)
  data     = defaultdict(list)
  g2       = defaultdict(dict)
  tmparry  = []

  filelist = glob.glob('12flav_'+volume+'/Wflow_'+'*'+volume+'*')
  for filename in filelist:
    tmparry = re.split('_', filename)
    beta = tmparry[4]
    betal[volume].append(beta)
    ftext = [i for i in open(filename) if i.startswith('WFLOW ' + str(t))==True][0]
    tmparry = re.split('\s+',ftext)
    data[beta].append(float(tmparry[4]))

  for b in betal[volume]:
    t2E = np.mean(data[b])
    g2[(volume,b)] = (128*math.pi**2*t2E)/(3*(3-1)*(1+dc))
  x = []
  y = []
  for b in betal[volume]:
    x.append(b)
    y.append(g2[(volume,b)])
  ax = np.array(x)
  ay = np.array(y)
  plt.plot(ax, ay,linestyle='None', marker='.')
plt.show()
