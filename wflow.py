#!/usr/bin/python

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import math
import scipy.stats
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
  t = round(((c*L)**2/8)/(0.02))*0.02
  betal    = defaultdict(list)
  data     = defaultdict(list)
  g2       = defaultdict(dict)
  g2_err   = defaultdict(dict)
  tmparry  = []

  filelist = glob.glob('12flav_'+volume+'/wflow/Wflow_'+'*'+volume+'*')
  for filename in filelist:
    tmparry = re.split('_', filename)
    beta = tmparry[4]
    betal[volume].append(beta)
    ftext = [i for i in open(filename) if i.startswith('WFLOW ' + str(t))==True][0]
    tmparry = re.split('\s+',ftext)
    data[beta].append(float(tmparry[4]))

  for b in betal[volume]:
    t2E = np.mean(data[b])
    t2E_err = scipy.stats.sem(data[b])
    g2[(volume,b)] = (128*math.pi**2*t2E)/(3*(3**2-1)*(1+dc))
    g2_err[(volume,b)] = (128*math.pi**2*t2E_err)/(3*(3**2-1)*(1+dc))

  x = []
  y = []
  e = []
  for b in betal[volume]:
    x.append(b)
    y.append(g2[(volume,b)])
    e.append(g2_err[(volume,b)])
  ax = np.array(x)
  ay = np.array(y)
  ae = np.array(e)
  plt.errorbar(ax, ay, yerr=ae, linestyle='None', marker='.', label=volume)
plt.legend()
plt.savefig("wflow.png", format='png')
