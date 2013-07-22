#!/usr/bin/python

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import math
import scipy.stats
import scipy
import glob
import re
import sys
import getopt
import argparse
from collections import defaultdict
from scipy import optimize

s  = 2
dc = 0

#fitfunc = lambda c, x: -1 / (c[0] + c[1]*(12.0/x) + c[2]*(12.0/x)**2 + c[3]*(12.0/x)**3 - x/12.0)
fitfunc = lambda c, x: c[0] + c[1]*x + c[2]*x**2 + c[3]*x**3
errfunc = lambda c, x, y, err: (fitfunc(c, x) - y) / err
c_in = [-0.01, -0.01, 10., 10.]   # Order-of-magnitude initial guesses

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
  tmpbetal = []

  filelist = glob.glob('12flav_'+volume+'/wflow/Wflow_'+'*'+volume+'*')
  for filename in filelist:
    tmparry = re.split('_', filename)
    beta = tmparry[4]
    tmpbetal.append(beta)
    ftext = [i for i in open(filename) if i.startswith('WFLOW ' + str(t))==True][0]
    tmparry = re.split('\s+',ftext)
    data[beta].append(float(tmparry[4]))

  betal[volume]=(set(tmpbetal))
  for b in betal[volume]:
    t2E = np.mean(data[b])
    t2E_err = scipy.stats.sem(data[b])
    g2[(volume,b)] = (128*math.pi**2*t2E)/(3*(3**2-1)*(1+dc))
    g2_err[(volume,b)] = (128*math.pi**2*t2E_err)/(3*(3**2-1)*(1+dc))

  x = []
  y = []
  e = []
  for b in sorted(betal[volume]):
    x.append(float(b))
    y.append(float(g2[(volume,b)]))
    e.append(float(g2_err[(volume,b)]))
  ax = np.array(x)
  ay = np.array(y)
  ae = np.array(e)
  coeff, success = optimize.leastsq(errfunc, c_in[:], args=(ax, ay, ae))
  plt.errorbar(ax, ay, yerr=ae, linestyle='None', marker='.', label=volume)
  plt.plot(ax, fitfunc(coeff,ax))
plt.legend()
plt.savefig("plots/wflow.png", format='png')


