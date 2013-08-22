#!/usr/bin/python

import sys
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
from collections import defaultdict
from scipy import optimize
import time
import datetime

#fitfunc = lambda c, x: -1 / (c[0] + c[1]*(12.0/x) + c[2]*(12.0/x)**2 + c[3]*(12.0/x)**3 - x/12.0)
fitfunc = lambda c, x: c[3] + c[2]*x + c[1]*x**2 + c[0]*x**3
errfunc = lambda c, x, y, err: (fitfunc(c, x) - y) / err
c_in = [1.0, 1.0, 1.0, 1.0]   # Order-of-magnitude initial guesses

parser = argparse.ArgumentParser(description='Wilson Flow Crossing.')
parser.add_argument('s', metavar='scale factor', type=float, nargs=1, help='Enter the scale factor between the two volumes.')
parser.add_argument('c', metavar='constant', type=float, nargs=1, help='Enter the smearing constant that defines the extent of the smearing you are interested in.')
parser.add_argument('t0', metavar='offset', type=float, nargs=1, help='Enter the smearing offset t_0.')
parser.add_argument('sv', metavar='small volume', type=str, nargs=1, help='Enter the small volume.')
parser.add_argument('lv', metavar='large volume', type=str, nargs=1, help='Enter the large volume.')
args = parser.parse_args()
argdict = vars(args)
s    = argdict['s'].pop()
c    = argdict['c'].pop()
t0   = argdict['t0'].pop()
svol = argdict['sv'].pop()
lvol = argdict['lv'].pop()

d_c = {
    0.2  : -0.005,
    0.25 : -0.012,
    0.3  : -0.03,
    0.35 : -0.045,
    0.4  : -0.07,
    }
dc = d_c[c]

betal    = defaultdict(list)
g2       = defaultdict(dict)
g2_err   = defaultdict(dict)
dg2      = defaultdict(dict)
coeff    = defaultdict(dict)
pcolor   = {'66': 'y', '88': 'b', '1212': 'g', '1616': 'r', '2424': 'c', '3232': 'm', '4848': 'y'}

fig1=plt.figure()
for vol in (svol,lvol):
  t_start = time.time()
  print "Working on %s:" % vol
  L = int(vol[:len(vol)/2])
  t = ((c*L)**2)/8
  tau = t + t0
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
    data[beta].append(f(tau))

  betal[vol]=(set(tmpbetal))
  for b in betal[vol]:
    E = np.mean(data[b])
    E_err = scipy.stats.sem(data[b])*math.sqrt(10) #temporary hack
    g2[(vol,b)] = (128*math.pi**2*E*t**2)/(3*(3**2-1)*(1+dc))
    g2_err[(vol,b)] = (128*math.pi**2*E_err*t**2)/(3*(3**2-1)*(1+dc))

    x = []
    y = []
    e = []
  for b in sorted(betal[vol]):
    x.append(float(b))
    y.append(float(g2[(vol,b)]))
    e.append(float(g2_err[(vol,b)]))
  ax = np.array(x)
  ay = np.array(y)
  ae = np.array(e)
  splt1 = fig1.add_subplot(111)
  splt1.errorbar(ax, ay, yerr=ae, linestyle='None', marker='.', color=pcolor[vol], label=vol)
  coeff[vol], success = optimize.leastsq(errfunc, c_in[:], args=(ax, ay, ae))
  splt1.plot(ax, fitfunc(coeff[vol],ax))
  splt1.legend(loc=3)
  t_finish = time.time()
  elapsed  = t_finish - t_start
  print "Took %s to parse\n.\t.\t.\t." % str(datetime.timedelta(seconds=elapsed))
fig1.savefig("plots/"+str(s)+"_"+str(c)+"_"+str(t0)+"_"+svol+"-"+lvol+"_wflow.png", format='png')

# - - - - - - - - - - - - - - - - - -

roots = np.roots(coeff[lvol]-coeff[svol])
possible = []
for i in roots:
  if ((i < float(sorted(betal[svol])[0])) or (i > float(sorted(betal[svol])[-1]))):
    next
  else:
    possible.append(i) 

if len(possible) == 1:
  beta_tmp = possible.pop()
  x = float(1)/float((int(svol[:len(svol)/2]))**2)
  crossing = np.polyval(coeff[lvol],beta_tmp).real
else:
  print "Did not find one crossing"
  sys.exit()

# - - - - - - - - - - - - - - - - - -
# find 4 corners of error region
err_crossing = []

# small + large +
for vol in (svol,lvol):
  x = []
  y = []
  e = []
  for b in sorted(betal[vol]):
    x.append(float(b))
    y.append(float(g2[(vol,b)]))
    e.append(float(g2_err[(vol,b)]))
  ax = np.array(x)
  ay = np.array(y)
  ae = np.array(e)
  ay = ay + ae
  coeff[vol], success = optimize.leastsq(errfunc, c_in[:], args=(ax, ay, ae))
roots = np.roots(coeff[lvol]-coeff[svol])
possible = []
for i in roots:
  if ((i < float(sorted(betal[svol])[0])) or (i > float(sorted(betal[svol])[-1]))):
    next
  else:
    possible.append(i) 
if len(possible) == 1:
  beta_tmp = possible.pop()
  x = float(1)/float((int(svol[:len(svol)/2]))**2)
  err_crossing.append(np.polyval(coeff[lvol],beta_tmp).real)
else:
  print "Did not find one crossing between the small volume plus errors and the large volume plus errors"


# small - large -
for vol in (svol,lvol):
  x = []
  y = []
  e = []
  for b in sorted(betal[vol]):
    x.append(float(b))
    y.append(float(g2[(vol,b)]))
    e.append(float(g2_err[(vol,b)]))
  ax = np.array(x)
  ay = np.array(y)
  ae = np.array(e)
  ay = ay - ae
  coeff[vol], success = optimize.leastsq(errfunc, c_in[:], args=(ax, ay, ae))
roots = np.roots(coeff[lvol]-coeff[svol])
possible = []
for i in roots:
  if ((i < float(sorted(betal[svol])[0])) or (i > float(sorted(betal[svol])[-1]))):
    next
  else:
    possible.append(i) 
if len(possible) == 1:
  beta_tmp = possible.pop()
  x = float(1)/float((int(svol[:len(svol)/2]))**2)
  err_crossing.append(np.polyval(coeff[lvol],beta_tmp).real)
else:
  print "Did not find one crossing between the small volume minus errors and the large volume minus errors"

# small - large +
for vol in (svol,lvol):
  x = []
  y = []
  e = []
  for b in sorted(betal[vol]):
    x.append(float(b))
    y.append(float(g2[(vol,b)]))
    e.append(float(g2_err[(vol,b)]))
  ax = np.array(x)
  ay = np.array(y)
  ae = np.array(e)
  if vol == svol:
    ay = ay - ae
  elif vol == lvol:
    ay = ay + ae
  coeff[vol], success = optimize.leastsq(errfunc, c_in[:], args=(ax, ay, ae))
roots = np.roots(coeff[lvol]-coeff[svol])
possible = []
for i in roots:
  if ((i < float(sorted(betal[svol])[0])) or (i > float(sorted(betal[svol])[-1]))):
    next
  else:
    possible.append(i) 
if len(possible) == 1:
  beta_tmp = possible.pop()
  x = float(1)/float((int(svol[:len(svol)/2]))**2)
  err_crossing.append(np.polyval(coeff[lvol],beta_tmp).real)
else:
  print "Did not find one crossing between the small volume minus errors and the large volume minus errors"

# small + large -
for vol in (svol,lvol):
  x = []
  y = []
  e = []
  for b in sorted(betal[vol]):
    x.append(float(b))
    y.append(float(g2[(vol,b)]))
    e.append(float(g2_err[(vol,b)]))
  ax = np.array(x)
  ay = np.array(y)
  ae = np.array(e)
  if vol == svol:
    ay = ay + ae
  elif vol == lvol:
    ay = ay - ae
  coeff[vol], success = optimize.leastsq(errfunc, c_in[:], args=(ax, ay, ae))
  roots = np.roots(coeff[lvol]-coeff[svol])
possible = []
for i in roots:
  if ((i < float(sorted(betal[svol])[0])) or (i > float(sorted(betal[svol])[-1]))):
    next
  else:
    possible.append(i) 
if len(possible) == 1:
  beta_tmp = possible.pop()
  x = float(1)/float((int(svol[:len(svol)/2]))**2)
  err_crossing.append(np.polyval(coeff[lvol],beta_tmp).real)
else:
  print "Did not find one crossing between the small volume minus errors and the large volume minus errors"

err_plus  = max(err_crossing)
err_minus = min(err_crossing)

print x, beta_tmp.real, crossing, err_minus, err_plus
