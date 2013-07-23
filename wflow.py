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
fitfunc = lambda c, x: c[3] + c[2]*x + c[1]*x**2 + c[0]*x**3
errfunc = lambda c, x, y, err: (fitfunc(c, x) - y) / err
c_in = [1.0, 1.0, 1.0, 1.0]   # Order-of-magnitude initial guesses

parser = argparse.ArgumentParser(description='Wilson Flow Matching.')
parser.add_argument('c', metavar='constant', type=float, nargs=1, help='Enter the smearing constant that defines the extent of the smearing you are interested in.')
parser.add_argument('volumes', metavar='volume', type=str, nargs= '+', help='Enter the volume(s) you are intersted in.')
args = parser.parse_args()
argdict = vars(args)
c = argdict['c'].pop()
vol = defaultdict(list)
vol = argdict['volumes']

betal    = defaultdict(list)
g2       = defaultdict(dict)
g2_err   = defaultdict(dict)
dg2      = defaultdict(dict)
coeff    = defaultdict(dict)
pcolor   = {'88': 'b', '1212': 'g', '1616': 'r', '2424': 'c', '3232': 'm', '4848': 'y'}

fig1=plt.figure()
for volume in vol:
  L = int(volume[:len(volume)/2])
  t = round(((c*L)**2/8)/(0.02))*0.02
  data     = defaultdict(list)
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
  print volume
  print betal[volume]
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
  splt1 = fig1.add_subplot(111)
  splt1.errorbar(ax, ay, yerr=ae, linestyle='None', marker='.', color=pcolor[volume], label=volume)
  if ((volume != '4848')and(volume != '3232')):
    coeff[volume], success = optimize.leastsq(errfunc, c_in[:], args=(ax, ay, ae))
    splt1.plot(ax, fitfunc(coeff[volume],ax))
splt1.legend()
fig1.savefig("plots/wflow.png", format='png')

g0_3232_s2 = []
if (('88' in vol) and ('1616' in vol) and ('3232' in vol)):
  print "Matching 88, 1616, 3232"
  for b in sorted(betal['3232']):
    tempg = g2[('1616',b)]
    g0rnd = round(tempg,3)
    g0_3232_s2.append(g0rnd)
    print b + "\t" + str(g0rnd)
    dg2[('2','1616',g0rnd)] = (g2[('3232',b)]-g2[('1616',b)])/(2*np.log(2))
    x = []
    y = []
    e = []
    for bb in sorted(betal['88']):
      x.append(float(bb))
      yval = float(g2[('88',bb)])-float(tempg)
      y.append(yval)
      e.append(float(g2_err[('88',bb)]))
    nx = np.array(x)
    ny = np.array(y)
    ne = np.array(e)
    scoeff, success = optimize.leastsq(errfunc, c_in[:], args=(nx, ny, ne))
    roots = np.roots(scoeff)
    possible = []
    for i in roots:
      if ((i < float(sorted(betal['88'])[0])) or (i > float(sorted(betal['88'])[-1]))):
        next
      else:
        possible.append(i) 
    if len(possible) == 1:
      beta_tmp = possible.pop()
      dg2[('2','88',g0rnd)] = (np.polyval(coeff['1616'],beta_tmp)-np.polyval(coeff['88'],beta_tmp))/(2*np.log(2))
    else:
      print "Did not find one crossing"
    x = []
    y = []
    e = []
    for bb in sorted(betal['1212']):
      x.append(float(bb))
      yval = float(g2[('1212',bb)])-float(tempg)
      y.append(yval)
      e.append(float(g2_err[('1212',b)]))
    nx = np.array(x)
    ny = np.array(y)
    ne = np.array(e)
    scoeff, success = optimize.leastsq(errfunc, c_in[:], args=(nx, ny, ne))
    roots = np.roots(scoeff)
    possible = []
    for i in roots:
      if ((i < float(sorted(betal['1212'])[0])) or (i > float(sorted(betal['1212'])[-1]))):
        next
      else:
        possible.append(i) 
    if len(possible) == 1:
      beta_tmp = possible.pop()
      dg2[('2','1212',g0rnd)] = (np.polyval(coeff['2424'],beta_tmp)-np.polyval(coeff['1212'],beta_tmp))/(2*np.log(2))
    else:
      print "Did not find one crossing"

g0_4848_s2 = []
if (('1212' in vol) and ('2424' in vol) and ('4848' in vol)):
  print "Matching 1212, 2424, 4848"
  for b in sorted(betal['4848']):
    tempg = g2[('2424',b)]
    g0rnd = round(tempg,3)
    g0_4848_s2.append(g0rnd)
    print b + "\t" + str(g0rnd)
    dg2[('2','2424',g0rnd)] = (g2[('4848',b)]-g2[('2424',b)])/(2*np.log(2))
    x = []
    y = []
    e = []
    for bb in sorted(betal['1212']):
      x.append(float(bb))
      yval = float(g2[('1212',bb)])-float(tempg)
      y.append(yval)
      e.append(float(g2_err[('1212',bb)]))
    nx = np.array(x)
    ny = np.array(y)
    ne = np.array(e)
    scoeff, success = optimize.leastsq(errfunc, c_in[:], args=(nx, ny, ne))
    roots = np.roots(scoeff)
    possible = []
    for i in roots:
      if ((i < float(sorted(betal['1212'])[0])) or (i > float(sorted(betal['1212'])[-1]))):
        next
      else:
        possible.append(i) 
    if len(possible) == 1:
      beta_tmp = possible.pop()
      dg2[('2','1212',g0rnd)] = (np.polyval(coeff['2424'],beta_tmp)-np.polyval(coeff['1212'],beta_tmp))/(2*np.log(2))
    else:
      print "Did not find one crossing"
    x = []
    y = []
    e = []
    for bb in sorted(betal['88']):
      x.append(float(bb))
      yval = float(g2[('88',bb)])-float(tempg)
      y.append(yval)
      e.append(float(g2_err[('88',bb)]))
    nx = np.array(x)
    ny = np.array(y)
    ne = np.array(e)
    scoeff, success = optimize.leastsq(errfunc, c_in[:], args=(nx, ny, ne))
    roots = np.roots(scoeff)
    possible = []
    for i in roots:
      if ((i < float(sorted(betal['88'])[0])) or (i > float(sorted(betal['88'])[-1]))):
        next
      else:
        possible.append(i) 
    if len(possible) == 1:
      beta_tmp = possible.pop()
      dg2[('2','88',g0rnd)] = (np.polyval(coeff['1616'],beta_tmp)-np.polyval(coeff['88'],beta_tmp))/(2*np.log(2))
    else:
      print "Did not find one crossing"


if (('88' in vol) and ('1212' in vol) and ('1616' in vol) and ('2424' in vol) and ('3232' in vol) and ('4848' in vol)):
  print "Matching 88, 1212, 1616, 2424, 3232, 4848"


print "plot 3232 s=2"
for g0 in g0_3232_s2:
  x = []
  y = []
  x.append(1.0/16**2)
  y.append(dg2[('2', '1616',g0)])
  x.append(1.0/12**2)
  y.append(dg2[('2', '1212',g0)])
  x.append(1.0/8**2)
  y.append(dg2[('2', '88',g0)])

  ax = np.array(x)
  ay = np.array(y)
  fig2 = plt.figure()
  splt2 = fig2.add_subplot(111)
  splt2.plot(ax, ay)
  fig2.savefig('plots/3232_'+ str(g0)+'.png', format='png')

print "plot four.png s=2"
for g0 in g0_4848_s2:
  x = []
  y = []
  x.append(1.0/24**2)
  y.append(dg2[('2', '2424',g0)])
  #x.append(1.0/16**2)
  #y.append(dg2[('2', '1616',g0)])
  x.append(1.0/12**2)
  y.append(dg2[('2', '1212',g0)])
  x.append(1.0/8**2)
  y.append(dg2[('2', '88',g0)])

  ax = np.array(x)
  ay = np.array(y)
  fig3 = plt.figure()
  splt3 = fig3.add_subplot(111)
  splt3.plot(ax, ay)
  fig3.savefig('plots/4848_'+ str(g0)+'.png', format='png')
