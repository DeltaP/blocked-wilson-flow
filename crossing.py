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
import time
import datetime
from collections import defaultdict
from scipy.optimize import curve_fit
from scipy import optimize 

def fitFunc(x,c0,c1,c2,c3):
  return 1/x - c0 - c1*x - c2*x**2 - c3*(x)**3

def fitFunc2(x,c0,c1,c2,c3):
  return 1 - c0*x - c1*x**2 - c2*x**3 - c3*(x)**4

def invfitFunc(x,c0,c1,c2,c3):
  return 1/(1/x - c0 - c1*x - c2*x**2 - c3*x**3)

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
ayerr    = defaultdict(dict)
yfit     = defaultdict(dict)
xfit     = defaultdict(dict)
g2       = defaultdict(dict)
g2_err   = defaultdict(dict)
dg2      = defaultdict(dict)
coeff    = defaultdict(dict)
pcolor   = {'66': 'y', '88': 'b', '1212': 'g', '1616': 'r','1818': 'y', '2424': 'c', '3232': 'm', '4848': 'y'}

fig1=plt.figure(1)
fig2=plt.figure(2)
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
    if float(beta) < 4.0:    #temporary for test
      continue               #temporary for test
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
  index = np.arange(0,len(x),1)
  ax = np.array(x)[::-1]
  ay = np.array(y)[::-1]
  ae = np.array(e)[::-1]
  ix = 12.0/ax
  iy = 1.0/ay
  ie = ae*iy/ay
  splt1 = fig1.add_subplot(111)
  splt1.errorbar(ix, iy, yerr=ie, linestyle='None', marker='.', color=pcolor[vol], label=vol)
  splt1.legend(loc=3)
  splt1.set_title('Inverted to perform fits')
  splt1.set_xlabel('12/beta')
  splt1.set_ylabel('1/g^2')
  splt2 = fig2.add_subplot(111)
  splt2.errorbar(ax, ay, yerr=ae, linestyle='None', marker='.', color=pcolor[vol], label=vol)
  splt2.legend(loc=3)
  splt2.set_title('Data with inverted fits')
  splt2.set_xlabel('beta')
  splt2.set_ylabel('g^2')
  fitParams, fitCovariances = curve_fit(fitFunc, ix, iy, sigma=ie)
  coeff[vol] = fitParams
  coeff[vol] = np.insert(coeff[vol],0,1)
  
  yerr = []
  for i in index:
    xvec = np.array([1,ix[i],ix[i]**2,ix[i]**3])
    yerr.append(math.sqrt(np.dot(np.dot(xvec,fitCovariances),xvec)))
  ayerr[vol] = np.array(yerr)

  yerr2p = []
  yerr2m = []
  inv_yfit = invfitFunc(ix, fitParams[0], fitParams[1], fitParams[2], fitParams[3])
  inv_yerr = ay * ayerr[vol] / iy

  xfit[vol] = ix
  yfit[vol] = fitFunc(ix, fitParams[0], fitParams[1], fitParams[2], fitParams[3])
  splt1.plot(ix, yfit[vol],color=pcolor[vol])
  splt1.fill_between(ix, yfit[vol]-ayerr[vol],yfit[vol]+ayerr[vol],color=pcolor[vol])
  splt2.plot(ax, inv_yfit,color=pcolor[vol])
  splt2.fill_between(ax, inv_yfit-inv_yerr,inv_yfit+inv_yerr,color=pcolor[vol])

  t_finish = time.time()
  elapsed  = t_finish - t_start
  print "Took %s to parse\n.\t.\t.\t." % str(datetime.timedelta(seconds=elapsed))
fig1.savefig("plots/crossing/wflow/"+str(s)+"_"+str(c)+"_"+str(t0)+"_"+svol+"-"+lvol+"_invrt.png", format='png')
fig2.savefig("plots/crossing/wflow/"+str(s)+"_"+str(c)+"_"+str(t0)+"_"+svol+"-"+lvol+"_wflow.png", format='png')

# - - - - - - - - - - - - - - - - - -

roots = np.roots(coeff[lvol][::-1]-coeff[svol][::-1])
possible = []
for i in roots:
  im = i.imag
  if ((im != 0) or (i < xfit[svol][0]) or (i < xfit[lvol][0]) or (i > xfit[svol][-1]) or (i > xfit[lvol][-1])):
    next
  else:
    possible.append(i.real) 

if len(possible) == 1:
  crossing_x = possible.pop()
  il = float(1)/float((int(svol[:len(svol)/2]))**2)
  crossing_y = 1/(fitFunc2(crossing_x,coeff[svol][1],coeff[svol][2],coeff[svol][3],coeff[svol][4])/crossing_x)
  crossing_x = 12.0/crossing_x
else:
  print "Did not find one crossing:"
  print possible
  #sys.exit()

# - - - - - - - - - - - - - - - - - -
# find 4 corners of error region
err_crossing = []

# small + large +
ys=yfit[svol]+ayerr[svol]
fitParams, fitCovariances = curve_fit(fitFunc, xfit[svol], yfit[svol]+ayerr[svol])
scoeff = fitParams
scoeff = np.insert(scoeff,0,1)

yl=yfit[lvol]+ayerr[lvol]
fitParams, fitCovariances = curve_fit(fitFunc, xfit[lvol], yfit[lvol]+ayerr[lvol])
lcoeff = fitParams
lcoeff = np.insert(lcoeff,0,1)

roots = np.roots(scoeff[::-1]-lcoeff[::-1])
possible = []
for i in roots:
  im = i.imag
  if ((im != 0) or (i < xfit[svol][0]) or (i < xfit[lvol][0]) or (i > xfit[svol][-1]) or (i > xfit[lvol][-1])):
    next
  else:
    possible.append(i.real) 
if len(possible) == 1:
  beta_tmp = possible.pop()
  il = float(1)/float((int(svol[:len(svol)/2]))**2)
  err_crossing.append(1/(fitFunc2(beta_tmp,coeff[svol][1],coeff[svol][2],coeff[svol][3],coeff[svol][4])/beta_tmp))
else:
  print "Did not find one crossing between the small volume plus errors and the large volume plus errors:"
  print possible
  for p in possible:
    print 1/(fitFunc2(p,coeff[svol][1],coeff[svol][2],coeff[svol][3],coeff[svol][4])/p)


# small - large -
fitParams, fitCovariances = curve_fit(fitFunc, xfit[svol], yfit[svol]-ayerr[svol])
scoeff = fitParams
scoeff = np.insert(scoeff,0,1)

fitParams, fitCovariances = curve_fit(fitFunc, xfit[lvol], yfit[lvol]-ayerr[lvol])
lcoeff = fitParams
lcoeff = np.insert(lcoeff,0,1)

roots = np.roots(scoeff[::-1]-lcoeff[::-1])
possible = []
for i in roots:
  im = i.imag
  if ((im != 0) or (i < xfit[svol][0]) or (i < xfit[lvol][0]) or (i > xfit[svol][-1]) or (i > xfit[lvol][-1])):
    next
  else:
    possible.append(i.real) 
if len(possible) == 1:
  beta_tmp = possible.pop()
  il = float(1)/float((int(svol[:len(svol)/2]))**2)
  err_crossing.append(1/(fitFunc2(beta_tmp,coeff[svol][1],coeff[svol][2],coeff[svol][3],coeff[svol][4])/beta_tmp))
else:
  print "Did not find one crossing between the small volume minus errors and the large volume minus errors:"
  print possible
  for p in possible:
    print 1/(fitFunc2(p,coeff[svol][1],coeff[svol][2],coeff[svol][3],coeff[svol][4])/p)

# small - large +
fitParams, fitCovariances = curve_fit(fitFunc, xfit[svol], yfit[svol]-ayerr[svol])
scoeff = fitParams
scoeff = np.insert(scoeff,0,1)

fitParams, fitCovariances = curve_fit(fitFunc, xfit[lvol], yfit[lvol]+ayerr[lvol])
lcoeff = fitParams
lcoeff = np.insert(lcoeff,0,1)

roots = np.roots(scoeff[::-1]-lcoeff[::-1])
possible = []
for i in roots:
  im = i.imag
  if ((im != 0) or (i < xfit[svol][0]) or (i < xfit[lvol][0]) or (i > xfit[svol][-1]) or (i > xfit[lvol][-1])):
    next
  else:
    possible.append(i.real) 
if len(possible) == 1:
  beta_tmp = possible.pop()
  il = float(1)/float((int(svol[:len(svol)/2]))**2)
  err_crossing.append(1/(fitFunc2(beta_tmp,coeff[svol][1],coeff[svol][2],coeff[svol][3],coeff[svol][4])/beta_tmp))
else:
  print "Did not find one crossing between the small volume minus errors and the large volume plus errors:"
  print possible
  for p in possible:
    print 1/(fitFunc2(p,coeff[svol][1],coeff[svol][2],coeff[svol][3],coeff[svol][4])/p)

# small + large -
fitParams, fitCovariances = curve_fit(fitFunc, xfit[svol], yfit[svol]+ayerr[svol])
scoeff = fitParams
scoeff = np.insert(scoeff,0,1)

fitParams, fitCovariances = curve_fit(fitFunc, xfit[lvol], yfit[lvol]-ayerr[lvol])
lcoeff = fitParams
lcoeff = np.insert(lcoeff,0,1)

roots = np.roots(scoeff[::-1]-lcoeff[::-1])
possible = []
for i in roots:
  im = i.imag
  if ((im != 0) or (i < xfit[svol][0]) or (i < xfit[lvol][0]) or (i > xfit[svol][-1]) or (i > xfit[lvol][-1])):
    next
  else:
    possible.append(i.real) 
if len(possible) == 1:
  beta_tmp = possible.pop()
  il = float(1)/float((int(svol[:len(svol)/2]))**2)
  err_crossing.append(1/(fitFunc2(beta_tmp,coeff[svol][1],coeff[svol][2],coeff[svol][3],coeff[svol][4])/beta_tmp))
else:
  print "Did not find one crossing between the small volume plus errors and the large volume minus errors:"
  print possible
  for p in possible:
    print 1/(fitFunc2(p,coeff[svol][1],coeff[svol][2],coeff[svol][3],coeff[svol][4])/p)

err_plus  = max(err_crossing)
err_minus = min(err_crossing)

print il, crossing_x.real, crossing_y, err_minus, err_plus
