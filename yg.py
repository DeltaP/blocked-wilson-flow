#!/usr/bin/python

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import math
import scipy.stats
import scipy
from collections import defaultdict
from scipy import optimize

fitfunc = lambda c, x: c[2] + c[1]*x + c[0]*x**2
errfunc = lambda c, x, y, err: (fitfunc(c, x) - y) / err
c_in = [1.0, 1.0, 1.0]   # Order-of-magnitude initial guesses

print "Scheme 7"
sch7 = pd.read_table('y_g/3lm_7_33',sep='\s+')

x7 = sch7['beta']
y7 = sch7['db']
em7 = sch7['db_md']
ep7 = sch7['db_pd']
e7 = ((em7-y7).abs()+(ep7-y7).abs())/2.0
print e7

plt.errorbar(x7,y7,e7, linestyle='None', marker='.', label='Sch 7', color='r')
coeff7, success = optimize.leastsq(errfunc, c_in[:], args=(x7, y7, e7))
plt.plot(x7, fitfunc(coeff7,x7), color='r')

roots = np.roots(coeff7)
possible = []
for i in roots:
  if ((i > 4.0) and (i < 7.5)):
    possible.append(i)
if len(possible) == 1:
  der = np.polyder(coeff7)
  slope = np.polyval(der,i)
  slope_err = np.polyval(der,i)*0.02
  print "Slope at %f is %f" % (i,slope)
  yg7 = -math.log(1-slope)/math.log(2)
  yg7e = -math.log(1-slope_err)/math.log(2)
  print "y_g is %f" % yg7 
  print "y_g is %f +/- %f" % (yg7, yg7e) 

print "Scheme 9"
sch9 = pd.read_table('y_g/3lm_9_33',sep='\s+')

x9 = sch9['beta']
y9 = sch9['db']
em9 = sch9['db_md']
ep9 = sch9['db_pd']
e9 = ((em9-y9).abs()+(ep9-y9).abs())/2.0
print e9

plt.errorbar(x9,y9,e9, linestyle='None', marker='.', label='Sch 9', color='g')
coeff9, success = optimize.leastsq(errfunc, c_in[:], args=(x9, y9, e9))
plt.plot(x9, fitfunc(coeff9,x9), color='g')

roots = np.roots(coeff9)
possible = []
for i in roots:
  if ((i > 4.0) and (i < 7.5)):
    possible.append(i)
if len(possible) == 1:
  der = np.polyder(coeff9)
  slope = np.polyval(der,i)
  print "Slope at %1f is %f" % (i,slope)
  yg9 = -math.log(1-slope)/math.log(2)
  print "y_g is %f +/- %f" % (yg9, yg9*0.02) 

print "Scheme 11"
sch11 = pd.read_table('y_g/3lm_11_33',sep='\s+')

x11 = sch11['beta']
y11 = sch11['db']
em11 = sch11['db_md']
ep11 = sch11['db_pd']
e11 = ((em11-y11).abs()+(ep11-y11).abs())/2.0
print e11

plt.errorbar(x11,y11,e11, linestyle='None', marker='.', label='Sch 11', color='b')
coeff11, success = optimize.leastsq(errfunc, c_in[:], args=(x11, y11, e11))
plt.plot(x11, fitfunc(coeff11,x11), color='b')
plt.legend()
plt.savefig('y_g/fit.png')

roots = np.roots(coeff11)
possible = []
for i in roots:
  if ((i > 4.0) and (i < 7.5)):
    possible.append(i)
if len(possible) == 1:
  der = np.polyder(coeff11)
  slope = np.polyval(der,i)
  print "Slope at %f is %f" % (i,slope)
  yg11 = -math.log(1-slope)/math.log(2)
  print "y_g is %f +/- %f" % (yg11, yg11*0.02) 

