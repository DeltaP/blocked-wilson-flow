#!/usr/bin/python

from pandas import Series, DataFrame
import pandas as pd
import glob
import re
import numpy

flavor = 4
mass = 0.0

vol = []
b = []
m_t = []
w_t = []
val = []

#tup_vol = (1224, 1632, 2448)
tup_vol = 1224, 1632
for v in tup_vol:
  filelist = glob.glob(str(flavor)+'flav_'+str(v)+'/BlockedWflow_low_'+str(v)+'_*_0.0.*')
  for filename in filelist:
    print 'Reading filename:  '+filename
    f = open(filename, 'r')
    junk, start, vv, beta, junk, mass, mont_t = re.split('_', filename)
    ftext = f.readlines()
    for line in ftext:
      if re.match('^WFLOW.*', line):
        line=line.strip()
        junk, smear_t, junk, junk, wilson_flow, junk, junk, junk = re.split('\s+', line)
        vol.append(v)
        b.append(beta)
        m_t.append(mont_t)
        w_t.append(smear_t)
        val.append(wilson_flow)
zipped = zip(vol, beta, m_t, w_t)
index = MultiIndex.from_tuples(zipped, names=['volume', 'beta', 'montecarlo_time, smearing_time'])
data = Series(val, index=index)
print data
        
        


