#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import glob
import re
import sys
import getopt
import argparse
from collections import defaultdict

def uniq(seq, idfun=None):  
    # order preserving 
    if idfun is None: 
        def idfun(x): return x 
    seen = {} 
    result = [] 
    for item in seq: 
        marker = idfun(item) 
        # in old Python versions: 
        # if seen.has_key(marker) 
        # but in new ones: 
        if marker in seen: continue 
        seen[marker] = 1 
        result.append(item) 
    return result

class AbsAction(argparse.Action):
  def __call__(self, parser, namespace, values, option_string=None):
    if len(values) % 3 == 0:
      # If valid, store the values.
      setattr(namespace, self.dest, values)
    else:
      # Otherwise, invoke a parser error with a message.
      parser.error('scheme, volume, and blocking steps must be supplied in groups of three')

parser = argparse.ArgumentParser(description='Produce plots of blocked observable(s) as a function of beta.')
parser.add_argument('tsmear', metavar='tsmear', type=float, nargs=1, help='Enter the Wflow smearing time you are interested in.')
parser.add_argument('rest', metavar='scheme volume blocking', nargs = '+', action = AbsAction,  help='Enter tripplets of schemes, volumes, and blocking steps.')
args = parser.parse_args()
argdict = vars(args)
tsmear = argdict['tsmear'].pop()
count = 0
sch = defaultdict(list)
blk = defaultdict(list)
vol = defaultdict(list)
while argdict['rest']:
  blk[count] = argdict['rest'].pop()
  vol[count] = argdict['rest'].pop()
  sch[count] = argdict['rest'].pop()
  count += 1

value = defaultdict(list)
for i in range(count):
  betal=[]
  grab = '12flav_' + vol[i] + '/WMCRG' + sch[i] + '*' + vol[i] + '*' + '_0.0.*'
  print grab
  filelist = glob.glob(grab)
  for filename in filelist:
    junk, junk, junk, junk, beta, junk, junk = re.split('_', filename)
    betal.append(beta)
    f = open(filename, 'r')
    ftext = f.readlines()
    for line in ftext:
      if re.match('^LOOPS.*', line):
        line=line.strip()
        junk, smear_t, loop, junk, bk, junk, val = re.split('\s+', line)
        if float(smear_t) == tsmear:
          if loop == '0':
            if bk == blk[i]:
              value[beta].append(float(val))
  x = []
  y = []
  for b in sorted(uniq(betal)):
    np.append(x, b)
    np.append(y, np.mean(value[b]))
    plt.plot(x,y)
plt.show()
