#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import glob
import re
import sys
import getopt
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser(description='Produce plots of blocked observable(s) as a function of beta.')
parser.add_argument('tsmear', metavar='tsmear', type=float, nargs=1, help='Enter the Wflow smearing time you are interested in.')
parser.add_argument('rest', metavar='rest', type=str, nargs='+', help='Enter pairs of volumes and blocking steps.')
args = parser.parse_args()
argdict = vars(args)
tsmear = argdict['tsmear'].pop()
count = 0
blk = defaultdict(list)
vol = defaultdict(list)
while argdict['rest']:
  blk[count] = argdict['rest'].pop()
  vol[count] = argdict['rest'].pop()
  count += 1
print tsmear
print vol
print blk
