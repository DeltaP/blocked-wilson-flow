#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import glob
import re
import sys
import getopt
import argparse
from collections import defaultdict

s = 2
dc = 0
c = 0.3

parser = argparse.ArgumentParser(description='Wilson Flow Matching.')
parser.add_argument('volume', metavar='volume', type=int, nargs=1, help='Enter the volume you are intersted in.')
args = parser.parse_args()
argdict = vars(args)
volume = argdict['volume'].pop()

t2E = 

g2 = (128*pi**2*t2E) / (3*(3-1)*(1+dc)
