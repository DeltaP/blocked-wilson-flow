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

c = 0.3

parser = argparse.ArgumentParser(description='Wilson Flow Matching.')
parser.add_argument('volume', metavar='volume', type=int, nargs=1, help='Enter the volume you are intersted in.')
parser.add_argument('tag', metavar='tag', type=str, nargs=1, help='Enter the base name of the run you are intersted in.')
parser.add_argument('tsmear', metavar='tsmear', type=float, nargs=1, help='Enter the Wflow smearing time you are interested in.')
args = parser.parse_args()
argdict = vars(args)
volume = argdict['volume'].pop()
tag = argdict['tag'].pop()
tsmear = argdict['tsmear'].pop()
