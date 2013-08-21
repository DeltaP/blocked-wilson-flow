import argparse
import glob
import re

parser = argparse.ArgumentParser(description='Strip Wflow Data')
parser.add_argument('v', metavar='volume', type=str, nargs='+', help='Enter the volumes that you would like to strip.')

args     = parser.parse_args()
argdict  = vars(args)
vol_list = argdict['v']

for vol in vol_list:
  print "Stripping data for volume:", vol
  filelist = glob.glob('12flav_'+vol+'/wflow/Wflow_'+'*'+vol+'*')
  for filename in filelist:
    split1   = re.split('/', filename)
    split2   = re.split('_', split1[2])
    outfile  = split1[0]+'/'+split1[1]+'/dat/wflow_'+split2[2]+'_'+split2[3]+'_'+split2[5]
    f = open(outfile, 'w')
    f.write('t td E t^2E der(t^2E) check topo\n')
    for l in open(filename):
      if l.startswith('WFLOW'):
        f.write(l.lstrip('WFLOW '))
