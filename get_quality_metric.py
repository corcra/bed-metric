#!/usr/bin/python
# * ----------------------------------*
# |  Script to get 'quality metric'   |
# |  for data in bed format.          |
# |                                   |
# |  Stephanie Hyland, Nov 2013.      |
# |  Tri-I CBM sh985@cornell.edu      |
# * ----------------------------------*

import random
import gzip
import sys
import os
import argparse
from subprocess import call
from copy import copy

random.seed()

parser = argparse.ArgumentParser(description='Calculate quality metric!')
parser.add_argument('file',help='bed.gz file! (must be sorted - use sort-bed)',metavar='DATA.bed.gz')
parser.add_argument('--plot',help='plot #unique v #total with (predicted) asymptote!',action='store_true')

args = parser.parse_args()
data_path = args.file
plot_true = args.plot
if plot_true:
    plot = '1'
else:
    plot = '0'

data_name = os.path.basename(data_path)
record_path = data_path+'.record'
recordfile = open(record_path,'w')

master_dict = dict()
loc_buffer = {'+': 0, '-': 0}
n_uniq = {'+': 0, '-': 0}
n_total = {'+': 0, '-': 0}

frac_list = [float(frac)/40 for frac in range(1,40)]
for frac in frac_list:
    master_dict[frac] = (copy(loc_buffer),copy(n_uniq),copy(n_total))

i=0
for line in gzip.open(data_path):
    if i%1000000==0:
        print i
    i+=1

    if 'rRNA' in line:
        continue

    for frac in frac_list:
        r = random.random()
        if r<=frac:
            strand = line.split()[5]
            master_dict[frac][2][strand] +=1            # increasing the total count for this strand, for this frac
            loc = int(line.split()[1])
            if loc != master_dict[frac][0][strand]:     # is it the same as what was recorded previously? if no, register a unique hit!
                master_dict[frac][1][strand]+=1         # increasing the unique count for this strand, for this frac
                master_dict[frac][0][strand] = loc          # recording new location

print 'Summary:'
for frac in sorted(master_dict):
    n_uniq = master_dict[frac][1]
    n_total = master_dict[frac][2]
    n_uniq_both=0
    n_total_both=0
    for strand in n_uniq:
        print frac, strand, n_uniq[strand], n_total[strand] 
        n_uniq_both = n_uniq[strand] + n_uniq_both
        n_total_both = n_total[strand] + n_total_both
    recordfile.write(data_name+'\t'+str(frac)+'\t'+str(n_uniq_both)+'\t'+str(n_total_both)+'\n')

recordfile.close()

script_dir = os.path.dirname(os.path.realpath(sys.argv[0]))
R_script_path = script_dir+'/calc_quality_metric.r'
Rscript = 'R --slave --file='+R_script_path+' --args '+record_path+' '+plot+' '+data_path+' '+data_name
call(Rscript,shell=True)
