#!/usr/bin/python

####Dependencies
## ms http://home.uchicago.edu/rhudson1/source/mksamples.html
##
##

from __future__ import print_function
import random
import os
import math
import sys
import numpy
import csv
import numpy as np
from scipy import stats

print("\nUsage run_ms.model_fitting.py Nsims OutPrefix ObsS tsweep alpha")

##### Variable input
try:
    nsims = sys.argv[1]
except:
    nsims = raw_input("Number of simulations ")
try:
    out = sys.argv[2]
except:
    out = raw_input("Out prefix: ")

try:
    obs_s = sys.argv[3]
except:
    obs_s = raw_input("Observed S: ")

try:
    tsweep = sys.argv[4]
except:
    tsweep = raw_input("Estimated tsweep: ")

try:
    alpha = sys.argv[5]
except:
    alpha = raw_input("Exponential growth (alpha): ")


nsims = int(nsims)
obs_s = int(obs_s)
tsweep = float(tsweep)
alpha = float(alpha)

####Define fixed parameters

### D. melanogaster Ne in Zambia (Kapopoulou 2018)
ne = 3160475
print("D. melanogaster Ne in Zambia (Kapopoulou 2018): ", ne)

#####Mutation rate per site per generation (Laurent et al 2011, used in Kapopoulou 2018)
mutrate = 1.3e-9
print("Mutation rate per site per generation (Laurent et al 2011): ", mutrate)

### Frequency of In(2R)Mal - 1.47%
freq = 0.0147
print("Frequency of In(2R)Mal: ", freq)

##### Sequence lenght
#seql = 6000000
#print("Lenght of In(2R)Mal: ", seql)


##### Theta
#theta = 4*ne*mutrate*freq*seql
#print("Theta: ", theta)

#### N0 = Ne * frequency of In(2R)Mal
n0 = ne*freq
print ("N0: ", n0)

## proportion in reduction of population size in the sweep
### since n1 = 1, x1 must be 1/n0
x1 = 1/n0
n1 = x1*n0
print ("N1: ", n1)


## Parameters for goodness of fit test
### Observed S
print ("Observed S: ", obs_s)

### Observed S
print ("Estimated t sweep (4Ne gen): ", tsweep)

#### Observed
#obs_txt = open('%s' % sumstfile , 'r')
#obs = obs_txt.read().split(' ')
#obs_pi = float(obs[0])
#obs_s = float(obs[1])
#obs_tajd = float(obs[2])
#
###### Stablish thresholds for acceptance / prior on S
#
#thr = 0.05
#max_pi = obs_pi + thr*obs_pi
#min_pi = obs_pi - thr*obs_pi
#max_s = int(obs_s + thr*obs_s)
#min_s = int(obs_s - thr*obs_s)
#max_tajd = obs_tajd + abs(thr*obs_tajd)
#min_tajd = obs_tajd - abs(thr*obs_tajd)
#
#
##print obs
##print obs_pi
##print obs_s
##print obs_tajd

#print("Empirical summary statistics are Pi = %f, S = %d and Tajima's D = %f" % (obs_pi, obs_s, obs_tajd))



#########Simulation parameters
#nreps = 1000
nreps = nsims
nsam = 9



#### sweep at time t
os.system("./ms %d %d -s %d -eN %f %f | ./sample_stats |  awk '{print $2,$4,$6}' >> gof_sim_%d_%f_%s_sweep.txt" % (nsam, nreps, obs_s, tsweep, x1, obs_s,tsweep,out))

#### constant population size

os.system("./ms %d %d -s %d | ./sample_stats |  awk '{print $2,$4,$6}' >> gof_sim_%d_%s_const.txt" % (nsam, nreps, obs_s, obs_s, out))

#### exponential growth
os.system("./ms  %d %d -s %d -G %f | ./sample_stats |  awk '{print $2,$4,$6}' >> gof_sim_%d_%f_%s_expgrowth.txt" % (nsam, nreps, obs_s, alpha , obs_s, alpha, out))


 



