#!/usr/bin/python
#
##Changes from previous version:
##    - Saves rejected simulations to a file
##    - Acceptance based in mean stat from simulations within 5% of the observed mean (1% gets acceptance rate terribly low)
##    - Output stats on acceptance and rejection
##    - IMPORTANT: run the model with S instead of theta (put a prior on S too, limits by 1%)
##    - Flushes every line - solve memory flaws 
##    - Fixed Tajima's D issue
##    - Takes a third argument with the name of the sumst data file
##    - Lenght of In(2R)Mal set up to 9.5Mb - to use stats calculated in a single window approaching the lenght of the inversion
##
##
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

print("\nUsage run_ms.py NTotalSims OutputLabel SumStFile")

##### Variable input
try:
    nsims = sys.argv[1]
except:
    nsims = raw_input("Number of total samples: ")
try:
    out = sys.argv[2]
except:
    out = raw_input("Out prefix: ")

try:
    sumstfile = sys.argv[3]
except:
    sumstfile = raw_input("Observed SumSt File: ")

nsims = int(nsims)

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
seql = 9500000
print("Lenght of In(2R)Mal: ", seql)


##### Theta
theta = 4*ne*mutrate*freq*seql
print("Theta: ", theta)

#### N0 = Ne * frequency of In(2R)Mal
n0 = ne*freq
print ("N0: ", n0)

## proportion in reduction of population size in the bottleneck
### since n1 = 1, x1 must be 1/n0
x1 = 1/n0
n1 = x1*n0
print ("N1: ", n1)

#### Observed
obs_txt = open('%s' % sumstfile , 'r')
obs = obs_txt.read().split(' ')
obs_pi = float(obs[0])
obs_s = float(obs[1])
obs_tajd = float(obs[2])

##### Stablish thresholds for acceptance / prior on S

thr = 0.05
max_pi = obs_pi + thr*obs_pi
min_pi = obs_pi - thr*obs_pi
max_s = int(obs_s + thr*obs_s)
min_s = int(obs_s - thr*obs_s)
max_tajd = obs_tajd + abs(thr*obs_tajd)
min_tajd = obs_tajd - abs(thr*obs_tajd)


##print obs
##print obs_pi
##print obs_s
##print obs_tajd

print("Empirical summary statistics are Pi = %f, S = %d and Tajima's D = %f" % (obs_pi, obs_s, obs_tajd))



#########Simulation parameters
nreps = 1000
nsam = 9


print("Generating %s simulations of %s samples and %s repeats" % (nsims, nsam, nreps))

#outF = open(outfile, 'w')
#outF.write("#t mean_pi se_pi mean_S se_S mean_TajD se_TajD \n")

accF = open("acc_sim_%s.txt" % out, "w")
accF.write("#t mean_pi se_pi mean_S se_S mean_TajD se_TajD \n")

rejF = open("rej_sim_%s.txt" % out, "w")
rejF.write("#t mean_pi se_pi mean_S se_S mean_TajD se_TajD rejected \n")

i=0
acc=0
rejpi=0
rejtd=0
rejboth=0

########Generate the simulations
while i < nsims:
    ####t prior following a uniform distribution from 0.0 to 4Ne
    #t = random.uniform(0.0,1.0)
    #t = t+.1
    t = 0.00
    s = random.randint(min_s, max_s)
    ####generate a random distribution of S
    
        
##    with open("random_distrib_S_%s.txt" % out, "w") as sF:
##    for n in nreps:
##        s = random.randint(min_s,max_s)
##        sF.write(str(s)+ '\n')
##    sF.close()   

    #### Inizialize variables
    sim_pi = []
    sim_s = []
    sim_tajd = []

    #### Run ms with theta
    #os.system("./ms %d %d -t %f -eN %f %f | ./sample_stats |  awk '{print $2,$4,$6}' > sumstats_%s.tab" % (nsam, nreps, theta, t, x1, out))

    ### Run ms with fixed S
    #os.system("./ms %d %d -s %f -eN %f %f | ./sample_stats |  awk '{print $2,$4,$6}' > sumstats_%s.tab" % (nsam, nreps, s, t, x1, out))


    ### Run ms with constant population size
    os.system("./ms %d %d -s %f | ./sample_stats |  awk '{print $2,$4,$6}' > sumstats_%s.tab" % (nsam, nreps, s, out))


    with open("sumstats_%s.tab" % out) as sumst:
        reader = csv.reader(sumst, delimiter=' ')
        for row in reader:
            sim_pi.append(float(row[0]))
            sim_s.append(float(row[1]))
            sim_tajd.append(float(row[2]))
    sumst.close()
    #os.system("rm sumstats_%f.tab" % t)

    mean_pi = np.mean(sim_pi)
    se_pi = stats.sem(sim_pi)
    #ci_pi = stats.t.interval(0.95,len(sim_pi)-1, loc=mean_pi,scale=se_pi)
    mean_s = np.mean(sim_s)
    se_s = stats.sem(sim_s)
    #ci_s = stats.t.interval(0.95,len(sim_s)-1, loc=mean_s,scale=se_s)
    mean_tajd = np.mean(sim_tajd)
    se_tajd = stats.sem(sim_tajd)
    #ci_tajd = stats.t.interval(0.95,len(sim_tajd)-1, loc=mean_tajd,scale=se_tajd)

    if(min_pi < mean_pi < max_pi and min_tajd < mean_tajd < max_tajd) :
        print("Accepted simulations for t=%f" % t)
        #i = i+1
        acc = acc+1
        accF.write("%f %f %f %f %f %f %f\n" % (t,mean_pi,se_pi,mean_s,se_s,mean_tajd,se_tajd))
	accF.flush()
    elif((min_pi > mean_pi or mean_pi > max_pi) and (min_tajd > mean_tajd or mean_tajd > max_tajd)) :
        print("Rejected simulations for t=%f, both Pi and TajD off-limits" % t)
        rejboth=rejboth+1
        rejF.write("%f %f %f %f %f %f %f 'rej_both' \n" % (t,mean_pi,se_pi,mean_s,se_s,mean_tajd,se_tajd))
	rejF.flush()
    elif(min_pi < mean_pi < max_pi) :
        print("Rejected simulations for t=%f, TajD off-limits" % t)
        rejpi = rejpi + 1
        rejF.write("%f %f %f %f %f %f %f 'rej_tajd' \n" % (t,mean_pi,se_pi,mean_s,se_s,mean_tajd,se_tajd))
	rejF.flush()
    else :
        print("Rejected simulations for t=%f, Pi off-limits" % t)
        rejtd = rejtd+1
        rejF.write("%f %f %f %f %f %f %f 'rej_Pi' \n" % (t,mean_pi,se_pi,mean_s,se_s,mean_tajd,se_tajd))
	rejF.flush()

    i = i+1	


accF.close()
rejF.close()

os.system("rm sumstats_%s.tab" % out)


####Write summary
sumF = open("sum_sim_%s.txt" % out, "w")
tot_sims = i +rejpi+rejtd+rejboth
acc_rate = i/tot_sims

sumF.write("Total simulations: %d \n \
Accepted simulations: %d \n \
Acceptance rate: %f \n \
Rejected, Pi off-limits: %d \n \
Rejected, TajD off-limits: %d \n \
Rejected, both Pi and TajD off-limits: %d " % (tot_sims,i,acc_rate,rejpi,rejtd,rejboth))

sumF.close()
