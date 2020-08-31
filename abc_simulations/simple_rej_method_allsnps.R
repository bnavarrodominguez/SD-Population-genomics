wd = "~/Documentos/investigacion/ur_larracuente_presgraves/SD_SNPs/ms_abc_simulations/ms_sample_rejection/"
setwd(wd)
#dir.create('model_selection_plots/')

######Load package
#library("abc")
library(readr)
library(dplyr)

######Pop parameters
source("population_parameters.R")

##### Load empirical (observed) summary statistics
emp <- read_delim("sims_10kacc_allsnps_in2rmal_overall/sumst_ZI_in2rmal_overall_allsnps.txt", 
                  " ", escape_double = FALSE, col_names = FALSE, 
                  trim_ws = TRUE)
colnames(emp) = c("pi.m","s.m","tajD.m")




###########Load data from simulations
acc <- read_delim("sims_10kacc_allsnps_in2rmal_overall/sims_tol005_allsnps_posterior.txt",
                    " ", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE,comment = "#")
colnames(acc) = c("t", "pi.m", "pi.sd", "s.m", "s.sd" ,"tajD.m", "TajD.sd")

rej <- read_delim("sims_10kacc_allsnps_in2rmal_overall/sims_tol005_allsnps_rejected.txt",
                  " ", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE,comment = "#")
colnames(rej) = c("t", "pi.m", "pi.sd", "s.m", "s.sd" ,"tajD.m", "TajD.sd")


post.pr = nrow(acc)/sum(nrow(acc)+nrow(rej))
post.pr
post.pr*100

post = acc


#######Check the posterior
par(mfrow=c(2,2))

hist(post$pi.m)
abline(v = emp$pi.m,col='red')

hist(post$s.m)
abline(v = emp$s.m,col='red')

hist(post$tajD.m)
abline(v = emp$tajD.m,col='red')


####t-sweep
#par(mfrow=c(1,1))
d =density(post$t)
plot(d)

tmode = d$x[which.max(d$y)]
tmode
abline(v = tmode,col='blue')



tsweep = tmode * 4 *n0/10
tsweep

min(post$t)* 4 *n0/10
max(post$t)* 4 *n0/10



nrow(acc)+nrow(rej)






##############Goodness of fit
gof <- read_delim("sims_10kacc_allsnps_in2rmal_overall/gof_sim_2905_0.096480_allsnps.txt", 
                  " ", escape_double = FALSE, col_names = FALSE, 
                  comment = "#", trim_ws = TRUE)
colnames(gof) = c('pi','s','tajD')
par(mfrow=c(1,2))
hist(gof$pi)
abline(v=emp$pi.m,col='red')

hist(gof$tajD)
abline(v=emp$tajD.m,col='red')

p_pi = ecdf(gof$pi)
p_pi(emp$pi.m)

plot(ecdf(gof$pi))
abline(v = emp$pi.m,col='red')

p_tajd = ecdf(gof$tajD)
p_tajd(emp$tajD.m)

plot(ecdf(gof$tajD))
abline(v = emp$tajD.m,col='red')


##########
