######Load package
library(readr)
library(dplyr)

######Pop parameters
source("population_parameters.R")

##### Load empirical (observed) summary statistics
emp <- read_delim("simulations_data/sumst_ZI_in2rmal_noncoding_private.txt", 
                  " ", escape_double = FALSE, col_names = FALSE, 
                  trim_ws = TRUE)
colnames(emp) = c("pi.m","s.m","tajD.m")




###########Load data from simulations
acc <- read_delim("simulations_data/posterior_10ksamples_noncoding_private.txt",
                    " ", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE,comment = "#")
colnames(acc) = c("t", "pi.m", "pi.sd", "s.m", "s.sd" ,"tajD.m", "TajD.sd")

rej <- read_delim("simulations_data/rejected_noncoding_private.txt",
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





