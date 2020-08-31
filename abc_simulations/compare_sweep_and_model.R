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
emp



###########Load data from simulations
acc <- read_delim("sims_10kacc_allsnps_in2rmal_overall/sims_tol005_allsnps_posterior.txt",
                    " ", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE,comment = "#")
colnames(acc) = c("t", "pi.m", "pi.sd", "s.m", "s.sd" ,"tajD.m", "TajD.sd")

rej <- read_delim("sims_10kacc_allsnps_in2rmal_overall/sims_tol005_allsnps_rejected.txt",
                  " ", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE,comment = "#")
colnames(rej) = c("t", "pi.m", "pi.sd", "s.m", "s.sd" ,"tajD.m", "TajD.sd",'post',"sobra")

rej$sobra=NULL
tot = nrow(acc)+nrow(rej)
post.pr = nrow(acc)/tot
post.pr
post.pr*100

acc$post = 'acc'

sweep = bind_rows(acc,rej)
sweep$model = 'sweep'

###########Load data from simulations (constant model)
acc2 <- read_delim("sims_allsnps_constantsize_in2rmal_overall/sims_const_size_allsnps_accepted.txt",
                  " ", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE,comment = "#")
colnames(acc2) = c("t", "pi.m", "pi.sd", "s.m", "s.sd" ,"tajD.m", "TajD.sd")

rej2 <- read_delim("sims_allsnps_constantsize_in2rmal_overall/sims_const_size_allsnps_rejected.txt",
                  " ", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE,comment = "#")
colnames(rej2) = c("t", "pi.m", "pi.sd", "s.m", "s.sd" ,"tajD.m", "TajD.sd",'post',"sobra")

rej2$sobra=NULL
post.pr2 = nrow(acc2)/sum(nrow(acc2)+nrow(rej2))
post.pr2*100

rej2 = head(rej2,n = nrow(sweep))
#acc$post = 'acc'


const = bind_rows(acc2,rej2)
const$model = 'const'


datos =bind_rows(sweep,const)

par(mfrow=c(1,2))
boxplot(datos$pi.m~datos$model)
abline(h=emp$pi.m,col='red')

boxplot(datos$tajD.m~datos$model)
abline(h=emp$tajD.m,col='red')

# datos$color = ifelse(datos$post=='acc','red','black')
# datos$shape = ifelse(datos$model=='sweep',1,2)
# plot(datos$pi.m,datos$tajD.m,col=datos$color,pch=datos$shape)




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
