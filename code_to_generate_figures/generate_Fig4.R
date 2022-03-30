library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(cowplot)



source('abc_all_snps.R')
all_post = post
all_emp = emp
rm(list = setdiff(ls(), c("comp_models",'all_post','all_emp')))

source('abc_private_snps.R')
nosh_post = post
nosh_emp = emp
rm(list = setdiff(ls(), c("comp_models",'all_post','all_emp','nosh_post','nosh_emp')))

source('population_parameters.R')



########t sweep
all_post$snps = 'all'
nosh_post$snps = 'noshared'

d=density(all_post$t)
t_all = d$x[which.max(d$y)]
t_all

d=density(nosh_post$t)
t_nosh = d$x[which.max(d$y)]
t_nosh


######## credibility intervals
quantile(all_post$t, probs=c(0.025,0.975))
quantile(nosh_post$t, probs=c(0.025,0.975))


post = bind_rows(all_post,nosh_post)
post

tdens = ggplot(filter(post,t<0.5),aes(x = t,color=snps))+
  geom_density()+
  scale_color_manual(values = c('darkorange','darkorange4'))+
  #geom_vline(xintercept = t_all,color='darkorange')+
  #geom_vline(xintercept = t_nosh,color='darkorange4')+
  #scale_x_continuous(limits = c(0,0.2))+
  annotate('text', label= paste('All SNPs\nt =',round(t_all,4),"\n~", round(t_all* 4 *n0/10,0)," y.a."), x = 0.12,y = 50,color='darkorange',)+
  annotate('text', label= paste('Private SNPs\nt =',round(t_nosh,4),"\n~", round(t_nosh* 4 *n0/10,0)," y.a."), x = 0.04,y = 50,color='darkorange4',)+
  theme_minimal(base_size = 18)+
  labs(x=expression('t (4N'[e]*' generations)'),y='Density')+
  theme(legend.position = 'none')



tdens = tdens +
  scale_x_continuous(limits = c(0,0.2),sec.axis = sec_axis(trans = ~.* 4 *n0/10, name = 't (years ago)'))+
  geom_vline(xintercept = t_all ,linetype='dashed',color='darkorange')+
  geom_vline(xintercept = t_nosh,linetype='dashed',color='darkorange4')

tdens

