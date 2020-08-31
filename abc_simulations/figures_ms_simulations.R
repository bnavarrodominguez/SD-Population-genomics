library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(cowplot)

source('compare_sweep_and_model.R')
comp_models = datos
rm(list=setdiff(ls(), "comp_models"))

source('simple_rej_method_allsnps.R')
all_post = post
all_emp = emp
rm(list = setdiff(ls(), c("comp_models",'all_post','all_emp')))

source('simple_rej_method_noshared.R')
nosh_post = post
nosh_emp = emp
rm(list = setdiff(ls(), c("comp_models",'all_post','all_emp','nosh_post','nosh_emp')))

source('population_parameters.R')


#####boxplot models

head(comp_models)

pi = comp_models %>% select(model, pi.m)

a = ggplot(pi,aes(x = model,y = pi.m))+
  geom_boxplot()+
  #scale_color_manual(values=c('darkgreen','darkorange'))+
  labs(x=NULL, y= "Pi")+
  scale_x_discrete(labels=c("Const","Sweep"))+
  geom_hline(yintercept = all_emp$pi.m, color='red')+
  theme_minimal(base_size = 16)+
  theme(legend.position = 'none')
a

tajd = comp_models %>% select(model, tajD.m)

b = ggplot(tajd,aes(x = model,y = tajD.m))+
  geom_boxplot()+
  #scale_color_manual(values=c('darkgreen','darkorange'))+
  labs(x=NULL, y= "Tajima's D")+
  scale_x_discrete(labels=c("Const","Sweep"))+
  geom_hline(yintercept = all_emp$tajD.m, color='red')+
  theme_minimal(base_size = 16)+
  theme(legend.position = 'none')

b



require(cowplot)
plot_grid(a,b,labels = 'auto')





########t sweep
all_post$snps = 'all'
nosh_post$snps = 'noshared'

d=density(all_post$t)
t_all = d$x[which.max(d$y)]
t_all

d=density(nosh_post$t)
t_nosh = d$x[which.max(d$y)]
t_nosh



post = bind_rows(all_post,nosh_post)

tdens = ggplot(post,aes(x = t,color=snps))+
  geom_density()+
  scale_color_manual(values = c('darkorange','darkorange4'))+
  #geom_vline(xintercept = t_all,color='darkorange')+
  #geom_vline(xintercept = t_nosh,color='darkorange4')+
  scale_x_continuous(limits = c(0,0.2))+
  annotate('text', label= paste('All SNPs','t=0.0965',sep = '\n'), x = 0.12,y = 50,color='darkorange')+
  annotate('text', label= paste('Private SNPs','t=0.0737',sep = '\n'), x = 0.05,y = 50,color='darkorange4')+
  theme_minimal(base_size = 18)+
  labs(x=expression(t[sweep]),y='Density')+
  theme(legend.position = 'none')

tdens 


ggplot(post,aes(x = pi.m, color=snps,fill=snps))+
  geom_histogram(alpha=.25)+
  scale_color_manual(values = c('darkorange','darkorange4'))+
  scale_fill_manual(values = c('darkorange','darkorange4'))+
  geom_vline(xintercept = all_emp$pi.m,color='darkorange')



