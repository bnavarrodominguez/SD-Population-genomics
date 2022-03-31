#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#############Load libraries
require(dplyr)
require(stringr)
require(ggplot2)
require(ggforce)
require(ggpubr)
require(viridis)
require(cowplot)




priv_sites <- read.delim(args[1], header=FALSE)
shar_sites <- read.delim(args[2], header=FALSE)


priv_sites %>% group_by(V1) %>% dplyr::summarise(tot= n())
shar_sites %>% group_by(V1) %>% dplyr::summarise(tot=n())

shar_sites$id = interaction(shar_sites$V1,shar_sites$V2,sep = '_',drop = TRUE)
priv_sites$id = interaction(priv_sites$V1,priv_sites$V2,sep = '_',drop = TRUE)
 

interest_marks_v3 <- read.csv("interest_marks_v3.csv")

priv_sites$class = 'priv'
shar_sites$class = 'shar'

datos = rbind(priv_sites,shar_sites)
head(datos)
colnames(datos) = c("chr","snp","ref","alt","idk",'class')
head(datos)

table(datos$class)


#######Distrib chromosome
require(ggplot2)

# dens = ggplot(datos,aes(x = snp,fill=class,color=class))+
#   geom_density(alpha=.4)+
#   facet_grid(.~chr)+
#   scale_color_manual(values = colores)+
#   theme_bw()
# 
# dens
# dens +   coord_cartesian(ylim=c(0,1e-7))

head(datos)


dfr = datos

w=100000
dfr$wind <- cut(dfr$snp,breaks=seq(from=0,to=25286936+1,by=w))
head(dfr)

dfr1 = dfr %>% mutate(start=as.integer(str_extract(str_replace_all(wind,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                        end=as.integer(str_extract(str_replace_all(wind,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                        mid=start+((end-start)/2)) %>% 
  group_by(chr,mid,class) %>% 
  tally()


counts = ggplot(dfr1,aes(x=mid,y=n,color=class))+
  geom_line()+
  facet_grid(.~chr)+
  theme_minimal()

counts


#####Pattern within In(2R)Mal
head(dfr1)

dfr2 = filter(dfr1, chr==2, mid>8855602,mid< 18774479)

bkpts = c(8855602,14591030,17749309,18774475)

ggplot(dfr2, aes(x = mid,y = n,color=class))+
  geom_line()+
  geom_vline(xintercept = bkpts,linetype = 'dotted')+
  theme_minimal()


dfr3 = mutate(dfr2, inv.part = case_when(
  mid >= 17749309  & mid < 18774475 ~ "double",
  mid >= 14591030  & mid < 17749309 ~ "distal",
  mid >= 8855602 & mid < 14591030 ~ "proximal"
) )



                  
head(dfr3)

ggplot(dfr3, aes(x = mid,y = n,color=class))+
  geom_line()+
  geom_vline(xintercept = bkpts,linetype = 'dotted')+
  theme_minimal()

comp = list(c('double','distal'),c("distal","proximal"),c("double","proximal"))

ggplot(dfr3,aes(x = inv.part,y = n))+
  geom_boxplot()+
 facet_grid(.~class)+
  stat_compare_means(comparisons = comp)
  





#############Runs
head(datos)

datos = datos[order(datos$chr, datos$snp),]
head(datos)

runs = datos %>% 
  mutate(rleid = cumsum(class != lag(class, default = ""))) %>%
  group_by(chr,rleid) %>%
  dplyr::summarise(class = first(class), x = min(snp), y = max(snp), n.snps = n(), dist.xy = y-x) %>%
  select(-rleid)
runs

list.runs = split(runs,f = runs$chr)
names(list.runs) = c("a","b")

list.runs$a$dist.xx = ave(list.runs$a$x,list.runs$a$chr, FUN=function(x) c(diff(x),0))
list.runs$a$dist.next = list.runs$a$dist.xx - list.runs$a$dist.xy

list.runs$b$dist.xx = ave(list.runs$b$x,list.runs$b$chr, FUN=function(x) c(diff(x),0))
list.runs$b$dist.next = list.runs$b$dist.xx - list.runs$b$dist.xy

runs = do.call(rbind, list.runs)
runs


runs$dist.xx = NULL

#######Being:
## x = start of the run
## y = end of the run
## n.snps = number of SNPs in the run
## dist.xy = distance between first and last SNP of the cluster
## dist.next = distance to next cluster (of the opposite kind -> if P: distance to next S, if S: distance to next P)



###############################################
#####Count runs of 1 snp by categorie
dd = filter(runs,dist.xy==0)
table(dd$class,dd$chr)
inv = filter(runs, chr==2, x > 8855602, y < 18774479,dist.xy == 0)
table(inv$class,inv$chr)

####Discard runs with only one SNPs 
runs = subset(runs,runs$n.snps>1)



########Summary tables

head(runs)
#runs = subset(runs,runs$n.snps>1)


sum.runs = runs %>% group_by(chr, class) %>% 
  dplyr::summarise(mean.nsnps = mean(n.snps),sd.nsnps = sd(n.snps),
            mean.lenght = mean(dist.xy),sd.length = sd(dist.xy),
            mean.distnext = mean(dist.next),sd.distnext = sd(dist.next),
            n.runs = n())
sum.runs = as.data.frame(sum.runs)

sum.runs


#####Only in In(2R)Mal
head(runs)
inv.runs = filter(runs, chr==2, x > 8855602, y < 18774479)
write.csv(inv.runs,"snp_runs_In2RMal.csv")
head(inv.runs)


sum.inv.runs = inv.runs %>% group_by(chr, class) %>% 
  dplyr::summarise(mean.nsnps = mean(n.snps),sd.nsnps = sd(n.snps),
            mean.lenght = mean(dist.xy),sd.length = sd(dist.xy),
            mean.distnext = mean(dist.next),sd.distnext = sd(dist.next),
            n.runs = n())

sum.inv.runs = as.data.frame(sum.inv.runs)
sum.inv.runs

sum.runs$region = 'wca'
sum.inv.runs$region = 'in2rmal'

sum.runs = rbind(sum.runs,sum.inv.runs)
sum.runs



write.csv(sum.runs,"runs_summary_in2rmal.csv",row.names = FALSE)




##########
head(runs)


require(ggpubr)
require(ggforce)

in2rmal = interest_marks_v3[1,c(2,3)]

dat = subset(runs,runs$chr == '2' & runs$x>min(in2rmal) & runs$y<max(in2rmal) & runs$dist.xy>0)

distr_inv = ggplot(dat,aes(x = x,y = dist.xy,color=class))+
  geom_point(alpha=.5)+
  scale_color_manual(values = c("darkorange",'darkblue'))+
  theme_minimal()

hist_runs = ggplot(dat,aes(x = dist.xy,fill=class,color=class))+
  geom_histogram(alpha=.5,binwidth = 100,position = position_dodge())+
  scale_fill_manual(values = c("darkorange",'darkblue'))+
  scale_color_manual(values = c("darkorange",'darkblue'))+
  coord_cartesian(xlim = c(0,5000))


#######Get histogram table
brks = c(seq(from=0, to= 10000, by=1000),
         seq(from = 10000, 100000, by = 10000),
         seq(from = 100000, to = 300000, by = 100000 ))
brks = unique(brks)


dat$dist.xy_int = cut(dat$dist.xy,breaks = brks)
myhist = data.frame(table(dat$dist.xy_int,dat$class))

require(reshape2)
myhist = dcast(myhist, Var1~Var2)


write.table(myhist,'shared_priv_hist_noSingle.txt')



formula = y ~ x
runs_dist = ggplot(dat,aes(x=n.snps,y = dist.xy,color=class))+
  geom_point(alpha=.5)+
  scale_color_manual(values = c("darkorange",'darkblue'))+
  geom_smooth(method=lm, formula = formula, se=FALSE, fullrange=TRUE,alpha=.5)+
  scale_y_continuous(breaks = seq(0,300000,by = 50000),labels = c('0',"50Kb","100Kb","150Kb","200Kb","250Kb","300Kb"))+
  theme_bw()

runs_dist = runs_dist + stat_poly_eq(aes(label = paste(..rr.label..)), 
                         formula = formula, parse = TRUE, size = 3)



# row1 = plot_grid(a,labels = 'a')
# row2 = plot_grid(b,c,ncol = 2,labels = c("b","c"))
# plot_grid(row1,row2,nrow = 2)


ggscatterhist(
  dat, x = "x", y = "dist.xy",
  color = "class", size = .5, alpha = 0.6,
  palette = c("darkorange",'darkblue'),
  margin.params = list(fill = "class", color = "black", size = 0.2)
)
