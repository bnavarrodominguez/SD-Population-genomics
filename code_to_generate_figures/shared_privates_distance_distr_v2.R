#############Load libraries
require(dplyr)
require(stringr)
require(ggplot2)
require(ggforce)
require(ggpubr)
require(viridis)
require(cowplot)


priv_sites <- read.delim("shared_private_snps/zisd_private_snps.txt", header=FALSE)
shar_sites <- read.delim("shared_private_snps/zisd_shared_snps.txt", header=FALSE)


priv_sites %>% group_by(V1) %>% dplyr::summarise(tot= n())
shar_sites %>% group_by(V1) %>% dplyr::summarise(tot=n())

shar_sites$id = interaction(shar_sites$V1,shar_sites$V2,sep = '_',drop = TRUE)
priv_sites$id = interaction(priv_sites$V1,priv_sites$V2,sep = '_',drop = TRUE)


regions <- read.csv("../Figure_2_and_3/chromosome_regions.csv")

priv_sites$class = 'priv'
shar_sites$class = 'shar'

datos = rbind(priv_sites,shar_sites)
head(datos)
colnames(datos) = c("chr","snp","ref","alt","idk",'class')
head(datos)

table(datos$class)


#### Add marks
regions

datos$chromstate = ifelse(datos$chr==1 & datos$snp > regions$start[6] & datos$snp | datos$chr==2 & datos$snp < regions$end[5],"het","eu" )
datos$in2rmal = ifelse(datos$chr==2 & datos$snp > regions$start[1] &datos$snp< regions$end[1],"in2rmal",'st' )
datos$sweep = ifelse(datos$chr==1 & datos$snp > regions$start[7] & datos$snp & datos$snp < regions$end[7] 
                     | datos$chr==2 & datos$snp > regions$start[8] & datos$snp & datos$snp < regions$end[7] ,"sweep","st" )

head(datos)


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



in2rmal = regions[1,c(2,3)]

dat = subset(runs,runs$chr == '2' & runs$x>min(in2rmal) & runs$y<max(in2rmal) & runs$dist.xy>0)

distr_inv = ggplot(dat,aes(x = x,y = dist.xy,color=class))+
  geom_point(alpha=.5)+
  scale_color_manual(values = c("darkorange",'darkblue'))+
  theme_minimal()

distr_inv

hist_runs = ggplot(dat,aes(x = dist.xy,fill=class,color=class))+
  geom_histogram(alpha=.5,binwidth = 100,position = "dodge")+
  scale_fill_manual(values = c("darkorange",'darkblue'))+
  scale_color_manual(values = c("darkorange",'darkblue'))+
  coord_cartesian(xlim = c(0,10000))+
  theme()

hist_runs

