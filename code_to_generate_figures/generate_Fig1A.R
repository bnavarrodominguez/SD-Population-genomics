library(tidyverse)
library(cowplot)


### load depth data 
temp = list.files(path = "metasembly_100bp/",pattern="*2R_meta-contig.txt",full.names = TRUE)

##do not load the recombinant SD-ZI
temp = temp[c(1:6,8:16)]
temp

myfiles = lapply(temp, function(x){read.delim(gzfile(x,'rt'),sep = '\t',header = FALSE)})

nms = gsub('metasembly_100bp//','',temp)
nms = gsub('.2R_meta-contig.txt', "",nms)
nms
names(myfiles) = nms

datos = do.call(rbind,myfiles)
head(datos)

colnames(datos) = c("contig","start",'stop','dp')

datos$line = sub("\\..*", "", rownames(datos))
head(datos)
datos$geno = ifelse(grepl('ERR',datos$line),'iso-1','sd-zi')

head(datos)

####cut to Rsp major locus
w = 180000

### average depth chrom 2 for normalization
ave.line <- read.table("metasembly_100bp/ave.depth.chrom2.txt", quote="\"", comment.char="")
colnames(ave.line) = c("line","ave.line")


datos = left_join(datos,ave.line) %>% mutate("norm.dp" = dp/ave.line)




# p = ggplot(filter(datos,stop<w),aes(x = start,y = norm.dp,color=geno))+
#   geom_point()+
#   facet_grid(.~contig,scales = 'free')+
#   coord_cartesian(ylim=c(0,10))
# 
# p

#### Load repeats location near Rsp major

repeat_bed = read.csv('meta-assembly_r6-PBcR-BLASR.V2.rsp_maj.txt',sep = '\t',header = FALSE,comment.char = "#")

head(repeat_bed)


repeat_bed = repeat_bed %>% filter(grepl(pattern = 'rsp',ignore.case = TRUE,x = V9)) %>%
  select('rsp'=V9,'contig'=V1,"start"=V4,"stop"=V5) %>%
  mutate('loc'='Rsp-major')


colnames(repeat_bed)=c('rsp','contig','start','stop','loc')
head(repeat_bed)
repeat_bed$start = as.numeric(repeat_bed$start)
repeat_bed$stop = as.numeric(repeat_bed$stop)
str(repeat_bed)

### define colors
cols = c("iso-1"="darkslategrey",'sd-zi'="darkorange")




# ggplot(datos,aes(x=start))+
#   geom_rect(data=repeat_bed, aes(xmin = start, xmax = stop), ymin = -Inf, ymax = Inf, alpha = 0.1,fill="dodgerblue")+
#   geom_point(aes(y = norm.dp,color=geno,fill=geno), alpha=.8,size=1)+
#   geom_line(aes(y = norm.dp,color=geno,fill=geno))+
#   scale_color_manual(values = cols)+
#   scale_fill_manual(values=cols)+
#   coord_cartesian(ylim=c(0,12))+
#   theme_minimal()+
#   ylab('Relative Depth')+
#   #xlab("utg_2021")+
#   theme(legend.position = 'none')


#### get average across all SD lines (exclude B6-Z1138 because low coverage in 2R)
sumdat = datos %>% 
  filter(line != 'B6-Z1138') %>%
  group_by(contig,start,stop,geno) %>%
  summarise(mean.norm.dp = mean(norm.dp),stdev=sd(norm.dp))

head(sumdat)

### plot depth in windows

dp = ggplot(filter(sumdat,stop<w),aes(x=start))+
  geom_rect(data=repeat_bed, aes(xmin = start, xmax = stop), ymin = -Inf, ymax = Inf, alpha = 0.4,fill="gray")+
  geom_ribbon(aes(ymax=mean.norm.dp+stdev,ymin=mean.norm.dp-stdev,fill=geno),alpha=.5)+
  geom_point(aes(y = mean.norm.dp,color=geno,fill=geno), alpha=.8,size=1)+
  #geom_line(aes(y = mean.norm.dp,color=geno,fill=geno))+
  scale_color_manual(values = cols,labels=c('Iso-1','SD-Mal'))+
  scale_fill_manual(values=cols,labels=c('Iso-1','SD-Mal'))+
  #coord_cartesian(ylim=c(0,15))+
  ylab('Relative Depth')+
  #xlab("utg_2021")+
  theme_minimal(base_size = 20)+
  #theme(legend.position = c(0.8, 0.2))+
  theme(legend.position = c(0.95,0.95),
        legend.title = element_blank(),
        legend.text = element_text(face='italic'))
dp

#################plot repeat coordinates from gff

source('meta_gff_plot.R')

require(cowplot)

plot_grid(dp,pgff,ncol = 1,align = "v")

dpm = dp + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())

dpm

pgff

pgff = pgff+
  xlab(label = 'Rsp Major')+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(face = 'italic'),
        axis.text.y = element_text(face = 'italic'))

pgff

plot_grid(dpm,pgff,ncol = 1,align = "v",rel_heights = c(5,4))

