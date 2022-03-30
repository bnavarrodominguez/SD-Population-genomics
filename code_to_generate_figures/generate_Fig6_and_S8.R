###Load McClintock merged data
source("functions.R")
library(ggplot2)
library(gridExtra)
library(dplyr)
library(ggpubr)
require(readr)
require(stringr)
require(tidyverse)


############## load data

datplus = load.TEdata(file_mergedTEs = "te_calls.chr2.merged.tab3",file_TEfam = "te_fam_superfam.csv")
marks <- read_csv("../Figure_2_and_3/chromosome_regions.csv")

### Remove calls from France
datos = subset(datplus,datplus$strain%in%levels(datplus$strain)[c(3:5)],drop = TRUE)
datos=droplevels(datos)

####Remove redundancy - (TE insertions that are too close)

nr_datos = reduce_redundancy(datos,100)

dic = data.frame('strain' = levels(nr_datos$strain),'population'=c("In(2L)t",'In(2R)NS','SD-Mal'))
dic

nr_datos = merge(nr_datos,dic)

##define colors for plots
colores=c("dodgerblue4","dodgerblue","darkorange")
names(colores)=levels(nr_datos$population)


pirna = c("2R",6255432,6499291,"pirna_cluster","42AB")
marks = rbind(marks,pirna)
marks$Contig = marks$chrom
marks$start = as.numeric(marks$start)
marks$end = as.numeric(marks$end)

inv = marks %>% filter(feature == "inversion") %>%
  mutate('Contig' = chrom,
         'population' = levels(nr_datos$population)[c(3,1,2)],
         'mdpt' = (start+end)/2,
         'yaxis' = c(45,40,40))



#intpts = marks %>% filter(name %in% c("RanGAP","42AB") ) %>%
intpts = marks %>% filter(name %in% c("RanGAP") ) %>%
  mutate('Contig' = chrom,
         'mdpt' = (start+end)/2,
         'yaxis' = c(45))



sum_nr_datos =  nr_datos %>% 
  mutate(wind = cut(start,breaks = seq(0,max(nr_datos$start),by = 100000))) %>% 
  group_by(population,Contig,wind) %>%
  summarise('te_count' = dplyr:::n()) %>%
  mutate(start=as.integer(str_extract(str_replace_all(wind,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
         end=as.integer(str_extract(str_replace_all(wind,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
         mid=start+((end-start)/2))


hist3 = ggplot(sum_nr_datos,aes(x=start,y = te_count,color=population,fill=population))+
  geom_rect(data = subset(marks,marks$name=="heterochromatin"),
            inherit.aes = FALSE,aes(xmin = start, xmax = end,ymin=-Inf,ymax=+Inf),
            alpha=.2)+
  geom_vline(data = inv, aes(xintercept=start), linetype='dotted')+
  geom_vline(data = inv, aes(xintercept=end), linetype='dotted')+
  geom_segment(data = inv,aes(x = start,xend=end,y = yaxis,yend=yaxis),color='black')+
  geom_label(data = inv,aes(x = mdpt, y = yaxis,label=name),fontface = "italic",fill='white',color='black')+
  #geom_bar(stat="identity", width=.5,alpha=.5, position = position_dodge(10000))+
  geom_point(alpha=.8)+
  geom_smooth(span=.4)+
  geom_vline(data = intpts,aes(xintercept=mdpt),linetype='dotted')+
  geom_label(data = intpts,aes(x = mdpt, y = 45,label='Sd-RanGAP'),fontface = "italic", fill='white',color='black')+
  facet_grid(.~Contig,scales = 'free_x')+
  # ggtitle('Non-redundant TE calls produced by at least 2 different methods (600 bp merging window) \n
  #         Frequency considering same TE in 100bp ; binwidth=10Kb')+
  scale_fill_manual(values = colores)+
  scale_color_manual(values = colores)+
  theme_minimal(base_size = 14)

#hist3



############## Figure 6B
head(nr_datos)
colnames(nr_datos) = c('pop','family','chrom','start','end','strand','individual',
                    'insertion','Nprograms','programs','evidence','orig.ins.sites',
                    'class','subclass','poplab')
nr_datos$pop2 = sub("\\_.*", "", nr_datos$pop)

##Chromosome lenghts
clength = c("2L"=23513712,"2R"=25286936)


###Add chromatin state
het2l = min(subset(marks,marks$chrom=='2L'&marks$feature=='heterochromatin')[c(2,3)])
het2r = max(subset(marks,marks$chrom=='2R'&marks$feature=='heterochromatin')[c(2,3)])

nr_datos$chromstate = ifelse(nr_datos$chrom=='2L'&nr_datos$end<het2l | 
                            nr_datos$chrom=='2R'&nr_datos$start>het2r,
                          "eu","het")

###Add sweep region
sweep2l = subset(marks,marks$name=='2Lsweep')[c(2,3)]
sweep2r = subset(marks,marks$name=='2Rsweep')[c(2,3)]

nr_datos$region = ifelse(nr_datos$chromstate=='het',
                     'het',
                     ifelse(nr_datos$chrom=='2L'& nr_datos$start>sweep2l$start & nr_datos$end<sweep2l$end | 
                              nr_datos$chrom=='2R'& nr_datos$pop2 == 'zi' & nr_datos$start>sweep2r$start & nr_datos$end<sweep2r$end, 
                            "eu.sweep","eu.nosweep"))
head(nr_datos)



#Region lenght
marks$length = abs(marks$end-marks$start)+1
marks

reglen=data.frame('chrom' = c("2R","2R","2R",
                              "2R","2R",
                              "2L","2L","2L",
                              "2L","2L","2L"), 
                  'pop2' = c('zi','zi','zi',
                            'fr','fr',
                            'zi','zi','zi',
                            'fr','fr','fr'),
                  'region'=c('het','eu.sweep','eu.nosweep',
                             'het','eu.nosweep',
                             'het','eu.sweep','eu.nosweep',
                             'het','eu.sweep','eu.nosweep'), 
                  "rlen" = c(marks$length[5],marks$length[8],clength[2]-marks$end[8],
                             marks$length[5],clength[2]-marks$length[5],
                             marks$length[6],marks$length[7],marks$start[2]-1,
                             marks$length[6],marks$length[7],marks$start[2]-1))
reglen
head(nr_datos)


df = merge(nr_datos,reglen,by = c("pop2","chrom","region"))

head(df)



df$n = 1

detach("package:dplyr", unload=TRUE)
require('dplyr')

sum.df = df %>%
  group_by(chrom,pop2,pop,individual,region,rlen) %>%
  summarise(ins.sum = sum(n, na.rm = TRUE))
head(sum.df)
sum.df$ins.per.kb = sum.df$ins.sum/sum.df$rlen*1000
head(sum.df)

dat = sum.df
colores2 = colores
names(colores2) = levels(dat$pop)

######Normalize by the TE load in2L

####Calculate the normalization factor (number of TE insertions in euchromatin 2L distal to RanGAP)

nf = subset(df,df$chromstate!="het" & df$chrom=="2L" & df$region=="eu.nosweep")
nf

nf$n = 1

### need to detach dplyr and load again for some reason
detach("package:dplyr", unload=TRUE)
require('dplyr')

nf = nf %>%
  group_by(chrom,pop,individual) %>%
  summarise(ins.sum = sum(n, na.rm = TRUE))


nf = as.data.frame(nf)
head(nf)

nf$nf = nf$ins.sum/reglen$rlen[8]

nf


################### Normalize by 2L (whole euchromatin)
head(df)
df1 = subset(df,df$chrom=="2R" & df$chromstate=="eu")

#### All of 2R euchromatin
sumdf1 = df1 %>%
  group_by(pop,individual) %>%
  summarise(ins.sum = sum(n, na.rm = TRUE))

head(sumdf1)
rlen1 = clength[2]-marks$end[5]
sumdf1$ins.sum.bp = sumdf1$ins.sum/rlen1
head(sumdf1)
head(nf)
head(sumdf1)
dat1 = merge(sumdf1,nf,by = c("individual","pop"))
head(dat1)
dat1$region = "euchr.2R"
head(dat1)
dat1$norm.ins = dat1$ins.sum.bp/dat1$nf


#### 3) Distal from sweep

df3 = subset(df,df$chrom=="2R"&df$region=="eu.nosweep")

sumdf3 = df3 %>%
  group_by(pop,individual) %>%
  summarise(ins.sum = sum(n, na.rm = TRUE))

head(sumdf3)
rlen3 = clength[2]-marks$length[8]-marks$length[5]
sumdf3$ins.sum.bp = sumdf3$ins.sum/rlen3


head(nf)

dat3 = merge(sumdf3,nf,by = c("individual","pop"))
head(dat3)

dat3$norm.ins = dat3$ins.sum.bp/dat3$nf
dat3$region = "distal"


##### In2RMal

in2rmal = data.frame("chrom"=rep("2R",3),
                     "start"=c(14591034,8855601,8855601),
                     "stop"=c(18774475,15616195,18774475),
                     "feat"=c("distal","proximal","double"))

in2rmal$length = in2rmal$stop-in2rmal$start

df4 = subset(df,df$chrom=="2R"&df$start>in2rmal$start[3] & df$end < in2rmal$stop[3])

sumdf4 = df4 %>%
  group_by(pop,individual) %>%
  summarise(ins.sum = sum(n, na.rm = TRUE))

head(sumdf4)
rlen4 = in2rmal$length[3]
sumdf4$ins.sum.bp = sumdf4$ins.sum/rlen4


head(nf)

dat4 = merge(sumdf4,nf,by = c("individual","pop"))
head(dat4)

dat4$norm.ins = dat4$ins.sum.bp/dat4$nf
dat4$region="double.inversion"


###### plot all together

lista_dats_paper = list(dat1,dat4,dat3)

ms.datos = do.call(rbind,lista_dats_paper)
head(ms.datos)
levels(ms.datos$pop)

levels(as.factor(ms.datos$region))

ms.datos$xaxis = factor(ms.datos$region,levels = c("euchr.2R","double.inversion","distal"))

# New facet label names for supp variable
supp.labs <- c("Euchromatin", "In(2R)Mal","Distal")
names(supp.labs) <- levels(ms.datos$xaxis)


ms.boxplot = ggplot(ms.datos,aes(x = pop,y = norm.ins,fill=pop,color=pop))+
  geom_boxplot(alpha=0.5)+
  #geom_jitter(width = 0.25)+
  facet_grid(~xaxis,scales = "free",labeller = labeller(xaxis = supp.labs))+
  scale_color_manual(values = colores2)+
  scale_fill_manual(values = colores2)+
  scale_x_discrete(labels = c("Standard","SD+ In(2R)NS","SD In(2R)Mal"))+
  theme_minimal()+
  #ggtitle("Ratio 2R/2L - ratios of counts/bp")+
  xlab("Genotype")+
  ylab("Ratio 2R/2L of TE insertions/bp")

ms.boxplot

levels(ms.datos$pop)

my_comparisons = list(c("zi_dpgp2lt","zi_dpgp2rns"),c("zi_dpgp2lt","zi_sd"),c("zi_dpgp2rns","zi_sd"))
ms.boxplot = ms.boxplot+stat_compare_means(test = "kruskal.test",comparisons = my_comparisons,label = "p.signif")
#ms.boxplot

######### Format figure


hist3 = hist3 +
  scale_x_continuous(breaks = seq(0,25000000,by = 5000000),
                     labels = c("0Mb","5Mb","10Mb","15Mb","20Mb","25Mb"),
                     expand = c(0,0))+
  #scale_y_continuous(limits = c(0,50))+
  ylab("TE insertions/100kb")+
  xlab('')+
  theme(legend.position = "top",
        legend.justification = 'right',
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.text = element_text(face = 'italic'),
        strip.text.x = element_text(face = "italic"),
        legend.title = element_blank())+
  ggtitle("")


ms.boxplot = ms.boxplot+
  xlab("")+
  scale_x_discrete(labels = c("In(2L)t","In(2R)NS","SD-Mal"))+
  theme_minimal(base_size = 14)+
  theme(legend.position = "none",
        axis.text.x = element_text(face = 'italic',angle = 90,vjust=.8, hjust=0.8),
        strip.text.x = element_text(face='italic'))
ms.boxplot

# theme_get()

require(cowplot)

plot_load = plot_grid(hist3,ms.boxplot,nrow = 2,
                      labels = 'AUTO')  


plot_load



######################### Generate S8
require(dplyr)
require(tidyr)
require(ggplot2)
require(ggpubr)

head(nr_datos)

counts = nr_datos %>% 
  mutate(chrom2 = ifelse(chrom =='2R','R','L')) %>%
  group_by(class, subclass,family,pop,chrom2) %>%
  summarise(nins = n()) %>%
  spread(chrom2,nins)


zisd = filter(counts,pop=='zi_sd')
counts
counts$diff = counts$R-counts$L
counts$label = ifelse(counts$R>30,as.character(counts$family),NA)



##########Compare SD+ and SD in In(2R)Mal
#nr_datos$chrom = ifelse(nr_datos$chrom=='2R','R','L')
head(nr_datos)
str(nr_datos)
nr_datos$in2rmal = ifelse(nr_datos$chrom=='2R' & nr_datos$start > 8855602 & nr_datos$end<18774479,'inv','st')

head(nr_datos)


mytes = c("M4DM","MDG1_LTR",'LINEJ1_DM','HOBO',"ROO_I")

counts3 = nr_datos %>%
  filter(chrom == '2R' & in2rmal=='inv' |chrom == '2R' & start> 21100000 ) %>%
  group_by(class, subclass,family,pop,in2rmal,individual) %>%
  summarise(nins = n()) %>%
  ungroup() %>%
  group_by(class, subclass,family,pop,in2rmal) %>%
  summarise(aver.ins = mean(nins)) %>%
  spread(pop,aver.ins) %>%
  mutate(tot_ins = sum(zi_dpgp2lt,zi_dpgp2rns,zi_sd),diff_ins = zi_sd-zi_dpgp2lt, telabs=ifelse(family%in%mytes,as.character(family),NA)) %>%
  filter(tot_ins>1)

counts3$reg = ifelse(counts3$in2rmal=='inv',"2R:8.85-18.77","2R:18.77-25.29")
counts3$reg = factor(counts3$reg, levels = c("2R:8.85-18.77","2R:18.77-25.29"), ordered = TRUE)


pp3 = ggpaired(counts3,cond1 = 'zi_dpgp2lt',
               cond2 = 'zi_sd',
               color = "condition",
               facet.by = c('reg','class'),
               label = 'telabs',
               repel = TRUE,
               label.rectangle = TRUE,font.label = list(size = 11, face = "italic") ) +
  scale_color_manual(values = c("dodgerblue4","darkorange"))+
  scale_x_discrete(labels=c("2R","SD-Mal"))+
  xlab('Genotype')+
  ylab('TE insertions')+
  theme_bw()



pp3+theme(strip.text = element_text(face = 'italic'),
          axis.text.x = element_text(face = 'italic'),
          legend.position = 'none')



