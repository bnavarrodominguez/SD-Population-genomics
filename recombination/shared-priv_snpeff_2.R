library(dplyr)
library(tidyr)
library(readr)


#####generated with script snpeff_snpsift_filtered_snps.R
snpeff_sd <- read.csv("snpeff_zisd_snpsift_filtered.txt")
snpeff_zi2lt <- read.csv("snpeff_zi2Lt_snpsift_filtered.txt")
snpeff_zi2rns <- read.csv("snpeff_zi2rns_snpsift_filtered.txt")

head(snpeff)
snpeff_sd$single_ef =  ifelse(grepl("missense_variant",snpeff_sd$V9),
                           "missense_variant", "synonymous_variant")
snpeff_zi2lt$single_ef =  ifelse(grepl("missense_variant",snpeff_zi2lt$V9),
                              "missense_variant", "synonymous_variant")
snpeff_zi2rns$single_ef =  ifelse(grepl("missense_variant",snpeff_zi2rns$V9),
                              "missense_variant", "synonymous_variant")


####load shared and private data
files = list.files("data_density_varpos_nohet_mac1/",pattern = "*.fromVCF.txt",full.names = TRUE)
snps_class = lapply(files,function(x){read.delim(x,header = FALSE)})
df_names = gsub("data_density_varpos_nohet_mac1//","",files)
df_names = gsub("_hardfilters2.sites.fromVCF.txt","",df_names)
names(snps_class) = df_names

snps_class$zi2lt_private$V5 = 'priv'
snps_class$zi2lt_sharedby2$V5 = 'shar'
snps_class$zi2rns_private$V5 = 'priv'
snps_class$zi2rns_sharedby2$V5 = 'shar'
snps_class$zisd_private$V5 = 'priv'
snps_class$zisd_shared$V5 = 'shar'

snps_class = lapply(snps_class, setNames, c('chrom','pos','ref','alt','class'))


zi2lt_class = bind_rows(snps_class[c(1,2)]) %>%
  mutate('chrom' = ifelse(chrom=='1','2L','2R'))

zi2rns_class = bind_rows(snps_class[c(3,4)]) %>%
  mutate('chrom' = ifelse(chrom=='1','2L','2R'))

zisd_class = bind_rows(snps_class[c(5,6)]) %>%
  mutate('chrom' = ifelse(chrom=='1','2L','2R'))

####Shared & private SNPs in In(2R)Mal
shr_priv_in2rmal = filter(zisd_class, chrom=='2R' & 
                            pos>min(in2rmal_int_brkpts$X2) & 
                            pos<max(in2rmal_int_brkpts$X3))
table(shr_priv_in2rmal$class)



###rename and subset SNPeff info
head(snpeff_sd)


snpeff_zi2lt = snpeff_zi2lt %>%
  select('chrom' = V1, 'pos' = V2, 'ref' = V4, 'alt' = V5, 'freq' = maxAF, 'eff' = single_ef, 'gene' = V12)
snpeff_zi2rns = snpeff_zi2rns %>% 
  select('chrom' = V1, 'pos' = V2, 'ref' = V4, 'alt' = V5, 'freq' = maxAF, 'eff' = single_ef, 'gene' = V12)
snpeff_sd = snpeff_sd %>% 
  select('chrom' = V1, 'pos' = V2, 'ref' = V4, 'alt' = V5, 'freq' = maxAF, 'eff' = single_ef, 'gene' = V12)

#####merge
head(snpeff_sd)
head(zi2rns_class)

dat_zi2lt = left_join(snpeff_zi2lt,zi2lt_class,c("chrom","pos","ref","alt")) %>%
  mutate('pop' = 'zi2lt')
dat_zi2rns = left_join(snpeff_zi2rns,zi2rns_class,c("chrom","pos","ref","alt")) %>%
  mutate('pop' = 'zi2rns')
dat_zisd = left_join(snpeff_sd,zisd_class,c("chrom","pos","ref","alt")) %>%
  mutate('pop' = 'zisd')

dat_all = bind_rows(dat_zi2lt,dat_zi2rns,dat_zisd)


head(dat_all)
dat_all = dat_all[complete.cases(dat_all),]


sum_all = dat_all %>% mutate('eff' = ifelse(eff=='missense_variant', 'nsyn','syn'),
                   'class2' = interaction(eff,class)) %>% 
  group_by(chrom,gene,pop) %>%
  summarise('pos' = mean(pos),
            'nsyn.shar' = sum(class2 == 'nsyn.shar'),
            'syn.shar'= sum(class2 == 'syn.shar'),
            'nsyn.priv' = sum(class2 == 'nsyn.priv'),
            'syn.priv'= sum(class2 == 'syn.priv'))




sum_all2 = sum_all %>% mutate('ns_shar' = nsyn.shar/syn.shar, 'ns_priv' = nsyn.priv/syn.priv) %>%
  select(chrom,gene,pop,pos,ns_shar,ns_priv) %>%
  gather(class,ns,-chrom,-gene,-pop,-pos) %>%
  filter(ns!=Inf)


library(ggplot2)
ggplot(sum_all2,aes(x = pos,y = ns,color=class))+
  geom_point(alpha=.5)+
  facet_grid(pop~chrom)+
  theme_bw()

head(sum_all)

sumall_genes = sum_all %>% 
  group_by()
  mutate(ns_overall = (nsyn.shar + nsyn.priv)/(syn.shar + syn.priv)) %>%
  filter(ns_overall != Inf) 

ggplot(sumall_genes,aes(x = pos,y = ns_overall))+
  geom_point(alpha=.5)+
  facet_grid(pop~chrom)+
  theme_bw()



##############All different Regions

marks <- read_csv("~/Documentos/investigacion/ur_larracuente_presgraves/SD_SNPs/interest_marks_v3.csv")
marks
head(dat_all)

dat_all2 = dat_all %>% mutate('reg'= ifelse(chrom=='2L' & pos > 22000975 | chrom=='2R' & pos < 5398184, 'het',
                     ifelse(chrom=='2L' & pop=='zi2lt' & pos > 2225744 & pos < 13154180, 'in2lt',
                            ifelse(chrom=='2R' & pop=='zi2rns' & pos > 15391154 & pos < 20276334, 'in2rns',
                                   ifelse(chrom=='2R' & pop=='zisd' & pos > 8855602 & pos < 18774479, 'in2rmal',
                                          'st')))))
#check it 
ggplot(dat_all2,aes(x = pos,y = pop,color=reg))+
  geom_point()+
  facet_grid(.~chrom)

head(dat_all2)

sumreg = dat_all2 %>% 
  mutate(eff = ifelse(eff == 'missense_variant','nsyn','syn'),
         class2 = interaction(class,eff)) %>%
  group_by(chrom,class2,pop,reg) %>%
  summarise(n = n()) %>%
  spread(key = class2,value = n) %>%
  mutate(ns_shared = shar.nsyn/shar.syn, ns_priv = priv.nsyn/priv.syn)


##############Inversions and same uninverted region
marks
head(dat_all)

dat_all3 = dat_all %>% 
  mutate('chromstate'= ifelse(chrom=='2L' & pos > 22000975 | chrom=='2R' & pos < 5398184, 'het', 'eu'), 
         'in2lt' = ifelse(chrom=='2L' &  pos > 2225744 & pos > 13154180, 'in2lt','st'),
         'in2rmal' = ifelse(chrom=='2R' & pos > 8855602 & pos < 18774479, 'in2rmal', 'st'),
         'in2rns' = ifelse(chrom=='2R' & pos > 15391154 & pos < 20276334, 'in2rns', 'st'))

sum_in2rmal = dat_all3 %>% 
  filter(chromstate == 'eu') %>%
  mutate(eff = ifelse(eff == 'missense_variant','nsyn','syn'),
         class2 = interaction(class,eff)) %>%
  group_by(chrom,class2,pop,in2rmal) %>%
  summarise(n = n()) %>%
  spread(key = class2,value = n) %>%
  mutate(ns_shared = shar.nsyn/shar.syn, ns_priv = priv.nsyn/priv.syn)

sum_in2lt = dat_all3 %>% 
  filter(chromstate == 'eu') %>%
  mutate(eff = ifelse(eff == 'missense_variant','nsyn','syn'),
         class2 = interaction(class,eff)) %>%
  group_by(chrom,class2,pop,in2lt) %>%
  summarise(n = n()) %>%
  spread(key = class2,value = n) %>%
  mutate(ns_shared = shar.nsyn/shar.syn, ns_priv = priv.nsyn/priv.syn)

sum_in2rns = dat_all3 %>% 
  filter(chromstate == 'eu') %>%
  mutate(eff = ifelse(eff == 'missense_variant','nsyn','syn'),
         class2 = interaction(class,eff)) %>%
  group_by(chrom,class2,pop,in2rns) %>%
  summarise(n = n()) %>%
  spread(key = class2,value = n) %>%
  mutate(ns_shared = shar.nsyn/shar.syn, ns_priv = priv.nsyn/priv.syn)


sum_in2lt
sum_in2rmal
sum_in2rns


####Regions in In(2R)Mal
in2rmal_int_brkpts <- read_delim("~/Documentos/investigacion/ur_larracuente_presgraves/SD_SNPs/ld/in2rmal_int_brkpts",
"\t", escape_double = FALSE, col_names = FALSE,trim_ws = TRUE)

head(dat_all3)
df = dat_all3 %>% filter(in2rmal == 'in2rmal',pop=='zisd') %>% mutate(eff2 = ifelse(eff=='missense_variant','nsyn','syn')) 
df1 = df %>% select(pos,class,eff2) %>% gather(variable, value, -pos,class,eff2)


head(df1)


ggplot(df1,aes(x = pos,fill=value))+
  geom_density(alpha=.5)+
  #facet_grid(variable~.,scales = 'free_y')+
  geom_vline(xintercept = in2rmal_int_brkpts$X3,linetype='dotted')+
  geom_vline(xintercept = in2rmal_int_brkpts$X2,linetype='dotted')+
  theme_minimal()+
  theme(legend.position = 'right')

ggplot(df1,aes(x = pos,y = value,fill=value))+
  geom_violin(alpha=.5)+
  facet_grid(variable~.,scales = 'free_y')+
  geom_vline(xintercept = in2rmal_int_brkpts$X3,linetype='dotted')+
  geom_vline(xintercept = in2rmal_int_brkpts$X2,linetype='dotted')+
  theme_minimal()+
  theme(legend.position = 'none')

df2 = df %>% select(pos,class,freq,eff2) %>% gather(variable, value, -pos,-freq,class,eff2)

ggplot(df2,aes(x = pos,y = freq,fill=value,color=value))+
  geom_point()+
  geom_smooth(alpha=.5)+
  facet_grid(variable~.,scales = 'free_y')+
  geom_vline(xintercept = in2rmal_int_brkpts$X3,linetype='dotted')+
  geom_vline(xintercept = in2rmal_int_brkpts$X2,linetype='dotted')+
  theme_minimal()+
  theme(legend.position = 'left')



#############Chisquared -test

sum_in2rmal


dt = ungroup(sum_in2rmal) %>% 
  filter(pop !='zi2rns' & in2rmal=='in2rmal') %>%
  mutate(geno = c('wt','sd')) %>%
  select(geno,priv.nsyn,shar.nsyn,priv.syn,shar.syn)



##1) All SNPs
dt

wt_ns =  dt$priv.nsyn[1]+dt$shar.nsyn[1]
wt_s = dt$priv.syn[1]+dt$shar.syn[1]
sd_ns = dt$priv.nsyn[2]+dt$shar.nsyn[2]
sd_s = dt$priv.syn[2]+dt$shar.syn[2]

tab = matrix(c(wt_ns,sd_ns,wt_s,sd_s),nrow = 2)
tab

chisq.test(tab,correct = FALSE)

###goodness of fit
sd = c(sd_ns,sd_s)
wt = prop.table(c(wt_ns,wt_s))

wt
chisq.test(x = sd,p = wt)

###p<0.05 --> not same distribution

##2) shared snps
shar = as.matrix(dt[,c(3,5)],nrow = 2)

chisq.test(shar)

sd = shar[2,]
sd
wt = prop.table(shar[1,])
wt

chisq.test(x = sd,p = wt)

###p > 0.05 same distribution

##3) priv snps
dt
priv = as.matrix(dt[,c(2,4)],nrow = 2)
priv

sd = priv[2,]
sd
wt = prop.table(priv[1,])
wt

chisq.test(x = sd,p = wt)

###p < 0.05 different distribution

#######Differences between shared and priv in SD+
mtx = matrix(dt[1,c(2,3,4,5)],nrow = 2)
mtx = unlist(mtx)
mtx = matrix(mtx,nrow = 2)
mtx
chisq.test(mtx)
