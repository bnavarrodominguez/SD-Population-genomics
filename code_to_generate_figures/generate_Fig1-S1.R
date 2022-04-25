library(tidyverse)
library(cowplot)


#### load data from mosdepth

temp = list.files(path = "metasembly_rsp_locs/",pattern="*bed.gz",full.names = TRUE)

##do not load the recombinant SD-ZI
# temp
# temp = temp[c(1:6,8:16)]
# temp

# #### do not load the recombinant SD-ZI or the low coverage SD
temp
temp = temp[c(1:6,8,10:16)]
temp



myfiles = lapply(temp, function(x){read.delim(gzfile(x,'rt'),sep = '\t',header = FALSE)})

nms = gsub('metasembly_rsp_locs//','',temp)
nms = gsub('.rsp.loc.dp.regions.bed.gz', "",nms)
nms
names(myfiles) = nms

datos = do.call(rbind,myfiles)
head(datos)

colnames(datos) = c("contig","start",'stop','repeat','dp')
datos$line = sub("\\..*", "", rownames(datos))
head(datos)
datos$geno = ifelse(grepl('ERR',datos$line),'iso-1','sd-zi')

head(datos)

table(datos$contig)


#### load average depth for chromosome 2 for normalization
ave.line <- read.table("metasembly_100bp/ave.depth.chrom2.txt", quote="\"", comment.char="")
colnames(ave.line) = c("line","ave.line")

datos = left_join(datos,ave.line) %>% mutate("norm.dp" = dp/ave.line)


#### name Rsp loci
rsp_loc = ifelse(datos$contig == "2R_meta-contig_PBcR-BLASR-r6" & datos$stop <180000, 
                 "Rsp-major", 
                 ifelse(datos$contig=="X","Rsp-like",
                        ifelse(datos$contig=="3L","3L",
                               ifelse(datos$contig=='2R_meta-contig_PBcR-BLASR-r6' & datos$start >180000 & datos$stop < 24000000, "Rsp-minor", 
                                      ifelse(datos$contig=='2R_meta-contig_PBcR-BLASR-r6' & datos$start >24000000,
                                                          "60A", ifelse(datos$contig =="G2_rsp_contig_pilon", 
                                                         "Rsp-proximal","other") )))))
datos$rsp_loc = as.character(rsp_loc)



dat = datos %>% 
  group_by(geno,line,rsp_loc) %>% 
  summarise('tot.loc'=sum(norm.dp)) %>%
  filter(rsp_loc != "Rsp-like")

dat$rsp_loc = factor(dat$rsp_loc,levels = c("Rsp-proximal","Rsp-major","Rsp-minor","60A","3L"),ordered = TRUE)

cols = c("iso-1"="darkslategrey",'sd-zi'="darkorange")

ploc = ggplot(dat,aes(x=rsp_loc,y = tot.loc,color=geno))+
  geom_boxplot()+
  scale_color_manual(values = cols,labels=c("Iso-1","ZI[SD-Mal]"))+
  theme_minimal(base_size = 14)

ploc = ploc+
  ylab('Total Depth')+
  #xlab(label = "Rsp locus")+
  theme_bw(base_size = 14)+
  theme(legend.position = c(0.9, 0.9),
        #legend.position = "left",
        axis.text.x = element_text(face = "italic"),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(face='italic'))

ploc

