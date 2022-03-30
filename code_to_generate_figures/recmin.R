require(dplyr)
require(stringr)
require(tidyverse)
require(ggplot2)


####Read data
snps <- read.delim("recmin/maf2_338_variants.simple.vcf")
str(snps)

in2rmal = c(18774475,14591034,15616195,8855601)


####Binary table/matrix
makebin = function(x) {ifelse(x==snps$A2.C7.Z1125,0,1)}

bin = select(snps,-MELREF, -CHROM,-num_A,-num_T,-num_G,-num_C, pos = POS)

head(bin)
bin = apply(bin[c(2:ncol(bin))],2,makebin)
bin = as.data.frame(bin)
bin$pos = snps$POS
head(bin)

pos = data.frame("pos1" = c(in2rmal[1],bin$pos))
pos$pos2 = c(bin$pos,in2rmal[2])


# ####Prepare and export for RecMin
# bin = as.matrix(bin)
# tbin = t(bin)
# colnames(tbin) = tbin[10,]
# tbin = tbin[-10,]
# 
# write.table(tbin,'bin_data_refZI125.txt',col.names = TRUE,row.names = FALSE,sep = '')


####Visualize haplotype blocks
hdat = gather(bin, key = 'ind',value = 'nt',1:9)

head(hdat)
str(hdat)

hdat$ind = sub('.*\\.', '', hdat$ind)
hdat$ind = sub('Z1',"ZI",hdat$ind)


order_tree <- read_delim("recmin/tree.txt",
"\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
colnames(order_tree) = c("ind","phylo")
order_tree = arrange(order_tree,desc(phylo))


hdat1 = left_join(hdat,order_tree,'ind') %>%
  arrange(desc(phylo)) %>%
  mutate(ind2 = factor(ind, levels=order_tree$ind))





hplot = ggplot(hdat1,aes(x = pos,y = ind2,color=nt,fill=nt ))+
  geom_line(size=5,alpha=.8)+
  xlim(c(min(in2rmal),max(in2rmal)))+
  xlab("2R")+
  ylab("SD")+
  #geom_vline(xintercept = in2rmal,linetype='dotted')+
  scale_color_gradient(low='darkorange',high = 'darkorange4')+
  theme_minimal()+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.line.x =element_blank(),
        axis.text.x=element_blank(),
        plot.margin = margin(t = 0,r = 0,b = 3,l = 0)
        )
hplot



###Add marks where recombination happened

loc <- read.csv("recmin/recmin_output.txt",header = FALSE)
str(loc)
loc$V3 = (loc$V1+loc$V2)/2

rm = ggplot(loc,aes(x = V1))+
  geom_point(y = 1,shape=2,size = 3) +
  xlim(c(min(in2rmal),max(in2rmal)))+
  ylim(c(0.9,1.1))+
  ylab("Rm")+
  xlab("2R")+
  theme_minimal()+
  #geom_vline(xintercept = in2rmal,linetype='dotted')+
  theme(axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())

rm




require(cowplot)
bin_hapl = plot_grid(hplot,rm,nrow=2,align = "v", rel_heights  = c(8,2))
bin_hapl
