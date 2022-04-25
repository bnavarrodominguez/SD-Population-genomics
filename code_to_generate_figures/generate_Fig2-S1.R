require(dplyr)
require(ggplot2)
require(reshape2)
source('windows_1kb_plots.R')

setwd("Suppl.Fig.S5/")
datos <- read.csv("ZI_1kb_windows_popgen.csv")
landmarks = read.csv("chromosome_regions.csv")


idvars = colnames(datos[c(1:10)])
idvars
mdatos = melt(datos,id.vars = idvars)

#datos <- read.csv("~/Documents/SD_SNPs/tee_pipeline/data/BOTH_populations_sepSDin_dummyGFF_10kb.out")
mdatos$chrom = ifelse(datos$chrom==1,"2L","2R")
levels(mdatos$variable)
require(stringr)
mdatos$variable = str_replace(mdatos$variable, "dxy_numsites_", "dxy.numsites_")

mdatos$population = as.factor(sub(".*_", "", mdatos$variable))
mdatos$parameter = as.factor(sub("_.*.", "", mdatos$variable))
levels(mdatos$parameter)
head(mdatos)
levels(mdatos$parameter)
levels(mdatos$population)

zambia = c("zi2Lt","zi2Rns","ziSD")
genlabels = c("ZI [In(2L)t]","ZI [In(2R)NS]", "ZI [SD-Mal]")
colores = c("dodgerblue4","dodgerblue","darkorange")
names(colores)=zambia

pi = piplot(df = mdatos,pop = zambia,genlabels = genlabels,colores = colores,marks = landmarks)+
  coord_cartesian(xlim=c(min(datos$window_start),max(datos$window_end)))+
  geom_vline(xintercept = c(19441959,19447317),linetype='dotted')+
  geom_vline(xintercept = c(19435000,19455000))


pi+ylab(expression(pi))+
  scale_x_continuous(labels = c("19.35","19.40","19.45","19.50","19.55"))+
  theme(legend.title = element_blank())




#### karyoploter
library(karyoploteR)

kp = plotKaryotype(genome="dm6",chromosomes = 'chr2L')

#zoom.region <- toGRanges(data.frame("chr2L", min(datos$window_start), max(datos$window_end)))
zoom.region <- toGRanges(data.frame("chr2L", 19435000, 19455000))
zoom.region

kp <- plotKaryotype(genome = 'dm6', chromosomes="chr2L", zoom=zoom.region)
#kpDataBackground(kp)
kpAddBaseNumbers(kp,tick.dist = 10000)

library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
dm.genes <- genes(txdb)

genes.data <- makeGenesDataFromTxDb(txdb = txdb,
                                    karyoplot=kp,
                                    plot.transcripts = TRUE, 
                                    plot.transcripts.structure = TRUE)
gdata = mergeTranscripts(genes.data = genes.data)

kpPlotGenes(kp, data=gdata,r0=0, r1=0.35)

