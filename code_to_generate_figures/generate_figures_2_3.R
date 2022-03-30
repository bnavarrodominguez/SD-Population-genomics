#### Load packages
require(cowplot)
require(ggplot2)
require(reshape2)
require(stringr)

####Source file with plot codes
source("popgen_windows.R")

###Load interest parts
marks <- read.csv("chromosome_regions.csv")
marks

### Load data
datos = read.csv("ZI_min8_filtered_dummyGFF_10kb.v2.out")

##### scale dxy (per nucleotide)
datos$dxy_zi2Lt.zi2Rns = datos$dxy_zi2Lt.zi2Rns/datos$dxy_numsites_zi2Lt.zi2Rns
datos$dxy_zi2Lt.ziSD = datos$dxy_zi2Lt.ziSD/datos$dxy_numsites_zi2Lt.ziSD
datos$dxy_zi2Rns.ziSD = datos$dxy_zi2Rns.ziSD/datos$dxy_numsites_zi2Rns.ziSD

require(reshape2)
head(datos)

idvars = colnames(datos[c(1:10)])
idvars
mdatos = melt(datos,id.vars = idvars)

mdatos$chrom = ifelse(datos$chrom==1,"2L","2R")
levels(mdatos$variable)

mdatos$variable = str_replace(mdatos$variable, "dxy_numsites_", "dxy.numsites_")

mdatos$population = as.factor(sub(".*_", "", mdatos$variable))
mdatos$parameter = as.factor(sub("_.*.", "", mdatos$variable))
levels(mdatos$parameter)
head(mdatos)
levels(mdatos$parameter)
levels(mdatos$population)

zambia = c("zi2Lt","zi2Rns","ziSD")
genlabels = c("In(2L)t","In(2R)NS", "SD-Mal")
colores = c("dodgerblue4","dodgerblue","darkorange")
names(colores)=zambia

# piplot(df = mdatos,pop = zambia,genlabels = genlabels,colores = colores,marks = marks)
# tajdplot(df = mdatos,pop = zambia,genlabels = genlabels,colores = colores,marks = marks)
# tajdscaledplot(df = mdatos,pop = zambia,genlabels = genlabels,colores = colores,marks = marks)



td = tajdplot(df = mdatos,pop = zambia,genlabels = genlabels,colores = colores,marks = marks)+
  geom_vline(data = subset(marks,marks$name=="RanGAP"),aes(xintercept = start),color ="black",linetype="dotted")+
  geom_label(data = subset(marks,marks$name=="RanGAP"), aes(x=(start+end)/2,y = 2.5,label="Sd-RanGAP"),fontface = "italic",color="black",inherit.aes = FALSE)+
  geom_segment(data = subset(marks,marks$name=="In(2R)Mal"),aes(x = start, xend = end,y = 2.5,yend=2.5),color = "black",inherit.aes = FALSE)+
  geom_label(data = subset(marks,marks$name=="In(2R)Mal"), aes(x=(start+end)/2,y = 2.5,label=name),fontface = "italic",color="black",inherit.aes = FALSE)+
  geom_segment(data = subset(marks,marks$name=="In(2R)NS"),aes(x = start, xend = end,y = 2,yend=2),color = "black",inherit.aes = FALSE)+
  geom_label(data = subset(marks,marks$name=="In(2R)NS"), aes(x=(start+end)/2,y = 2,label=name),fontface = "italic",color="black",inherit.aes = FALSE)+
  geom_segment(data = subset(marks,marks$name=="In(2L)t"),aes(x = start, xend = end,y = 2,yend=2),color = "black",inherit.aes = FALSE)+
  geom_label(data = subset(marks,marks$name=="In(2L)t"), aes(x=(start+end)/2,y = 2,label=name),fontface = "italic",color="black",inherit.aes = FALSE)+
  geom_vline(data = subset(marks,marks$feature == 'inversion'), aes(xintercept=start),linetype="dotted",color="darkgray") +
  geom_vline(data = subset(marks,marks$feature == 'inversion'), aes(xintercept=end),linetype="dotted",color="darkgray")

td


tde = tajdscaledplot(df = mdatos,pop = zambia,genlabels = genlabels,colores = colores,marks = marks)+
  geom_vline(data = subset(marks,marks$name=="RanGAP"),aes(xintercept = start),color ="black",linetype="dotted")+
  geom_label(data = subset(marks,marks$name=="RanGAP"), aes(x=(start+end)/2,y = 1.5,label="Sd-RanGAP"),fontface = "italic",color="black",inherit.aes = FALSE)+
  geom_segment(data = subset(marks,marks$name=="In(2R)Mal"),aes(x = start, xend = end,y = 1.5,yend=1.5),color = "black",inherit.aes = FALSE)+
  geom_label(data = subset(marks,marks$name=="In(2R)Mal"), aes(x=(start+end)/2,y = 1.5,label=name),fontface = "italic",color="black",inherit.aes = FALSE)+
  geom_segment(data = subset(marks,marks$name=="In(2R)NS"),aes(x = start, xend = end,y = 1,yend = 1),color = "black",inherit.aes = FALSE)+
  geom_label(data = subset(marks,marks$name=="In(2R)NS"), aes(x=(start+end)/2,y = 1,label=name),fontface = "italic",color="black",inherit.aes = FALSE)+
  geom_segment(data = subset(marks,marks$name=="In(2L)t"),aes(x = start, xend = end,y = 1,yend = 1),color = "black",inherit.aes = FALSE)+
  geom_label(data = subset(marks,marks$name=="In(2L)t"), aes(x=(start+end)/2,y = 1,label=name),fontface = "italic",color="black",inherit.aes = FALSE)+
  geom_vline(data = subset(marks,marks$feature == 'inversion'), aes(xintercept=start),linetype="dotted",color="darkgray") +
  geom_vline(data = subset(marks,marks$feature == 'inversion'), aes(xintercept=end),linetype="dotted",color="darkgray")

pi = piplot(df = mdatos,pop = zambia,genlabels = genlabels,colores = colores,marks = marks)+
  geom_vline(data = subset(marks,marks$name=="RanGAP"),aes(xintercept = start),color ="black",linetype="dotted")+
  geom_label(data = subset(marks,marks$name=="RanGAP"), aes(x=(start+end)/2,y = 0.025,label='Sd-RanGAP'),fontface = "italic",color="black",inherit.aes = FALSE)+
  geom_segment(data = subset(marks,marks$name=="In(2R)Mal"),aes(x = start, xend = end,y = 0.025,yend=0.025),color = "black",inherit.aes = FALSE)+
  geom_label(data = subset(marks,marks$name=="In(2R)Mal"), aes(x=(start+end)/2,y = 0.025,label=name),fontface = "italic",color="black",inherit.aes = FALSE)+
  geom_segment(data = subset(marks,marks$name=="In(2R)NS"),aes(x = start, xend = end,y = 0.0225,yend=0.0225),color = "black",inherit.aes = FALSE)+
  geom_label(data = subset(marks,marks$name=="In(2R)NS"), aes(x=(start+end)/2,y = 0.0225,label=name),fontface = "italic",color="black",inherit.aes = FALSE)+
  geom_segment(data = subset(marks,marks$name=="In(2L)t"),aes(x = start, xend = end,y = 0.0225,yend=0.0225),color = "black",inherit.aes = FALSE)+
  geom_label(data = subset(marks,marks$name=="In(2L)t"), aes(x=(start+end)/2,y = 0.0225,label=name),fontface = "italic",color="black",inherit.aes = FALSE)+
  geom_vline(data = subset(marks,marks$feature == 'inversion'), aes(xintercept=start),linetype="dotted",color="darkgray") +
  geom_vline(data = subset(marks,marks$feature == 'inversion'), aes(xintercept=end),linetype="dotted",color="darkgray")
pi

ss = ssplot(df = mdatos,pop = zambia,genlabels = genlabels,colores = colores,marks = marks)+
  geom_vline(data = subset(marks,marks$name=="RanGAP"),aes(xintercept = start),color ="black",linetype="dotted")+
  geom_label(data = subset(marks,marks$name=="RanGAP"), aes(x=(start+end)/2,y = 600,label=name),color="black",inherit.aes = FALSE)+
  geom_segment(data = subset(marks,marks$name=="In(2R)Mal"),aes(x = start, xend = end,y = 600,yend=600),color = "black",inherit.aes = FALSE)+
  geom_label(data = subset(marks,marks$name=="In(2R)Mal"), aes(x=(start+end)/2,y = 600,label=name),color="black",inherit.aes = FALSE)+
  geom_segment(data = subset(marks,marks$name=="In(2R)NS"),aes(x = start, xend = end,y = 550,yend=550),color = "black",inherit.aes = FALSE)+
  geom_label(data = subset(marks,marks$name=="In(2R)NS"), aes(x=(start+end)/2,y = 550,label=name),color="black",inherit.aes = FALSE)+
  geom_segment(data = subset(marks,marks$name=="In(2L)t"),aes(x = start, xend = end,y = 550,yend=550),color = "black",inherit.aes = FALSE)+
  geom_label(data = subset(marks,marks$name=="In(2L)t"), aes(x=(start+end)/2,y = 550,label=name),color="black",inherit.aes = FALSE)+
  geom_vline(data = subset(marks,marks$feature == 'inversion'), aes(xintercept=start),linetype="dotted",color="darkgray") +
  geom_vline(data = subset(marks,marks$feature == 'inversion'), aes(xintercept=end),linetype="dotted",color="darkgray")

ss





levels(mdatos$population)
#dfst = subset(mdatos,mdatos$chrom=="2L"&mdatos$population=='zi2Rns.ziSD'|mdatos$chrom=="2R"&mdatos$population=="zi2Lt.ziSD")
#dfst$population = "WT.SD"

fst_comp = c("zi2Lt.zi2Rns","zi2Lt.ziSD","zi2Rns.ziSD")
genlabels_pw = c("In(2L)t vs. In(2R)NS","In(2L)t vs. SD-Mal","In(2R)NS vs. SD-Mal" )

#require(viridis)
#colores_pw = inferno(6)
#barplot(1:length(colores_pw),col=colores_pw)
#colores_pw = c(colores_pw[4],colores_pw[2],colores_pw[5])
colores_pw = c('darkgrey',"#AE311E","#EDA096")
names(colores_pw)=c("zi2Lt.zi2Rns","zi2Lt.ziSD","zi2Rns.ziSD")


fst = fstplot(df = mdatos,pop = fst_comp,genlabels = genlabels_pw,colores = colores_pw,marks = marks)+
  geom_vline(data = subset(marks,marks$name=="RanGAP"),aes(xintercept = start),color ="black",linetype="dotted")+
  geom_label(data = subset(marks,marks$name=="RanGAP"), aes(x=(start+end)/2,y = 1.125,label='Sd-RanGAP'),fontface = "italic",color="black",inherit.aes = FALSE)+
  geom_segment(data = subset(marks,marks$name=="In(2R)Mal"),aes(x = start, xend = end,y = 1.125,yend=1.125),color = "black",inherit.aes = FALSE)+
  geom_label(data = subset(marks,marks$name=="In(2R)Mal"), aes(x=(start+end)/2,y = 1.125,label=name),fontface = "italic",color="black",inherit.aes = FALSE)+
  geom_segment(data = subset(marks,marks$name=="In(2R)NS"),aes(x = start, xend = end,y = 1,yend=1),color = "black",inherit.aes = FALSE)+
  geom_label(data = subset(marks,marks$name=="In(2R)NS"), aes(x=(start+end)/2,y = 1,label=name),fontface = "italic",color="black",inherit.aes = FALSE)+
  geom_segment(data = subset(marks,marks$name=="In(2L)t"),aes(x = start, xend = end,y = 1,yend=1),color = "black",inherit.aes = FALSE)+
  geom_label(data = subset(marks,marks$name=="In(2L)t"), aes(x=(start+end)/2,y = 1,label=name),fontface = "italic",color="black",inherit.aes = FALSE)+
  geom_vline(data = subset(marks,marks$feature == 'inversion'), aes(xintercept=start),linetype="dotted",color="darkgray") +
  geom_vline(data = subset(marks,marks$feature == 'inversion'), aes(xintercept=end),linetype="dotted",color="darkgray") 
fst


#######Dxy

dd = droplevels(subset(mdatos,mdatos$parameter=='dxy'))
dd = droplevels(subset(dd,dd$population%in%fst_comp))
levels(dd$population)

fst_comp


dxy = dxyplot(df = mdatos,pop = fst_comp,genlabels = genlabels_pw,colores = colores_pw,marks = marks)
dxy = dxy+
  geom_vline(data = subset(marks,marks$name=="RanGAP"),aes(xintercept = start),color ="black",linetype="dotted")+
  geom_label(data = subset(marks,marks$name=="RanGAP"), aes(x=(start+end)/2,y = 0.022,label="Sd-RanGAP"),fontface = "italic",color="black",inherit.aes = FALSE)+
  geom_segment(data = subset(marks,marks$name=="In(2R)Mal"),aes(x = start, xend = end,y = 0.022,yend=0.022),color = "black",inherit.aes = FALSE)+
  geom_label(data = subset(marks,marks$name=="In(2R)Mal"), aes(x=(start+end)/2,y = 0.022,label=name),fontface = "italic",color="black",inherit.aes = FALSE)+
  geom_segment(data = subset(marks,marks$name=="In(2R)NS"),aes(x = start, xend = end,y = 0.02,yend=0.02),color = "black",inherit.aes = FALSE)+
  geom_label(data = subset(marks,marks$name=="In(2R)NS"), aes(x=(start+end)/2,y = 0.02,label=name),fontface = "italic",color="black",inherit.aes = FALSE)+
  geom_segment(data = subset(marks,marks$name=="In(2L)t"),aes(x = start, xend = end,y = 0.02,yend=0.02),color = "black",inherit.aes = FALSE)+
  geom_label(data = subset(marks,marks$name=="In(2L)t"), aes(x=(start+end)/2,y = 0.02,label=name),fontface = "italic",color="black",inherit.aes = FALSE)+
  geom_vline(data = subset(marks,marks$feature == 'inversion'), aes(xintercept=start),linetype="dotted",color="darkgray") +
  geom_vline(data = subset(marks,marks$feature == 'inversion'), aes(xintercept=end),linetype="dotted",color="darkgray") 
dxy

