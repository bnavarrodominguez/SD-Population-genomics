###Load interest parts
#marks <- read.csv("interest_marks_v2.csv")

# pop = c("zi2Lt","ziSD")
# # pop = "zi2Rns.ziSD"
# genlabels = c("WT In(2L)t", "SD In(2R)ns")
# colores = c("dodgerblue4","seagreen3")
# df = mdatos

tajdplot = function(df,pop,genlabels,colores,marks){
  require(ggplot2)
  df = subset(df,df$parameter=="TajD")
  df = subset(df,df$population%in%pop,drop = TRUE)
  df$population = factor(df$population, levels = pop)
  ggplot(df,aes(x = window_start,y=value,color=population,fill=population))+
    geom_rect(data = subset(marks,marks$name=="heterochromatin"),
              inherit.aes = FALSE,aes(xmin = start, xmax = end,ymin=-Inf,ymax=+Inf),
              alpha=.2)+
    geom_point(alpha=.5,shape=20,size=1.5)+
    geom_smooth(method = 'loess',span=.05,se = FALSE)+
    geom_hline(yintercept = 0, color='black')+
    scale_color_manual(values = colores,labels = genlabels, name="Genotype")+
    scale_fill_manual(values = colores,labels = genlabels, name="Genotype")+
    xlab("")+
    ylab("Tajima's D")+
    scale_y_continuous(limits = c(-2.5,2.5))+
    scale_x_continuous(breaks = seq(0,25000000,by = 5000000),
                       labels = c("0Mb","5Mb","10Mb","15Mb","20Mb","25Mb"),
                       expand = c(0,0))+
    facet_grid(.~chrom,scales = "free_x")+
    #ggtitle(title)+
    theme_minimal(base_size = 18)+
    theme(legend.position = 'top',plot.margin=unit(c(1,1,1,1),"cm"))
}



tajdscaledplot = function(df,pop,genlabels,colores,marks){
  require(ggplot2)
  df = subset(df,df$parameter=="TajDscaled")
  df = subset(df,df$population%in%pop,drop = TRUE)
  df$population = factor(df$population, levels = pop)
  ggplot(df,aes(x = window_start,y=value,color=population,fill=population))+
    geom_rect(data = subset(marks,marks$name=="heterochromatin"),
              inherit.aes = FALSE,aes(xmin = start, xmax = end,ymin=-Inf,ymax=+Inf),
              alpha=.2)+
    geom_point(alpha=.5,shape=20,size=1.5)+
    geom_smooth(method = 'loess',span=.05,se = FALSE)+
    geom_hline(yintercept = 0, color='black')+
    scale_color_manual(values = colores,labels = genlabels, name="Genotype")+
    scale_fill_manual(values = colores,labels = genlabels, name="Genotype")+
    xlab("")+
    ylab("D/Dmin")+
    scale_y_continuous(limits = c(-1.5,1.5))+
    scale_x_continuous(breaks = seq(0,25000000,by = 5000000),
                       labels = c("0Mb","5Mb","10Mb","15Mb","20Mb","25Mb"),
                       expand = c(0,0))+
    facet_grid(.~chrom,scales = "free_x")+
    #ggtitle(title)+
    theme_minimal(base_size = 18)+
    theme(legend.position = 'top',plot.margin=unit(c(1,1,1,1),"cm"))
}

ssplot = function(df,pop,genlabels,colores,marks){
  require(ggplot2)
  df = subset(df,df$parameter=="S")
  df = subset(df,df$population%in%pop,drop = TRUE)
  df$population = factor(df$population, levels = pop)
  
  ggplot(df,aes(x = window_start,y=value,color=population,fill=population))+
    geom_rect(data = subset(marks,marks$name=="heterochromatin"),
              inherit.aes = FALSE,aes(xmin = start, xmax = end,ymin=-Inf,ymax=+Inf),
              alpha=.2)+
    geom_point(alpha=.5,shape=20,size=1.5)+
    geom_smooth(method = 'loess',span=.05,se = FALSE)+
    scale_color_manual(values = colores,labels = genlabels, name="Genotype")+
    scale_fill_manual(values = colores,labels = genlabels, name="Genotype")+
    xlab("")+
    ylab("S")+
    scale_y_continuous(limits = c(0,600))+
    scale_x_continuous(breaks = seq(0,25000000,by = 5000000),
                       labels = c("0Mb","5Mb","10Mb","15Mb","20Mb","25Mb"),
                       expand = c(0,0))+
    facet_grid(.~chrom,scales = "free_x")+
    
    #ggtitle(title)+
    theme_minimal(base_size = 18)+
    theme(legend.position = 'top',plot.margin=unit(c(1,1,1,1),"cm"))

}

piplot = function(df,pop,genlabels,colores,marks){
  require(ggplot2)
  df = subset(df,df$parameter=="pi")
  df = subset(df,df$population%in%pop)
  df$population = factor(df$population, levels = pop)
  
  ggplot(df,aes(x = window_start,y=value,color=population,fill=population))+
    geom_rect(data = subset(marks,marks$name=="heterochromatin"),
              inherit.aes = FALSE,aes(xmin = start, xmax = end,ymin=-Inf,ymax=+Inf),
              alpha=.2)+
    geom_point(alpha=.5,shape=20,size=1.5)+
    geom_smooth(method = 'loess',span=.05,se = FALSE)+
    scale_color_manual(values = colores,labels = genlabels, name="Genotype")+
    scale_fill_manual(values = colores,labels = genlabels, name="Genotype")+
    xlab("")+
    ylab("Pi")+
    scale_y_continuous(limits = c(0,0.025))+
    scale_x_continuous(breaks = seq(0,25000000,by = 5000000),
                       labels = c("0Mb","5Mb","10Mb","15Mb","20Mb","25Mb"),
                       expand = c(0,0))+
    facet_grid(.~chrom,scales = "free_x")+
    #ggtitle(title)+
    theme_minimal(base_size = 18)+
    theme(legend.position = 'top',plot.margin=unit(c(1,1,1,1),"cm"))

}



###########Fst
fstplot = function(df,pop,genlabels,colores,marks){
  require(ggplot2)
  require(ggplot2)
  df = subset(df,df$parameter=="Fst")
  df = subset(df,df$population%in%pop)
  df$population = factor(df$population, levels = pop)
  
  ggplot(df,aes(x = window_start,y=value,color=population,fill=population))+
    geom_point(shape=20,alpha=.5,size=1.5)+
    geom_smooth(method = 'loess',span=.05,se = FALSE)+
    scale_color_manual(values = colores,labels = genlabels, name="Genotypes")+
    scale_fill_manual(values = colores,labels = genlabels, name="Genotypes")+
    xlab("Window Start")+
    ylab("Fst")+
    scale_y_continuous(limits = c(0,1.25),breaks = seq(0,1,by = 0.25))+
    scale_x_continuous(breaks = seq(0,25000000,by = 5000000),
                       labels = c("0Mb","5Mb","10Mb","15Mb","20Mb","25Mb"),
                       expand = c(0,0))+
    facet_grid(.~chrom,scales = "free_x")+
    geom_rect(data = subset(marks,marks$name=="heterochromatin"),
              inherit.aes = FALSE,aes(xmin = start, xmax = end,ymin=-Inf,ymax=+Inf),
              alpha=.2)+
    #ggtitle(title)+
    theme_minimal(base_size = 18)+
    theme(legend.position = 'top',plot.margin=unit(c(1,1,1,1),"cm"))
}



###########dxy
dxyplot = function(df,pop,genlabels,colores,marks){
  # df = mdatos
  # pop = fst_comp
  # genlabels = genlabels
  # colores = colores_pw
  # marks = landmarks
  
  #head(df)
  #levels(df$parameter)
  
  require(ggplot2)
  require(ggplot2)
  df = subset(df,df$parameter=="dxy")
  df = subset(df,df$population%in%pop)
  df$population = factor(df$population, levels = pop)
  
  ggplot(df,aes(x = window_start,y=value,color=population,fill=population))+
    geom_point(shape=20,alpha=.5,size=1.5)+
    geom_smooth(method = 'loess',span=.05,se = FALSE)+
    facet_grid(.~chrom,scales = "free_x")+
    scale_color_manual(values = colores,labels = genlabels, name="Genotypes")+
    scale_fill_manual(values = colores,labels = genlabels, name="Genotypes")+
    xlab("Window Start")+
    ylab("Dxy")+
    #scale_y_continuous(limits = c(0,1.25),breaks = seq(0,1,by = 0.25))+
    scale_x_continuous(breaks = seq(0,25000000,by = 5000000),
                       labels = c("0Mb","5Mb","10Mb","15Mb","20Mb","25Mb"),
                       expand = c(0,0))+
    geom_rect(data = subset(marks,marks$name=="heterochromatin"),
              inherit.aes = FALSE,aes(xmin = start, xmax = end,ymin=-Inf,ymax=+Inf),
              alpha=.2)+
    #ggtitle(title)+
    theme_minimal(base_size = 18)+
    theme(legend.position = 'top',plot.margin=unit(c(1,1,1,1),"cm"))
}

