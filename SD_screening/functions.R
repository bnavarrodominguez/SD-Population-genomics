get.bedpe.chr2R = function(df,min.reads=10,min.dist=3000000){
  require(Sushi)
  #detach("package:dplyr", unload=TRUE)
  #detach("package:plyr", unload=TRUE)
  require(dplyr)
  require(stringr)
  
  chrom=2
  chromstart=1
  chromend=25286936
  
  df = subset(df,df$V1==chrom)
  df$dist = abs(df$V2-df$V5)
  df = subset(df,df$dist>min.dist)
  df$dist=NULL
  bedpe = df
  bedpe$V4 = ifelse(bedpe$V4=="=",bedpe$V1)
  ex.data = data("Sushi_5C.bedpe")
  #str(Sushi_5C.bedpe)
  colnames(bedpe)=colnames(Sushi_5C.bedpe)[c(1:6)]
  bedpe$name = NA
  #head(bedpe)
  # bedpe$score = log10(abs(bedpe$start1-bedpe$start2))
  
  dfr1=bedpe
  dfr1$distc = cut(bedpe$start1,include.lowest = TRUE,breaks = seq(1,chromend-200000,by = 100000))
  dfr1 <- dfr1 %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                          end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                          mid=start+((end-start)/2))
  
  dfr2=bedpe
  dfr2$distc = cut(bedpe$start2,include.lowest = TRUE,breaks = seq(1,chromend,by = 100000))
  dfr2 <- dfr2 %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                          end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                          mid=start+((end-start)/2))
  
  
  bedpe$start1 = dfr1$start
  bedpe$start2 = dfr2$start
  bedpe$end1 = dfr1$end
  bedpe$end2 = dfr2$end
  
  head(bedpe)
  bedpe2 = as.data.frame(bedpe %>% dplyr::group_by_all() %>% dplyr::summarise(score = n()))
  head(bedpe2)
  inv.coord = data.frame('start1'=c(14500000,8800000),"start2"=c(18700000,17700000))
  inv2rns = data.frame('start1'=c(15300000),"start2"=c(20200000))
  #fiveormoreeu = data.frame('start1' = c(5500000,12700000),'start2'=c(25100000,18900000))
  smol = c(12700000,18900000)
  # bedpe2$id = ifelse(bedpe2$start1 %in% inv.coord$start1 &
  #                       bedpe2$start2 %in% inv.coord$start2,
  #                     "In(2R)Mal",ifelse(bedpe2$start1<5300000,
  #                                        "heterochromatin","other"))
  bedpe2$id = ifelse(bedpe2$start1 %in% inv.coord$start1 & bedpe2$start2 %in% inv.coord$start2,
                     "In(2R)Mal",
                     ifelse(bedpe2$start1 %in% inv2rns$start1 & bedpe2$start2 %in% inv2rns$start2,
                            "In(2R)NS", "unknown"))
  
  bedpe2 = subset(bedpe2, bedpe2$score>min.reads)
  
  
  return(bedpe2)
}


my.bedpe.plot.nreads.chr2R = function(dat,mytitle="Title"){
  
  chrom=2
  chromstart=1
  chromend=25286936
  
  dat = subset(dat,dat$chrom1==chrom & dat$chrom2==chrom)
  dat$id2 = ifelse(dat$id=="In(2R)Mal",1,
                   ifelse(dat$id=="In(2R)NS",2,3))
  str(dat)
  
  # dat$id2 = if(length(unique(dat$id2))==1){
  #   runif(nrow(dat),min = 3,max = 4)
  # } else {
  #     dat$id2
  #   }
  
  str(dat)
  
  # colores=ifelse(1%in%dat$id2, colorRampPalette(c("darkorange","grey")),
  #                ifelse(2%in%dat$id2, colorRampPalette(c("blue","grey")), colorRampPalette(c("grey","grey"))))
  # 
  # colores
  # 
  #colores
  #plot.new()
  
  try(
    ifelse(
      length(unique(dat$id2))==1, 
      plotBedpe(dat,
                heights = dat$score,
                chrom = '2',
                chromstart = chromstart,
                chromend = chromend,
                plottype = 'ribbons'),
      plotBedpe(dat,
                heights = dat$score,
                chrom = '2',
                chromstart = chromstart,
                chromend = chromend,
                plottype = 'ribbons',
                colorby = dat$id2,
                colorbycol = topo.colors)
    ), silent=TRUE)
  
  
  labelgenome(chrom = "2R", chromstart,chromend,n=50,scale="Mb")
  #abline(v=c(14591034,18774475),col='black')
  #abline(v=c(8855601,15616195),col='black',lty=6)
  abline(v = c())
  #legend("topright",inset=0.01,legend=c("In(2R)Mal","In(2R)NS","heterochromatin","unknown"),col=c(topo.colors(4)),pch=19,bty='n',text.font=2)
  # legend("topright",legend = levels(as.factor(dat$id)),
  #        col = c(topo.colors(4)),
  #        pch=19,bty='n',text.font=2)
  axis(side=2,las=2,tcl=.2)
  title(mytitle)
}


get.bedpe.chr2L = function(df,min.reads=50,min.dist=3000000){
  require(Sushi)
  #detach("package:dplyr", unload=TRUE)
  #detach("package:plyr", unload=TRUE)
  require(dplyr)
  require(stringr)
  
  chrom='1'
  chromstart=1
  chromend=23513712
  
  df = subset(df,df$V1==chrom)
  df$dist = abs(df$V2-df$V5)
  df = subset(df,df$dist>min.dist)
  df$dist=NULL
  bedpe = df
  bedpe$V4 = ifelse(bedpe$V4=="=",bedpe$V1)
  ex.data = data("Sushi_5C.bedpe")
  #str(Sushi_5C.bedpe)
  colnames(bedpe)=colnames(Sushi_5C.bedpe)[c(1:6)]
  bedpe$name = NA
  #head(bedpe)
  # bedpe$score = log10(abs(bedpe$start1-bedpe$start2))
  
  dfr1=bedpe
  dfr1$distc = cut(bedpe$start1,include.lowest = TRUE,breaks = seq(1,chromend-200000,by = 100000))
  dfr1 <- dfr1 %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                          end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                          mid=start+((end-start)/2))
  
  dfr2=bedpe
  dfr2$distc = cut(bedpe$start2,include.lowest = TRUE,breaks = seq(1,chromend,by = 100000))
  dfr2 <- dfr2 %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                          end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                          mid=start+((end-start)/2))
  
  
  bedpe$start1 = dfr1$start
  bedpe$start2 = dfr2$start
  bedpe$end1 = dfr1$end
  bedpe$end2 = dfr2$end
  
  head(bedpe)
  bedpe2 = as.data.frame(bedpe %>% dplyr::group_by_all() %>% dplyr::summarise(score = n()))
  head(bedpe2)
  inv.coord = data.frame('start1'=c(2200000),"start2"=c(13100000))
  #inv2rns = data.frame('start1'=c(15300000),"start2"=c(20200000))
  #fiveormoreeu = data.frame('start1' = c(5500000,12700000),'start2'=c(25100000,18900000))
  smol = c(12700000,18900000)
  # bedpe2$id = ifelse(bedpe2$start1 %in% inv.coord$start1 &
  #                       bedpe2$start2 %in% inv.coord$start2,
  #                     "In(2R)Mal",ifelse(bedpe2$start1<5300000,
  #                                        "heterochromatin","other"))
  bedpe2$id = ifelse(bedpe2$start1 %in% inv.coord$start1 & bedpe2$start2 %in% inv.coord$start2,
                     "In(2L)t", "unknown")
  
  bedpe2 = subset(bedpe2, bedpe2$score>min.reads)
  bedpe2
  
  return(bedpe2)
}


my.bedpe.plot.nreads.chr2L = function(dat,mytitle="Title"){
  
  
  chrom='1'
  chromstart=1
  chromend=23513712
  
  dat = subset(dat,dat$chrom1==chrom & dat$chrom2==chrom)
  dat$id2 = ifelse(dat$id=="In(2L)t",1,2)
  dat
  
  # dat$id2 = if(length(unique(dat$id2))==1){
  #   runif(nrow(dat),min = 3,max = 4)
  # } else {
  #   dat$id2
  # }
  # 
  # dat
  
  #colores=ifelse(1%in%dat$id2, colorRampPalette(c("blue","grey")), colorRampPalette(c("grey","grey")))
  
  #dat
  #colores
  #plot.new()
  
  try(ifelse(
      length(unique(dat$id2))==1, 
      plotBedpe(dat,
                heights = dat$score,
                chrom = '1',
                chromstart = chromstart,
                chromend = chromend,
                plottype = 'ribbons'),
      plotBedpe(dat,
                heights = dat$score,
                chrom = '1',
                chromstart = chromstart,
                chromend = chromend,
                plottype = 'ribbons',
                colorby = dat$id2,
                colorbycol = topo.colors)
      ), 
    silent=TRUE)
  
  
  
  labelgenome(chrom = "2L", chromstart,chromend,n=50,scale="Mb")
  #abline(v=c(14591034,18774475),col='black')
  #abline(v=c(8855601,15616195),col='black',lty=6)
  #legend("topright",inset=0.01,legend=c("In(2R)Mal","In(2R)NS","heterochromatin","unknown"),col=c(topo.colors(4)),pch=19,bty='n',text.font=2)
  # legend("topright",legend = levels(as.factor(dat$id)),
  #        col = c(topo.colors(4)),
  #        pch=19,bty='n',text.font=2)
  axis(side=2,las=2,tcl=.2)
  title(mytitle)
}

get.bedpe.Sd  = function(df,min.reads=1,min.dist=3000){
  require(Sushi)
  #detach("package:dplyr", unload=TRUE)
  #detach("package:plyr", unload=TRUE)
  require(dplyr)
  require(stringr)
  
  df = subset(df,df$V1==1)
  df$dist = abs(df$V2-df$V5)
  df = subset(df,df$dist>min.dist)
  df$dist=NULL
  bedpe = df
  bedpe$V4 = ifelse(bedpe$V4=="=",bedpe$V1)
  ex.data = data("Sushi_5C.bedpe")
  #str(Sushi_5C.bedpe)
  colnames(bedpe)=colnames(Sushi_5C.bedpe)[c(1:6)]
  bedpe$name = NA
  #head(bedpe)
  # bedpe$score = log10(abs(bedpe$start1-bedpe$start2))
  
  dfr1=bedpe
  dfr1$distc = cut(bedpe$start1,include.lowest = TRUE,breaks = seq(19400000,19600000,by = 500))
  dfr1 <- dfr1 %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                          end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                          mid=start+((end-start)/2))
  
  dfr2=bedpe
  dfr2$distc = cut(bedpe$start2,include.lowest = TRUE,breaks = seq(19400000,19600000,by = 500))
  dfr2 <- dfr2 %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                          end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                          mid=start+((end-start)/2))
  
  
  bedpe$start1 = dfr1$start
  bedpe$start2 = dfr2$start
  bedpe$end1 = dfr1$end
  bedpe$end2 = dfr2$end
  
  head(bedpe)
  bedpe2 = as.data.frame(bedpe %>% dplyr::group_by_all() %>% dplyr::summarise(score = n()))
  bedpe2$name ='unknown'
  head(bedpe2)
  
  bedpe2 = subset(bedpe2,bedpe2$score>min.reads)
  bedpe2 = bedpe2[complete.cases(bedpe2),]
  
  return(bedpe2)
}

my.bedpe.plot.Sd = function(dat,mytitle="Title"){
  
  
  chrom='1'
  chromstart=19430000
  chromend=19460000
  
  # pbpe = plotBedpe(dat,1,chromstart,chromend,heights = dat$score,
  #                  plottype="loops",colorby=dat$score,colorbycol=SushiColors(3))
  pbpe = plotBedpe(dat,1,chromstart,chromend,heights = dat$score,
                   plottype="loops")
  
  
  labelgenome(chrom = "2L", chromstart,chromend,n=50,scale="bp")
  abline(v=c(19441959,19447317),col='black',lty=6)
  axis(side=2,las=2,tcl=.2)
  title(mytitle)
}
