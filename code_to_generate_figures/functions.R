load.TEdata = function(file_mergedTEs, file_TEfam){
  dat <- read.delim(file_mergedTEs)
  
  
  datplus = droplevels(subset(dat,dat$N_programs>1))
  tmp = as.data.frame(table(datplus$programs))
  ###Lots of these calls are produced by the same program twice or three times
  '%nin%' = Negate('%in%')
  
  excl = as.factor(c("ngs-te-mapper,ngs-te-mapper",
                     "ngs-te-mapper,ngs-te-mapper,ngs-te-mapper",
                     "ngs-te-mapper,ngs-te-mapper,ngs-te-mapper,ngs-te-mapper",
                     "popoolationte,popoolationte",
                     "popoolationte,popoolationte,popoolationte",
                     "popoolationte,popoolationte,popoolationte,popoolationte",
                     "retroseq,retroseq",
                     "retroseq,retroseq,retroseq",
                     "retroseq,retroseq,retroseq,retroseq",
                     'retroseq,retroseq,retroseq,retroseq,retroseq',
                     'retroseq,retroseq,retroseq,retroseq,retroseq,retroseq',
                     "temp,temp",
                     "temp,temp,temp",
                     "temp,temp,temp,temp",
                     "telocate,telocate",
                     "telocate,telocate,telocate",
                     "telocate,telocate,telocate,telocate",
                     'telocate,telocate,telocate,telocate,telocate',
                     'telocate,telocate,telocate,telocate,telocate,telocate',
                     'telocate,telocate,telocate,telocate,telocate,telocate,telocate'))
  datplus=droplevels(subset(datplus,datplus$programs%nin%excl))
  datplus = droplevels(subset(datplus,datplus$generation!="B1-C8-FR43-MP"))
  tmp = as.data.frame(table(datplus$programs))
  tmp
  te_fam_superfam <- read.csv(file_TEfam, header=FALSE)
  colnames(te_fam_superfam) = c('family','class','subclass')
  datplus = merge(datplus,te_fam_superfam,by.x = 'Family',by.y = 'family')
  return(datplus)
}
'%nin%' = Negate('%in%')

get.TEcounts = function(df,...){
  library(dplyr)
  df$n = 1
  sum.df = df %>%
    group_by_(...) %>%
    summarise(ins.sum = sum(n, na.rm = TRUE))
  head(sum.df)
  return(sum.df)
}


reduce_redundancy = function(dat,wind=100){
  dat = dat[order(dat$Family,dat$generation,dat$Contig,dat$start),]
  
  for(i in 1:(nrow(dat)-1)){
    if(dat$Family[i]==dat$Family[i+1] && dat$Contig[i]==dat$Contig[i+1] && dat$generation[i]==dat$generation[i+1] && abs(dat$start[i]-dat$start[i+1])<wind){
      dat$red_w100[i] = 'dup'
    } else {
      dat$red_w100[i] = NA
    }
  }
  
  table(dat$red_w100)
  
  dat$red_w100 = ifelse(is.na(dat$red_w100),'uniq','dup')
  table(dat$red_w100)
  
  dat = subset(dat,dat$red_w100 == 'uniq')
  dat$red_w100=NULL
  return(dat)
}


calc_freq_window = function(datos,wind){
  
  dat = datos
  dat =droplevels(dat)
  dat = dat[order(dat$Family,dat$strain,dat$Contig,dat$start),]
  #str(dat)
  
  
  dat$id = seq(1:nrow(dat))
  

  datalist = list()
  
  while(nrow(dat)>0){
    subs = subset(dat,dat$Family==dat$Family[1]&dat$Contig==dat$Contig[1]&dat$strain==dat$strain[1]&abs(dat$start-dat$start[1])<wind)
    subs = droplevels(subs)
    newline=data.frame(FAM=levels(subs$Family),CHR=levels(subs$Contig),MIN.ST=min(subs$start),POP=levels(subs$strain),CT=length(levels(subs$generation)))
    i=length(datalist)+1
    datalist[[i]]=newline
    dat = droplevels(subset(dat,dat$id%nin%subs$id))
  }
  
  out=do.call(rbind,datalist)
  
  pop_samp = unique(datos[c("strain",'generation')])
  nsamples = data.frame(table(pop_samp$strain))
  #str(nsamples)
  colnames(nsamples)=c("POP","N")
  out = merge(nsamples,out)
  out$Freq = out$CT/out$N
  
  return(out)
  
}


calc_freq = function(datos){
  
  dat = datos
  dat =droplevels(dat)
  dat = dat[order(dat$Family,dat$strain,dat$Contig,dat$start),]
  #str(dat)
  
  
  dat$id = seq(1:nrow(dat))
  
  
  datalist = list()
  
  while(nrow(dat)>0){
    subs = subset(dat,dat$Family==dat$Family[1]&dat$Contig==dat$Contig[1]&dat$strain==dat$strain[1]&dat$start==dat$start[1])
    subs = droplevels(subs)
    newline=data.frame(FAM=levels(subs$Family),CHR=levels(subs$Contig),MIN.ST=min(subs$start),POP=levels(subs$strain),CT=length(levels(subs$generation)))
    i=length(datalist)+1
    datalist[[i]]=newline
    dat = droplevels(subset(dat,dat$id%nin%subs$id))
  }
  
  out=do.call(rbind,datalist)
  
  pop_samp = unique(datos[c("strain",'generation')])
  nsamples = data.frame(table(pop_samp$strain))
  #str(nsamples)
  colnames(nsamples)=c("POP","N")
  out = merge(nsamples,out)
  out$Freq = out$CT/out$N
  
  return(out)
  
}
