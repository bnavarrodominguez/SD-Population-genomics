library('Sushi')
source('functions.R')


#### Detect Sd-RanGAP duplication

temp = list.files(pattern="*.bedpe",path = 'bedpe_Sd/',full.names = TRUE)
temp = as.list(temp)

myfiles = lapply(temp, function(x){read.delim(x, header = FALSE,stringsAsFactors = FALSE) })
nms = gsub(x = temp,pattern = ".sel.rangap.bedpe",replacement = "")
nms = gsub(x = nms,pattern = 'bedpe_Sd//',replacement = "")
names(myfiles) = nms


mybedpes = lapply(myfiles,function(x){get.bedpe.Sd(x,min.reads = 1,min.dist = 3000)})

titles = names(mybedpes)

pdf("screening_Sd.pdf",width = 12,height = 6,onefile = TRUE)
par(mfrow = c(2,2))
for(i in 1:length(mybedpes)){
  my.bedpe.plot.Sd(mybedpes[[i]],mytitle = titles[i])
  #plot.new()
}

dev.off()


#### Detect 2L and 2R inversions

library('Sushi')

#source('functions.R')

temp = list.files(pattern="*.bedpe",path = 'bedpe_inversions/',full.names = TRUE)
myfiles = lapply(temp, function(x){read.delim(x, header = FALSE) })

nms = gsub(x = temp,pattern = ".chr2.sel.3Mb.bedpe",replacement = "")
nms = gsub(nms,pattern = 'bedpe_inversions//',replacement = '')
names(myfiles) = nms



########Chromosome 2R

mybedpes = lapply(myfiles,function(x){get.bedpe.chr2R(x,min.reads = 5,min.dist = 3000000)})
titles = names(mybedpes)


pdf("screening_chr2R_inversion.pdf",width = 12,height = 8,onefile = TRUE)
par(mfrow=c(2,2))
for(i in 1:length(mybedpes)){
  my.bedpe.plot.nreads.chr2R(mybedpes[[i]],mytitle = titles[i])
  #plot.new()
}
dev.off()



############Chromosome 2L 


mybedpes = lapply(myfiles,function(x){get.bedpe.chr2L(x,min.reads = 5,min.dist = 3000000)})

titles = names(mybedpes)

pdf("screening_chr2L_inversion.pdf",width = 12,height = 8,onefile = TRUE)
par(mfrow=c(2,2))
for(i in 1:length(mybedpes)){
  my.bedpe.plot.nreads.chr2L(mybedpes[[i]],mytitle = titles[i])
  #plot.new()
}
dev.off()



