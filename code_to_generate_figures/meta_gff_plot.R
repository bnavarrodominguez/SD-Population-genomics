library(ggbio)
library(GenomicRanges)
library(dplyr)

### cut in xaxis
w = 180000


#asm_gff = read.csv('asm_assembly_100bp/Supplemental_file_S11_rsp_locations.txt',sep = '\t',header = FALSE)
rsp_gff = read.csv('meta-assembly_r6-PBcR-BLASR.V2.rsp_maj.txt',sep = '\t',header = FALSE)
asm_gff = read.csv('meta-assembly_r6-PBcR-BLASR.V2.fasta.out.gff',sep='\t',header = FALSE)


repeat_bed = read.csv('meta-assembly_r6-PBcR-BLASR.V2.rsp_maj.txt',sep = '\t',header = FALSE,comment.char = "#")

head(repeat_bed)

table(repeat_bed$V1)
table(repeat_bed$V9)


repeat_bed = repeat_bed %>% filter(grepl(pattern = 'rsp',ignore.case = TRUE,x = V9)) %>%
  select('rsp'=V9,'contig'=V1,"start"=V4,"stop"=V5) %>%
  mutate('loc'='Rsp-major')


colnames(repeat_bed)=c('rsp','contig','start','stop','loc')
head(repeat_bed)
repeat_bed$start = as.numeric(repeat_bed$start)
repeat_bed$stop = as.numeric(repeat_bed$stop)
str(repeat_bed)


colnames(repeat_bed)=c('rsp','contig','start','stop','loc')
head(repeat_bed)
repeat_bed$start = as.numeric(repeat_bed$start)
repeat_bed$stop = as.numeric(repeat_bed$stop)
str(repeat_bed)

df = asm_gff %>% select("chr"=V1,'start'=V4,'end'=V5,'strand'=V7,"score"=V6,"id"=V9) %>%
  filter(chr=="2R_meta-contig_PBcR-BLASR-r6",start<w)


ggplot(df,aes(x = start,xend=end,y = id,yend=id))+
  geom_segment()+
  scale_x_reverse()

df$id2 = sub('.*:', '', df$id)
df$id2 = sub(' .*', '', df$id2)


df$id3 = ifelse(grepl('rsp',df$id2,ignore.case = TRUE),"Rsp",as.character(df$id2))


te_class = data.frame("id3"= levels(as.factor(df$id3)))
te_class$class = c("LTR","LTR","LTR","Non-LTR","Non-LTR","Non-LTR","LTR","LTR","DNA","LTR","LTR","SAT","DNA" )
te_class = te_class %>% arrange(match(class,  rev(c("SAT", "LTR", "LTR","Non-LTR","DNA"))))


df = left_join(df,te_class) %>%
  #arrange(match(class, c("SAT", "LTR", "LTR","Non-LTR","DNA"))) %>%
  mutate(id4 = factor(id3,levels = te_class$id3,ordered = TRUE))
df$id4

pgff = ggplot()+
  geom_rect(data=repeat_bed, aes(xmin = start, xmax = stop), 
            ymin = -Inf, ymax = Inf, alpha = 0.4,fill="gray")+
  geom_segment(filter(df,start<w),aes(x = start, xend=end,y = id4,yend=id4), size=3)+
  #geom_rect(data=filter(df,id3=="Rsp"), aes(xmin = start, xmax = end), ymin = -Inf, ymax = Inf, alpha = 0.4,fill="gray")+
  #scale_y_reverse()+
  theme_minimal(base_size = 14)
pgff

