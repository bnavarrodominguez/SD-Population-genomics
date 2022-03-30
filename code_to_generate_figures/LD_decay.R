require(dplyr)
require(ggplot2)
require(Rmisc)
require(stringr)
require(cowplot)


stnd <- read.csv("plink_ld_data/uninverted_2R_8.87-18.77.ld", sep="") 
zi2lt <- read.csv("plink_ld_data/in2lt.ld", sep="")
zi2rns <- read.csv("plink_ld_data/in2rns.ld", sep="")
zi.sd <- read.csv("plink_ld_data/in2rmal.ld", sep="")

stnd$pop = "st"
zi2lt$pop = "zi2lt"
zi2rns$pop = "zi2rns"
zi.sd$pop = "zisd"

datos = bind_rows(stnd,zi2lt,zi2rns,zi.sd)
head(datos)
levels(as.factor(datos$pop))

datos$dist = abs(datos$BP_A - datos$BP_B)
max(datos$dist)


dfr = select(datos,rsq=R2,dist,pop)

head(dfr)

wind = 10000

dfr$distc <- cut(dfr$dist,breaks=seq(from=min(dfr$dist)-1,to=max(dfr$dist)+1,by=wind))
head(dfr)


dfr1 = summarySE(dfr,groupvars = c("distc","pop"),measurevar = "rsq")

dat <- dfr1 %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                        end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                        mid=start+((end-start)/2))



#######My plots

####Original plot
ld_decay_plot_windows = function(dd,cols){
  #dd = dat
  require(ggplot2)
  #dd$rsq = dd$R2
  pd= position_dodge(0.1)

  ggplot(dd,aes(x = start,y = rsq,color=pop))+
    geom_smooth(se = FALSE,span=.5,alpha=.5,size=.75) +
    geom_point(position = pd,alpha=.5,size=2)+
    geom_errorbar(position = pd,aes(ymin=rsq+ci,ymax=rsq-ci))+
    labs(x="Distance (Kb)",y=expression(LD~(r^{2})))+
    scale_y_continuous(limits = c(0,1))+
    scale_x_continuous(breaks=seq(0,100000,by = 25000),
                       labels = c(0,25,50,75,100))+
    scale_color_manual(values=cols,labels = c("Standard",
                                              expression(italic("In(2L)t")),
                                              expression(italic("In(2R)NS")),
                                              expression(italic("In(2R)Mal"))
            ))+
    theme_minimal()+
    theme(legend.position = c(.9, .9),
          legend.justification = c("right", "top"),
          legend.box.just = "right",
          legend.margin = margin(6, 6, 6, 6),
          legend.text.align = 0,
          legend.title = element_blank()
          )
  
          }



cols = c('gray',"dodgerblue4","dodgerblue","darkorange")
names(cols) = levels(as.factor(dat$pop))

dat = filter(dat,N>3)

ld = ld_decay_plot_windows(dd = filter(dat,start<=100000),cols = cols)
ld
