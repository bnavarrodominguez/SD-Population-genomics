source("LD_decay.R")
source("shared_and_private.R")
source("recmin.R")


######ld
ld  = ld+
  labs(x="Distance (Kb)",y=expression(italic(r^2)))



######gene conversion
head(dat)

# pr1kb = filter(dat,dist.xy<1000 & class=='priv' )
# sh1kb = filter(dat,dist.xy<1000 & class=='shar' )
# table(dat$class)



hist_runs = hist_runs+
  xlab('Lenght of run (bp)')+
  ylab('Number of runs')+
  scale_fill_manual(labels=c("Privates","Shared"),values = c("darkorange","dodgerblue4"))+
  scale_color_manual(labels=c("Privates","Shared"),values = c("darkorange","dodgerblue4"))+
  theme_minimal()+
  theme(
    legend.position = c(.95, .95),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(0, 6, 6, 6),
    legend.title = element_blank()
  )

hist_runs


######crossover
lbl = sub('ZI',"SD-ZI",order_tree$ind)[-10]
lbl

hplot

hplot = hplot+
  scale_y_discrete(labels=lbl)+
  theme(axis.text.y = element_text(face = "italic"))

bin_hapl = plot_grid(hplot,rm,nrow=2,align = "v", rel_heights  = c(8,2))




require(cowplot)

r1 = plot_grid(ld,hist_runs,labels = 'AUTO',align = 'h')
r2 = plot_grid(bin_hapl,labels = 'C')
plot_grid(r1,r2,nrow = 2,align = "v")


#r1 = plot_grid(ld,hist_runs,labels = 'auto',align = 'h')
r3 = plot_grid(bin_hapl,labels = 'c')
plot_grid(r1,r3,nrow = 2,align = "v")



