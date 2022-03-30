library(tidyverse)
library(cowplot)


##############Goodness of fit (all non-coding snps)
emp <- read_delim("../Figure_4/simulations_data/sumst_ZI_in2rmal_overall_noncoding.txt", 
                  " ", escape_double = FALSE, col_names = FALSE, 
                  trim_ws = TRUE)
colnames(emp) = c("pi.m","s.m","tajD.m")
emp



### load sweep

sw = read_delim("simulations_models/sweep_S.2237_t.0.088400.txt", " ",escape_double = FALSE, col_names = c('pi','s','tajD'), 
                comment = "#", trim_ws = TRUE)
sw$model = "sweep"

const = read_delim("simulations_models/const_S.2237.txt", " ",escape_double = FALSE, col_names = c('pi','s','tajD'), 
                comment = "#", trim_ws = TRUE)
const$model = "const"

expon_sd = read_delim("simulations_models/exp.growth_S.2237_alpha.0.263700.txt", " ",escape_double = FALSE, col_names = c('pi','s','tajD'), 
                   comment = "#", trim_ws = TRUE)
expon_sd$model = "exp_growth_sd"


##### 
datos = bind_rows(sw,const,expon_sd)

datos = datos %>% pivot_longer(-model, "parameter","value") %>% filter(parameter != "s")
emp2 = data.frame("parameter" = c("pi","s","tajD"), "value" = as.numeric(emp))  %>% filter(parameter != "s")

mod_hist = ggplot(datos, aes(x = value))+
  geom_histogram(fill="white", color="black")+
  facet_grid(model~parameter, scales = "free_x")+
  geom_vline(data = emp2,
             aes(xintercept = value),color="red")

mod_bp = ggplot(datos, aes(x = model, y=value))+
  geom_boxplot(fill="white", color="black")+
  facet_grid(parameter~., scales = "free_y")+
  geom_hline(data = emp2,
             aes(yintercept = value),color="blue",linetype = "dashed")


mod_hist
mod_bp

###### Empirical Cumulative Distribution Function


my_ecdf = function(df,emp) {
  par(mfrow=c(1,2))
  p_pi = ecdf(df$pi)
  x = p_pi(emp$pi.m)
  
  plot(p_pi)
  abline(v = emp$pi.m,col='red')
  
  p_tajD = ecdf(df$tajD)
  y = p_tajD(emp$tajD.m)
  
  plot(p_tajD)
  abline(v = emp$tajD.m,col='red')
  
  res = c(x,y)
  names(res) = c("p_pi","p_tajD")
  return(res)
}

##### 

prob_models = data.frame(
           "constant" = my_ecdf(const,emp),
           "expon_sd" = my_ecdf(expon_sd,emp),
           "sweep" = my_ecdf(sw,emp)
           )

prob_models


#### boxplots for paper

p = filter(datos, parameter == "pi")

probs = prob_models[1,]
probs = round(probs, 4)

pi_all = ggplot(p, aes(x = model, y=value))+
  geom_boxplot(fill="white", color="black")+
  facet_grid(scales = "free_y")+
  geom_hline(data = filter(emp2, parameter == "pi"),
             aes(yintercept = value),color="blue",linetype = "dashed")+
  scale_x_discrete(labels=c(
    paste("Constant\nfrequency\n(p =", probs[1],")",sep = ""),
    paste("Exp.growth\n\U03B1 = 0.26\n(p =", probs[2],")",sep = ""),
    paste("Selective\nsweep\n(p =", probs[3],")",sep = "")))+
  labs(x="",y="\U03C0")+
  theme_bw(base_size = 20)+
  theme()


d = filter(datos, parameter == "tajD")
probs = prob_models[2,]
probs = round(probs, 4)

d_all = ggplot(d, aes(x = model, y=value))+
  geom_boxplot(fill="white", color="black")+
  facet_grid(scales = "free_y")+
  geom_hline(data = filter(emp2, parameter == "tajD"),
             aes(yintercept = value),color="blue",linetype = "dashed")+
  scale_x_discrete(labels=c(
    paste("Constant\nfrequency\n(p =", probs[1],")",sep = ""),
    paste("Exp.growth\n\U03B1 = 0.26\n(p =", probs[2],")",sep = ""),
    paste("Selective\nsweep\n(p =", probs[3],")",sep = "")))+
  labs(x="",y="Tajima's D")+
  theme_bw(base_size = 20)+
  theme()

bp_all = plot_grid(pi_all, d_all)
bp_all

