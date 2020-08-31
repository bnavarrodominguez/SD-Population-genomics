#setwd("~/Documentos/investigacion/ur_larracuente_presgraves/SD_SNPs/ms_abc_simulations/ms_and_abc_R/")

####### Instantaneous population size changes
#### From ms manual (https://uchicago.app.box.com/s/l3e5uf13tikfjm7e1il1eujitlsjdx13/file/245765534272?sb=/details) 

#### Define fixed parameters

#### D. melanogaster effective population size
#ne = 10e6

### D. melanogaster Ne in Zambia (Kapopoulou 2018)
ne = 3160475


### Frequency of In(2R)Mal - 1.47%
f = 0.0147

#####Mutation rate perB site per generation (Laurent et al 2011, used in Kapopoulou 2018)
mutrate = 1.3*10^-9

##### Sequence lenght
l = 9500000

#### Theta = 4*Ne*mut.rate*f*l (modified from Corbett-Detig & Hartl 2012) 
theta = 4*ne*mutrate*f*l
theta


#### N0 = Ne * frequency of In(2R)Mal
n0 = ne*f
n0

### proportion in reduction of population size in the bottleneck
### since n1 = 1, x1 must be 1/n0
x1 = 1/n0
x1
x1*n0 #### check it equals 1


ngen_year = 36.5