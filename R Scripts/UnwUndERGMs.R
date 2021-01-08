## Unweighted/Undirected ERGMs
# 
#========================================================	
# ---
### title: Phyloscanner Unweighted/Undirected ERGM workflow
# author: Marie Gilbertson
# date: "06/19/2020"
#---
###  Preamble	
# 
# What this code does:
# 1. Phyloscanner transmission network ERGM analysis.
# 2. Specific to unweighted/undirected transmission network.
# Note: see also supplementary materials for Silk et al 2019 (dog dominance networks) for 
# more extensive ERGM fitting tutorials/examples.

##### Clear Environment
remove(list=ls())


#### load initial libraries ####
library(igraph)
library(latticeExtra) 
library(fields)
library(plyr)

##### set seed ######

ergm.seed <- 3128


#### PREPARE DATA ####

##### Read in attribute data ##### 
# cats.fiv includes:
# 1. sex
# 2. age
# 3. freeway clusters (I-75)

cats.fiv <- get(load("Attribute Data/FL FIV covariates.Rdata"))

# update single "kitten" to be in same category as other non-adults (i.e. subadults)
cats.fiv$Age.categorical[cats.fiv$Age.categorical=="Kitten"] <- "Yearling"


# pairwise genetic distance
relatedness <- get(load("Attribute Data/FL FIV_panther relatedness.Rdata"))
# can't have NA's, so convert diagonal to zero
relatedness[is.na(relatedness)] <- 0
kin <- relatedness


# pairwise geographic distances (log adjusted; Inf replaced by 0)
pair.dist <- get(load("Attribute Data/FL FIV pairwise dists_logkm.Rdata"))


# spatial overlap matrix
so.matrix <- get(load("Attribute Data/FL FIV pairwise overlap_UDOI.Rdata"))



##### Read in FIV transmission network (outcome); Phyloscanner igraph object #####

ph.g <- get(load("Phyloscanner_networks/FL FIV transmission network.Rdata"))

# convert to adjacency matrix 
tt.matrix <- as_adj(ph.g)
tt.matrix <- as.matrix(tt.matrix)



# verify that ordering of individuals is consistent with attribute data objects
identical(cats.fiv$SVRG.ID, rownames(tt.matrix))

# check if cats.fiv ordering is the same as edge attributes
identical(cats.fiv$SVRG.ID, rownames(kin))
identical(cats.fiv$SVRG.ID, rownames(pair.dist))
identical(cats.fiv$SVRG.ID, rownames(so.matrix))


#### convert tt.matrix for ERGM analysis ####
# now detach igraph and load ergm libraries
detach("package:igraph", unload=TRUE)
#load sna for network regressions
library(sna)
library(tnet)
library(ergm)

# convert to network format
tt.net <- as.network.matrix(tt.matrix, directed = F) # "directed = F" makes this an undirected network


# view network
plot.network(tt.net,
             displayisolates = T)

# view network density
network.density(tt.net) 

#### Now add attributes to transmission network ####
# cats.fiv includes:
# 1. sex
# 2. age (continuous and categorical)
# 3. freeway clusters (I-75)
# Other matrices:
# 4. relatedness
# 5. pairwise distances
# 6. spatial overlap

# add sex as a numeric attribute 
network::set.vertex.attribute(tt.net,"sex",as.integer(cats.fiv$bin.sex))

# add age as a categorical variable
network::set.vertex.attribute(tt.net,"age.cat",as.character(cats.fiv$Age.categorical))

# add freeway cluster category (N = 0, S = 1)
network::set.vertex.attribute(tt.net,"freeway",as.integer(cats.fiv$bin.freeway))

# add distance to urban areas as a continuous variable
network::set.vertex.attribute(tt.net, "near.dist", as.numeric(cats.fiv$near.dist))

# other covariates are edge-based or structural and will be handled within models


##### ERGM ANALYSIS #####

#### start by evaluating basic, intercept-style edge metric ####
mod.edges <- ergm(tt.net~edges, control = control.ergm(seed = ergm.seed)) 


#### now perform forward selection of structural covariates of interest ####

# Structural covariates of interest:
# 1. gwesp
# 2. altkstar
# 3. twopaths
# 4. isolates

## single variable models:

mod.S.1.1 <- ergm(tt.net~edges + gwesp(0.5, fixed=T), control = control.ergm(seed = ergm.seed))
mod.S.1.2 <- ergm(tt.net~edges + altkstar(0.5, fixed = T), control = control.ergm(seed = ergm.seed))
mod.S.1.3 <- ergm(tt.net~edges + twopath, control = control.ergm(seed = ergm.seed))
mod.S.1.4 <- ergm(tt.net~edges + isolates, control = control.ergm(seed = ergm.seed))
# twopath and isolates fail, so don't include in further model fitting

aics.S.1 <- AIC(mod.edges, 
                mod.S.1.1,
                mod.S.1.2
                #mod.S.1.3,
                # mod.S.1.4
)
aics.S.1[order(aics.S.1$AIC),]

# mod.S.1.1 (edges + gwesp(0.5, fixed = T)) is best model, so carry forward

# however, test out different parameterizations of gwesp and altkstar before moving on
mod.S.1.1_3 <- ergm(tt.net~edges + gwesp(0.3, fixed=T), control = control.ergm(seed = ergm.seed))
mod.S.1.1_4 <- ergm(tt.net~edges + gwesp(0.4, fixed=T), control = control.ergm(seed = ergm.seed))
mod.S.1.1_5 <- ergm(tt.net~edges + gwesp(0.5, fixed=T), control = control.ergm(seed = ergm.seed))
mod.S.1.1_6 <- ergm(tt.net~edges + gwesp(0.6, fixed=T), control = control.ergm(seed = ergm.seed))
mod.S.1.1_7 <- ergm(tt.net~edges + gwesp(0.7, fixed=T), control = control.ergm(seed = ergm.seed))
mod.S.1.1_8 <- ergm(tt.net~edges + gwesp(0.8, fixed=T), control = control.ergm(seed = ergm.seed))
mod.S.1.1_9 <- ergm(tt.net~edges + gwesp(0.9, fixed=T), control = control.ergm(seed = ergm.seed))
# 0.9 fit fails

aics.S.1_b <- AIC(mod.S.1.1_3,
                  mod.S.1.1_4,
                  mod.S.1.1_5,
                  mod.S.1.1_6,
                  mod.S.1.1_7,
                  mod.S.1.1_8
                  # mod.S.1.1_9
)
aics.S.1_b[order(aics.S.1_b$AIC),]

# All weights for GWESP are statistically significant. Typically, more
# extreme weights yield poor MCMC diagnostics; for now, keep gwesp(0.5), but test sensitivity
# in final model.

mod.S.1.2_3 <- ergm(tt.net~edges + altkstar(0.3, fixed=T), control = control.ergm(seed = ergm.seed))
mod.S.1.2_4 <- ergm(tt.net~edges + altkstar(0.4, fixed=T), control = control.ergm(seed = ergm.seed))
mod.S.1.2_5 <- ergm(tt.net~edges + altkstar(0.5, fixed=T), control = control.ergm(seed = ergm.seed))
mod.S.1.2_6 <- ergm(tt.net~edges + altkstar(0.6, fixed=T), control = control.ergm(seed = ergm.seed))
mod.S.1.2_7 <- ergm(tt.net~edges + altkstar(0.7, fixed=T), control = control.ergm(seed = ergm.seed))
mod.S.1.2_8 <- ergm(tt.net~edges + altkstar(0.8, fixed=T), control = control.ergm(seed = ergm.seed))
mod.S.1.2_9 <- ergm(tt.net~edges + altkstar(0.9, fixed=T), control = control.ergm(seed = ergm.seed))

# weights 0.3, 0.4, and 0.9 failed

aics.S.1_b <- AIC(#mod.S.1.2_3,
                  #mod.S.1.2_4,
                  mod.S.1.2_5,
                  mod.S.1.2_6,
                  mod.S.1.2_7,
                  mod.S.1.2_8
                  #mod.S.1.2_9
)
aics.S.1_b[order(aics.S.1_b$AIC),]
# all are virtually equivalent.



## Two structural variable model:
mod.S.2.1_a <- ergm(tt.net~edges + gwesp(0.5, fixed=T) + altkstar(0.5, fixed = T),  control = control.ergm(seed = ergm.seed))
mod.S.2.1_b <- ergm(tt.net~edges + gwesp(0.5, fixed=T) + altkstar(0.6, fixed = T),  control = control.ergm(seed = ergm.seed))
# 1_b failed

aics.S.2 <- AIC(mod.edges, 
                mod.S.1.1,
                mod.S.2.1_a
                # mod.S.2.1_b
)
aics.S.2[order(aics.S.2$AIC),]

# mod.S.1.1 is best model, but 1_a is only about ∆2 AIC values difference
# check the effect of altering weight of gwesp argument:

mod.S.2.1_a1 <- ergm(tt.net~edges + gwesp(0.5, fixed=T) + altkstar(0.5, fixed = T),  control = control.ergm(seed = ergm.seed))
mod.S.2.1_a2 <- ergm(tt.net~edges + gwesp(0.6, fixed=T) + altkstar(0.5, fixed = T),  control = control.ergm(seed = ergm.seed))
mod.S.2.1_a3 <- ergm(tt.net~edges + gwesp(0.7, fixed=T) + altkstar(0.5, fixed = T),  control = control.ergm(seed = ergm.seed))

mod.S.2.1_b1 <- ergm(tt.net~edges + gwesp(0.5, fixed=T) + altkstar(0.6, fixed = T),  control = control.ergm(seed = ergm.seed))
mod.S.2.1_b2 <- ergm(tt.net~edges + gwesp(0.6, fixed=T) + altkstar(0.6, fixed = T),  control = control.ergm(seed = ergm.seed))
mod.S.2.1_b3 <- ergm(tt.net~edges + gwesp(0.7, fixed=T) + altkstar(0.6, fixed = T),  control = control.ergm(seed = ergm.seed))

aics.S.2_b <- AIC(mod.edges, # edges
                  mod.S.1.1,
                  mod.S.2.1_a1,
                  mod.S.2.1_a2,
                  mod.S.2.1_a3
                  # mod.S.2.1_b1,
                  # mod.S.2.1_b2,
                  # mod.S.2.1_b3
)
aics.S.2_b[order(aics.S.2_b$AIC),]

summary(mod.S.2.1_a3)

# mod.S.2.1_a3 has best AIC, yet altkstar is not statistically significant. Also suspect that 
# gwesp will become harder to fit if keep increasing weight. Let's examine GOF to see if 
# addition of altkstar variable is justified by GOF.

m.gof_1 <-gof(mod.S.1.1 ~ degree  + distance  + triadcensus, control = control.gof.ergm(seed = ergm.seed))
m.gof_2 <-gof(mod.S.2.1_a3 ~ degree  + distance  + triadcensus, control = control.gof.ergm(seed = ergm.seed))

plot(m.gof_1)
plot(m.gof_2)

# Marginally better GOF for mod.S.2.1_a3, so carry that model forward.


###### Now add node and edge covariates, controlling for structure ######
# node covariates
# 1. sex
# 2. age (continuous and categorical)
# 3. freeway clusters (I-75)
# 4. distance to nearest urban area (in kilometers)
# edge covariates
# 5. relatedness
# 6. pairwise distances
# 7. spatial overlap



## single variable models:
mod.1.1 <- ergm(tt.net~edges + gwesp(0.7, fixed=T) + altkstar(0.5, fixed = T) + nodefactor("sex"),  control = control.ergm(seed = ergm.seed))
mod.1.2 <- ergm(tt.net~edges + gwesp(0.7, fixed=T) + altkstar(0.5, fixed = T) + nodemix("sex"),  control = control.ergm(seed = ergm.seed))
mod.1.3 <- ergm(tt.net~edges + gwesp(0.7, fixed=T) + altkstar(0.5, fixed = T) + nodefactor("age.cat", base = 2),  control = control.ergm(seed = ergm.seed))
mod.1.4 <- ergm(tt.net~edges + gwesp(0.7, fixed=T) + altkstar(0.5, fixed = T) + nodemix("age.cat"),  control = control.ergm(seed = ergm.seed))
mod.1.5 <- ergm(tt.net~edges + gwesp(0.7, fixed=T) + altkstar(0.5, fixed = T) + nodematch("freeway"),  control = control.ergm(seed = ergm.seed))
mod.1.6 <- ergm(tt.net~edges + gwesp(0.7, fixed=T) + altkstar(0.5, fixed = T) + nodecov("near.dist"),  control = control.ergm(seed = ergm.seed))
mod.1.7 <- ergm(tt.net~edges + gwesp(0.7, fixed=T) + altkstar(0.5, fixed = T) + edgecov(kin),  control = control.ergm(seed = ergm.seed))
mod.1.8 <- ergm(tt.net~edges + gwesp(0.7, fixed=T) + altkstar(0.5, fixed = T) + edgecov(pair.dist),  control = control.ergm(seed = ergm.seed))
mod.1.9 <- ergm(tt.net~edges + gwesp(0.7, fixed=T) + altkstar(0.5, fixed = T) + edgecov(so.matrix),  control = control.ergm(seed = ergm.seed))


aics.1 <- AIC(mod.edges, 
              mod.S.2.1_a3,
              mod.1.1,
              mod.1.2,
              mod.1.3,
              mod.1.4,
              mod.1.5,
              mod.1.6,
              mod.1.7,
              mod.1.8,
              mod.1.9
)
aics.1[order(aics.1$AIC),]
# mod.1.3 is best (greatest decrease in AIC from mod.S.2.1_a3), so carry this model forward


## add second dyad-independent variable
mod.2.1 <- ergm(tt.net~edges + gwesp(0.7, fixed=T) + altkstar(0.5, fixed = T) + nodefactor("age.cat", base = 2) + nodefactor("sex"),  control = control.ergm(seed = ergm.seed))
mod.2.2 <- ergm(tt.net~edges + gwesp(0.7, fixed=T) + altkstar(0.5, fixed = T) + nodefactor("age.cat", base = 2) + nodemix("sex"),  control = control.ergm(seed = ergm.seed))
mod.2.3 <- ergm(tt.net~edges + gwesp(0.7, fixed=T) + altkstar(0.5, fixed = T) + nodefactor("age.cat", base = 2) + nodemix("age.cat"),  control = control.ergm(seed = ergm.seed))
mod.2.4 <- ergm(tt.net~edges + gwesp(0.7, fixed=T) + altkstar(0.5, fixed = T) + nodefactor("age.cat", base = 2) + nodematch("freeway"),  control = control.ergm(seed = ergm.seed))
mod.2.5 <- ergm(tt.net~edges + gwesp(0.7, fixed=T) + altkstar(0.5, fixed = T) + nodefactor("age.cat", base = 2) + nodecov("near.dist"),  control = control.ergm(seed = ergm.seed))
mod.2.6 <- ergm(tt.net~edges + gwesp(0.7, fixed=T) + altkstar(0.5, fixed = T) + nodefactor("age.cat", base = 2) + edgecov(kin),  control = control.ergm(seed = ergm.seed))
mod.2.7 <- ergm(tt.net~edges + gwesp(0.7, fixed=T) + altkstar(0.5, fixed = T) + nodefactor("age.cat", base = 2) + edgecov(pair.dist),  control = control.ergm(seed = ergm.seed))
mod.2.8 <- ergm(tt.net~edges + gwesp(0.7, fixed=T) + altkstar(0.5, fixed = T) + nodefactor("age.cat", base = 2) + edgecov(so.matrix),  control = control.ergm(seed = ergm.seed))



aics.2 <- AIC(mod.edges,
              mod.S.2.1_a3,
              mod.1.3,
              mod.2.1, 
              mod.2.2,
              mod.2.3,
              mod.2.4,
              mod.2.5,
              mod.2.6,
              mod.2.7,
              mod.2.8
)
aics.2[order(aics.2$AIC),]
# mod.2.7 is best (greatest decrease in AIC from mod.1.3 and mod.S.2.1_a3), so carry this model forward




## add a third dyad-independent variable
mod.3.1 <- ergm(tt.net~edges + gwesp(0.7, fixed=T) + altkstar(0.5, fixed = T) + nodefactor("age.cat", base = 2) + edgecov(pair.dist) + nodefactor("sex"),  control = control.ergm(seed = ergm.seed))
mod.3.2 <- ergm(tt.net~edges + gwesp(0.7, fixed=T) + altkstar(0.5, fixed = T) + nodefactor("age.cat", base = 2) + edgecov(pair.dist) + nodemix("sex"),  control = control.ergm(seed = ergm.seed))
mod.3.3 <- ergm(tt.net~edges + gwesp(0.7, fixed=T) + altkstar(0.5, fixed = T) + nodefactor("age.cat", base = 2) + edgecov(pair.dist) +  nodemix("age.cat"),  control = control.ergm(seed = ergm.seed))
mod.3.4 <- ergm(tt.net~edges + gwesp(0.7, fixed=T) + altkstar(0.5, fixed = T) + nodefactor("age.cat", base = 2) + edgecov(pair.dist) + nodematch("freeway"),  control = control.ergm(seed = ergm.seed))
mod.3.5 <- ergm(tt.net~edges + gwesp(0.7, fixed=T) + altkstar(0.5, fixed = T) + nodefactor("age.cat", base = 2) + edgecov(pair.dist) + nodecov("near.dist"),  control = control.ergm(seed = ergm.seed))
mod.3.6 <- ergm(tt.net~edges + gwesp(0.7, fixed=T) + altkstar(0.5, fixed = T) + nodefactor("age.cat", base = 2) + edgecov(pair.dist) + edgecov(kin),  control = control.ergm(seed = ergm.seed))
mod.3.7 <- ergm(tt.net~edges + gwesp(0.7, fixed=T) + altkstar(0.5, fixed = T) + nodefactor("age.cat", base = 2) + edgecov(pair.dist) + edgecov(so.matrix),  control = control.ergm(seed = ergm.seed))


aics.3 <- AIC(mod.edges, 
              mod.S.2.1_a3,
              mod.1.3,
              mod.2.7,
              mod.3.1,
              mod.3.2,
              mod.3.3,
              mod.3.4,
              mod.3.5,
              mod.3.6,
              mod.3.7
)
aics.3[order(aics.3$AIC),]
# mod.2.7 is still best model, so stop forward selection


## do another evaluation of sensitivity of this final model to variation in gwesp and altkstar weighting
mod.2.7_a1 <- ergm(tt.net~edges + gwesp(0.5, fixed=T) + altkstar(0.5, fixed = T) + nodefactor("age.cat", base = 2) + edgecov(pair.dist),  control = control.ergm(seed = ergm.seed))
mod.2.7_a2 <- ergm(tt.net~edges + gwesp(0.6, fixed=T) + altkstar(0.5, fixed = T) + nodefactor("age.cat", base = 2) + edgecov(pair.dist),  control = control.ergm(seed = ergm.seed))
mod.2.7_a3 <- ergm(tt.net~edges + gwesp(0.7, fixed=T) + altkstar(0.5, fixed = T) + nodefactor("age.cat", base = 2) + edgecov(pair.dist),  control = control.ergm(seed = ergm.seed))
mod.2.7_a4 <- ergm(tt.net~edges + gwesp(0.8, fixed=T) + altkstar(0.5, fixed = T) + nodefactor("age.cat", base = 2) + edgecov(pair.dist),  control = control.ergm(seed = ergm.seed))

mod.2.7_b1 <- ergm(tt.net~edges + gwesp(0.5, fixed=T) + altkstar(0.6, fixed = T) + nodefactor("age.cat", base = 2) + edgecov(pair.dist),  control = control.ergm(seed = ergm.seed))
mod.2.7_b2 <- ergm(tt.net~edges + gwesp(0.6, fixed=T) + altkstar(0.6, fixed = T) + nodefactor("age.cat", base = 2) + edgecov(pair.dist),  control = control.ergm(seed = ergm.seed))
mod.2.7_b3 <- ergm(tt.net~edges + gwesp(0.7, fixed=T) + altkstar(0.6, fixed = T) + nodefactor("age.cat", base = 2) + edgecov(pair.dist),  control = control.ergm(seed = ergm.seed))
mod.2.7_b4 <- ergm(tt.net~edges + gwesp(0.8, fixed=T) + altkstar(0.6, fixed = T) + nodefactor("age.cat", base = 2) + edgecov(pair.dist),  control = control.ergm(seed = ergm.seed))

# b2 and b3 failed

aics.sens <- AIC(mod.edges, # edges
                 mod.2.7,
                 mod.2.7_a1,
                 mod.2.7_a2,
                 mod.2.7_a3,
                 mod.2.7_a4,
                 mod.2.7_b1,
                 # mod.2.7_b2,
                 # mod.2.7_b3,
                 mod.2.7_b4
)
aics.sens[order(aics.sens$AIC),]
# models a2, a3, a4, and b1 are virtually equivalent with original (mod.2.7 = mod.2.7_a3)
# therefore justified in carrying forward mod.2.7 as best final model


#### EVALUATE FIT AND MCMC DIAGNOSTICS ####

### assign best model
combo.mod.1 <- mod.2.7


### evaluate diagnostics
mcmc.diagnostics(combo.mod.1)
# looks ok, but verify with longer chain


combo.mod.2 <- ergm(tt.net~edges + gwesp(0.7, fixed=T) + altkstar(0.5, fixed = T) + nodefactor("age.cat", base = 2) + edgecov(pair.dist),  
                    control = control.ergm(seed = ergm.seed, MCMC.samplesize = 4096, MCMC.interval = 8192), verbose = T)

mcmc.diagnostics(combo.mod.2)
# looks ok

AIC(combo.mod.1, combo.mod.2)
summary(combo.mod.1)
summary(combo.mod.2)
# no significant differences in results with extension of chain



### goodness of fit ###

# goodness of fit for degree distribution, geodesic distance, and triad census (not including ESP as this was in the original model)
m.gof <-gof(combo.mod.2 ~ degree  + distance  + triadcensus, control = control.gof.ergm(seed = ergm.seed))

par(mfrow = c(2,2))
plot(m.gof)
par(mfrow = c(1,1))


# check structure of simulated networks
combo.mod <- combo.mod.2

sim1 <- simulate(combo.mod, nsim = 1)
plot.network(sim1,
             displayisolates = T)



summary(combo.mod)

# save best model
# save(combo.mod, file = "Phyloscanner_networks/FL FIV_bestERGM.Rdata")
