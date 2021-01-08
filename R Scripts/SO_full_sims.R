## SO_full_sims
# 
#========================================================	
# ---
### title: Spatial overlap networks, full simulations
# author: Marie Gilbertson
# date: "09/05/2020"
#---
# ##  Preamble	
# 
# What this code does:
# 1. Simulate outbreaks on SO-based networks (adapted from "FIV_full_sims.R")

##### Clear Environment #####
remove(list=ls())


#### load libraries ####
library(igraph)
library(dplyr) # sample_n function
library(spatstat)
library(DCG) # as.symmetricAdjacencyMatrix function
library(ergm)
library(intergraph)


#### load external functions ####
# function for simulating an overlap-based network
source("R Scripts/simulate_SOnet.R")
# function for simulating FeLV transmission on a network
source("R Scripts/trans_sim.R")
# function for processing outbreaks that included births (creates new columns for the "birthed" individuals)
source("R Scripts/post_process_outbreak_data.R")
# function for calculating proportions of population in each disease category when births are included in disease simulations
source("R Scripts/props affected_births included.R")
# function for extracting key results from each simulation
source("R Scripts/extract_results.R")

#### load LHS parameter sets ####
LHS.data <- get(load(file = "LHS parameter sets.Rdata"))

#### load degree distribution parameters ####
so.nbinom.params <- get(load("SO network sims_nbinom params.Rdata"))

params.to.run <- seq(1, 150)

#### MAIN SIMULATION LOOP ####
for(q in 1:length(params.to.run)){
  
  
  #### set parameter set number and pull parameters ####
  param.set.num <- params.to.run[q]
  print(paste("set ", param.set.num, sep = ""))
  
  
  # pull parameters
  params <- LHS.data[LHS.data$set.id==param.set.num,]
  
  
  ##### Set seed #####
  sim.seed <- 23498+param.set.num
  set.seed(sim.seed) # ensure unique but reproducible seeding
  
  
  #### set number of simulations ####
  nsims <- 50
  
  
  #### assign SO network simulation parameters ####
  pop.size <-  params$pop.size                  # size of simulated population 
  model.mu <- so.nbinom.params$nbinom.mean
  model.size <- so.nbinom.params$nbinom.size
  
  
  ### assign FeLV simulation parameters ####
  # transmission parameters
  c.rate <- params$c.rate                 # weekly probability of contact
  beta <- params$beta              # probability of transmission
  duration.yrs <- 2.5            # number of years to run simulation
  
  
  # weekly probability of death for progressives (1/progressive infection duration in weeks)
  # progressive infection duration is either: 1 = short; 2 = long
  if(params$progr.dur==1){
    progr.dur <- 1/18
  } else if(params$progr.dur==2){
    progr.dur <- 1/26
  } 
  
  # regressive disease parameters
  if(params$regr.dur.c==1){
    regr.dur.c <- 0.5
  }else if(params$regr.dur.c==2){
    regr.dur.c <- 1
  }
  
  # constant applied to beta to determine transmission rate of latent/regressive individuals 
  #(1 = 0, 2 = 0.1, 3 = 0.5, 4 = 1)
  if(params$c.beta==1){
    c.beta <- 0              
  }else if(params$c.beta==2){
    c.beta <- 0.1
  }else if(params$c.beta==3){
    c.beta <- 0.5
  }else if(params$c.beta==4){
    c.beta <- 1
  }
  
  # proportions of those infected that become progressive, latent/regressive, and abortive infections
  # proportions are: 1 = (25, 25, 50); can test other proportions as desired
  if(params$prop.outcomes==1){
    prop.outcomes <- c(1/4, 1/4, 1/2)
  } 
  
  # re-spawn ("birth") rate
  terr.repop <- params$terr.repop          # rate of territory repopulation after the death of the prior occupant (1/duration of wait in weeks)
  
  # vaccination paramters
  use.vax <- T               # True/False to use vaccination
  vax.rate <-  params$vax.rate             # vaccine rate (1/avg weeks to vaccinate) - this is actually the probability for binomial trials, but functioning as a rate here
  vax.eff <-  params$vax.eff            # vaccine efficacy (probability that an infectious contact is "interrupted" by vaccine)
  
  
  
  # prep for storing z-loop results
  full.sims.results <- data.frame(sim.num = numeric(),
                                  dur.time = numeric(),
                                  total.prog = numeric(),
                                  total.lr = numeric(),
                                  total.ab = numeric(),
                                  total.vax = numeric(),
                                  num.failed = numeric()
  )
  
  
  
  
  
  
  #### loop through simulations ####
  
  for(z in 1:nsims){
    print(z)
    sim.num <- z
    
    #### SIMULATE NETWORK ####
    net3 <- simulate_sonet(pop.size = pop.size,
                           model.mu = model.mu,
                           model.size = model.size
    )
    
    # net3
    # plot(net3)
    
    
    
    
    # net3 is already igraph object; rename for consistency
    g <- net3
    
    
    #### SIMULATE FELV ####
    # simulate FeLV outbreak
    m.list <- trans_sim(g = g, 
                         beta = beta, 
                         c.beta = c.beta,
                         c.rate = c.rate,                 
                         duration.yrs = duration.yrs,     
                         prop.outcomes = prop.outcomes,
                         progr.dur = progr.dur,
                         regr.dur.c = regr.dur.c,   
                         terr.repop = terr.repop,
                         use.vax = use.vax,               
                         vax.rate = vax.rate,
                         vax.eff = vax.eff
    )
    
    m <- m.list[[1]]
    
    m.new <- post_process_outbreak_data(m)
    
    p.results <- props_affected_births_included(m.new=m.new)
    
    # save results
    g.name <- paste("SO_Simulation_Results/g network_SO paramset ", param.set.num, "_sim ",sim.num, ".Rdata", sep = "")
    save(g, file = g.name)
    
    m.name <- paste("SO_Simulation_Results/m data_SO paramset ", param.set.num, "_sim ", sim.num, ".Rdata", sep = "")
    save(m, file = m.name)
    
    mnew.name <- paste("SO_Simulation_Results/mnew data_SO paramset ", param.set.num, "_sim ", sim.num, ".Rdata", sep = "")
    save(m.new, file = mnew.name)
    
    presults.name <- paste("SO_Simulation_Results/presults data_SO paramset ", param.set.num, "_sim ", sim.num, ".Rdata", sep = "")
    save(p.results, file = presults.name)
    
    
    
    ### plot network results for monitoring
    # note: just documents status of each "territory" at the end of simulation for original occupant; doesn't show status of new occupants
    # 0 means susceptible; 1 means infectious; 2 means regressive/latent infected; 3 means recovered/immune; 4 means dead from infection
    # 5 means recovered regressive; 6 means vaccinated
    
    # outcomes based on end of two years or outbreak end, whichever is earliest
    if(nrow(m.new)>=105){
      colors <- m.new[105,c(1:ncol(m))] # add one because row 1 = time 0 in m.new
    }else if(nrow(m.new)<105){
      colors <- m.new[nrow(m.new),c(1:ncol(m))]
    }
    
    colors[colors==0] <- "white"
    colors[colors==1] <- "red"
    colors[colors==2] <- "blue"
    colors[colors==3] <- "purple"
    colors[colors==4] <- "black"
    colors[colors==5] <- "lightblue"
    colors[colors==6] <- "yellow"
    
    g <- set_vertex_attr(g, "outcome", index = V(g), value = colors)
    
    lay <- layout.fruchterman.reingold(g)
    
    gplot.name <- paste("SO_Simulation_Results/Simulation_Figures/g plot_SO paramset ", param.set.num, "_sim ", sim.num, ".jpg", sep = "")
    jpeg(gplot.name)
    
    plot.igraph(g,
                vertex.color=V(g)$outcome,
                vertex.label=NA,
                vertex.size = 5,
                layout=lay,
                edge.color="dark grey",
                edge.width=1)
    
    dev.off()
    
    
    
    eplot.name <- paste("SO_Simulation_Results/Simulation_Figures/epi plot_SO paramset ", param.set.num, "_sim ", sim.num, ".jpg", sep = "")
    jpeg(eplot.name)
    
    plot(p.results$prop.i, col="red", type='l', main="Simulation of FeLV-like spread on a SO-based network",sub="",xlab="Time (weeks)",ylab="Population numbers",ylim=c(0,1),xlim=c(0, nrow(p.results)))   
    lines(p.results$prop.s, col="green")   #green is sum of #0  SUSCEPTIBLE
    lines(p.results$prop.lr, col="blue")     #blue is sum of #2  LATENT/REGRESSIVE
    lines(p.results$prop.r, col = "purple") #purple is sum of #3 ABORTIVE
    lines(p.results$prop.v, col = "yellow") # yellow is #6 VACCINATED

    dev.off()
    
    ##### Store primary outcomes of interest ######
    # store 1. duration of outbreak, 2. total progressive infections, 3. total latent/regressive infections, 4. total abortive infections, 5, total vaccinated, 6. number of "failed" epidemics (number that initiated with an isolate)
    
    temp.results <- extract_results(p.results.a = p.results,
                                    m.new.a = m.new,
                                    c.beta.a = c.beta,
                                    m.list.a = m.list,
                                    end.point = 2
    )
    
    
    full.sims.results[z,] <- temp.results 
    
  }
  
  full.sims.results$param.set <- paste("SO_set_", param.set.num, sep = "")
  
  # Save final results
  full.name <- paste("SO_Simulation_Results/full set results_SO paramset ", param.set.num, ".Rdata", sep = "")
  save(full.sims.results, file = full.name)
  
}

