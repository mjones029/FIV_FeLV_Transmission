# trans_sim
#
#========================================================	
# ---
### title: Function for simulating transmission of FeLV on networks
# author: Marie Gilbertson
# date: "03/31/2020"
#---
# ##  Preamble	
# 
# What this code does:
# 1. Simulates transmission of FeLV on previously generated networks (chain binomial process)




trans_sim <- function(g = g,                           # network on which to simulate (should be an igraph object)
                             beta = beta,                     # base transmission rate (applies to all progressive individuals)
                             c.beta = c.beta,                 # constant applied to beta to determine transmission rate of regressive individuals 
                             c.rate = c.rate,                 # weekly probability of contact
                             duration.yrs = duration.yrs,     # number of years to run simulation
                             prop.outcomes = prop.outcomes,   # proportion per each infected state
                             progr.dur = progr.dur,           # weekly probability of death for progressives (1/duration of infection)
                             regr.dur.c = regr.dur.c,         # modifier for regressive infection duration
                             terr.repop = terr.repop,         # rate of territory repopulation after the death of the prior occupant (1/duration of wait in weeks)
                             use.vax = use.vax,               # True/False to use vaccination
                             vax.rate = vax.rate,             # vaccine rate (1/avg weeks to vaccinate) - this is actually the probability for binomial trials, but functioning as a rate here
                             vax.eff = vax.eff                # vaccine efficacy (0 to 1) - probability in binomial trials
) {
  
  # assign pop.size value 
  pop.size <- length(V(g))
  
  
  
  
  
  ############# Network model ########################
  
  ER<-get.adjacency(g) #make an adjacency network; convert graph to an adjacency matrix or an edge list
  
  
  netwb <- ER
  
  
  #### assign FeLV model parameters ####
  beta<-beta 
  sigma <- c.beta*beta # transmission rate for regressively infected individuals
  regr.dur <- regr.dur.c*progr.dur # multiply "c" constant by the progressive duration distribution to get value for regressives
  duration <- 52*duration.yrs # 52 * number of years to run simulation
  
  
  use.vax <- use.vax
  if(use.vax==T){
    vax.start <- F
  }
  vax.rate <- vax.rate
  vax.eff <- vax.eff
  # empty vector to record the number of successful vaccination events
  total.vax <- c()
  
  #### run the disease simulation ####
  time<-1
  gp <- ncol(netwb)    # number of columns in network    
  
  # randomly select individual initiating infection in the population
  # keep record of how many iterations before a non-isolate is selected
  
  initiate <- c()
  repeat {
    
    stat<-sample(c(numeric(gp-1),1));   #choose first node to get infected. This shuffles the order of one 1 and the rest zeros (to make bigger network, need to add more 0's)
    #sample function can be used to return a random permutation of a vector
    #numeric creates a real vector of the specified length. The elements of the vector are all 0.
    #c combines values into a vector or list
    
    initiate <- c(initiate, which(stat==1))
    
    # exit loop when infection initiates in a non-isolate
    if (any(netwb[which(stat==1),]==1) & any(netwb[,which(stat==1)]==1)) break
  }
  
  # "initiate.dur" gives the number of iterations before a non-isolate was selected (+1)
  initiate.dur <- length(initiate)
  
  
  m <- matrix(0,duration+1,gp)  # make empty matrix for number of time steps long by number of individuals wide (this will be the progression of infection per time step per indiv.) 52 plus 1 more for time=0
  m[time,] <- stat
  
  # create matrix to accommodate duration of infection assignments
  n <- matrix(0,duration*100,gp) # need this to be extra long to accommodate duration of infection past end of monitoring
  
  # assign duration of infection for first infected individual
  recovery.all <- c()
  repeat {
    # repeatedly draw until recovery = 1, recording number of draws until that point
    recovery <- rbinom(1,1,progr.dur)
    recovery.all <- c(recovery.all, recovery)
    
    # exit loop when death = 1
    if (recovery > 0) break
  }
  dur.1 <- length(recovery.all)
  n[c(1:dur.1),which(stat==1)] <- 1
  
  # (optional)
  write(c(0,stat),file="test_sim.txt",ncolumns=(gp+1),append=F) #writes initial conditions to file
  
  # now that a non-isolate is infected, proceed with simulation timesteps
  while ((any(stat==1)|any(stat==2))&time<=duration) {

    #0 means susceptible; 1 means infectious; 2 means regressively/latently infected; 3 means recovered/immune; 4 means dead from infection; 5 means "recovered" latent/regressive individuals; 6 means vaccinated
    
    # for each (weekly) time step
    statc <- stat  #assign stat to copy of stat so can change 'statc' while still working off of 'stat'
    
    #### vaccination ####
    if(use.vax==T){ # if using vaccination...
      if(vax.start==F){ # check if vaccination has started, and if not...
        if(time>52){ # check if simulation is through the first year
          vax.start <- TRUE # if so, turn vaccination "on"
        }
      }
      
      if(vax.start==T){ # if vaccination has started...
        
        if(rbinom(1,1,vax.rate)==1){ # decide if anyone gets vaccinated this time step, based on vax rate...
          
          vax.ind <- sample(seq(1:length(stat)),1) # and if so, randomly select the identity of individual to vaccinate
          
          # disease/compartment status will only change if the selected individual was susceptible
          statc[vax.ind] <- ifelse(stat[vax.ind]==0, 6, stat[vax.ind]) # "6" denotes vaccinated individuals
          
          # add to vector to store total number EVER vaccinated (to account for some vax individuals potentially becoming infected)
          total.vax <- c(total.vax, ifelse(stat[vax.ind]==0, 1, 0))
        }
      }
    }
    
    #### infection of susceptibles ####
    for (i in (which(stat==0))){     #for every susceptible individual in original stat...
      if(rbinom(1,1,c.rate)==1){    #does that individual have a contact this time step? If so...
        for (j in (1:gp)[-i]){          #for each j in 1 to number groups, with exclusion of itself (because that is always 0)
          
          if ((netwb[i,j]>=1) & (((stat[j]==1)&(rbinom(1,1,beta)==1)) | ((stat[j]==2)&(rbinom(1,1,sigma)==1))) ) #is there an edge to an infected individual and does a transmission event take place?  (By looking down column of that indiv in netw for a 1)  
          {
            
            # does the infected individual become progressively, latently, or regressively infected?
            state <- sample(c(1:3), size = 1, prob = prop.outcomes) # 1 = progressive, 2 = regressive, 3 = abortive
            statc[i] <- state
            
            if(state==1){ 
              # assign duration of infection
              recovery.all <- c()
              repeat {
                # repeatedly draw until recovery = 1, recording number of draws until that point
                recovery <- rbinom(1,1,progr.dur)
                recovery.all <- c(recovery.all, recovery)
                
                # exit loop when death = 1
                if (recovery > 0) break
              }
              dur.temp <- length(recovery.all)
              n[c((time+1):(time+dur.temp)),i] <- 1
            }
            if(state==2){
              # assign duration of infection
              recovery.all <- c()
              repeat {
                # repeatedly draw until recovery = 1, recording number of draws until that point
                recovery <- rbinom(1,1,regr.dur)
                recovery.all <- c(recovery.all, recovery)
                
                # exit loop when death = 1
                if (recovery > 0) break
              }
              dur.temp <- length(recovery.all)
              n[c((time+1):(time+dur.temp)),i] <- 2
            }
            if(state==3){ # assign "duration" for abortive infections (for record keeing purposes)
              n[c((time+1):(duration+1)),i] <- 3 # this bit only matters for record keeping
            }
            break # assign outcome of infection to the copy and get out of loop
          }
        }
      }
    }
    
    #### infection of vaccinated individuals ####
    # loop through all vaccinated individuals (same basic loop as for susceptibles, but now accounting for vaccine efficacy)
    if(use.vax==T){ # if using vaccination...
      if(vax.start==T){ # & if vaccination is turned on (have to have separated out for instances where not using vaccination)
        
        for (i in (which(stat==6))){     #for every vaccinated individual in original stat...
          if(rbinom(1,1,c.rate)==1){    #does that individual have a contact this time step? If so...
            for (j in (1:gp)[-i]){          #for each j in 1 to number groups, with exclusion of itself (because that is always 0)
              
              if ((netwb[i,j]>=1) & (((stat[j]==1)&(rbinom(1,1,beta)==1)&(rbinom(1,1,vax.eff)==0)) | ((stat[j]==2)&(rbinom(1,1,sigma)==1)&(rbinom(1,1,vax.eff)==0))) ) #is there an edge to an infected individual and does a transmission event take place and does the vaccine fail to prevent infection? 
              {
                
                # does the infected individual become progressively, latently, or regressively infected?
                state <- sample(c(1:3), size = 1, prob = prop.outcomes) # 1 = progressive, 2 = regressive, 3 = abortive
                statc[i] <- state
                
                if(state==1){ 
                  # assign duration of infection
                  recovery.all <- c()
                  repeat {
                    # repeatedly draw until recovery = 1, recording number of draws until that point
                    recovery <- rbinom(1,1,progr.dur)
                    recovery.all <- c(recovery.all, recovery)
                    
                    # exit loop when death = 1
                    if (recovery > 0) break
                  }
                  dur.temp <- length(recovery.all)
                  n[c((time+1):(time+dur.temp)),i] <- 1
                }
                if(state==2){
                  # assign duration of infection
                  recovery.all <- c()
                  repeat {
                    # repeatedly draw until recovery = 1, recording number of draws until that point
                    recovery <- rbinom(1,1,regr.dur)
                    recovery.all <- c(recovery.all, recovery)
                    
                    # exit loop when death = 1
                    if (recovery > 0) break
                  }
                  dur.temp <- length(recovery.all)
                  n[c((time+1):(time+dur.temp)),i] <- 2
                }
                if(state==3){ # assign "duration" for regressives (for record keeing purposes)
                  n[c((time+1):(duration+1)),i] <- 3 # this bit only matters for record keeping
                }
                break # assign outcome of infection to the copy and get out of loop
              }
            }
          }
        }
      }
    }
    
    #### death of progressively infected individuals ####
    for (i in (which(stat==1)))     #for every progressively infectious individual in original stat
    {
      if (n[time+1,i]==0) # determine which individuals die from infection
      {
        statc[i] <- 4
        # Assign a "duration of death" to account for new "births" occupying empty territories
        
        # first, draw binomial trials to determine time until territory is repopulated
        # if/else loop to accommodate no birth/repopulation rate
        if(terr.repop>0){ # if territories are repopulable, assign duration of time until repopulation
          death.all <- c()
          
          repeat {
            # repeatedly draw until death = 1, recording number of draws until that point
            death <- rbinom(1,1,terr.repop)
            death.all <- c(death.all, death)
            
            # exit loop when death = 1
            if (death > 0) break
          }
          death.dur <- length(death.all) # assign time until repopulation of territory 
        }else{ # if territories are not repopulable, just assign out death duration through end of simulation
          death.dur <- duration+1 # will put 4s out too far, but does the job
        }
        
        # then assign that death duration in the n matrix
        n[c((time+1):(time+death.dur)),i] <- 4 
        
      }
    }
    
    #### recovery of regressively infected individuals ####
    for (i in (which(stat==2)))     #for every latently/regressively infectious individual in original stat
    {
      if (n[time+1,i]==0) # determine which individuals recover from infection
      {
        statc[i] <- 5
        n[c((time+1):(duration+1)),i] <- 5
      }
    }
    
    #### respawning ####
    # determine which "territories" are re-occupied by a new, susceptible individual ("births") after prior individual's death
    for (i in (which(stat==4)))     #for every dead individual in original stat
    {
      if (n[time+1,i]==0) # determine which territories are re-occupied
      {
        statc[i] <- 0
        n[c((time+1):(duration+1)),i] <- 0
        # THIS FORMUALTION MEANS THAT A SINGLE COLUMN COULD REPRESENT MORE THAN ONE INDIVIDUAL
        # MUST BE CAPTURED IN "POST-PROCESSING"
      }
    }
    
    # update recording
    stat <- statc
    m[time+1,] <- stat
    write(c(time,m[time+1,]),file="test_sim.txt",ncolumns=(gp+1),append=TRUE)
    time<-(time+1)
  }
  
  # assemble results
  m <- m[rowSums(m)>0,]
  m.list <- list(m, initiate.dur, total.vax)
  names(m.list) <- c("m", "initiate.dur", "total.vax")
  return(m.list)
  
}