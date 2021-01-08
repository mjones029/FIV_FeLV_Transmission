## FeLV_homogmix.R
# 
#========================================================	
# ---
### title: FeLV ODE model
# author: Marie Gilbertson
# date: "05/06/2020"
#---
# ##  Preamble	
# 
# What this code does:
# 1. Homogeneous mixing compartmental model (Gillespie algorithm) of FeLV transmission in Florida panthers
# Intended to compare to model predictions of FeLV transmission using FIV networks 
# (e.g. does FIV-based approach perform better than compartmental model in capturing empirically observed dynamics)

##### Clear Environment #####
remove(list=ls())


#### load libraries ####
library("plyr")     ## for rdply()
library("reshape2") ## for melt()
library("emdbook")  ## for lambertW()
library("ggplot2"); theme_set(theme_bw())


#### load LHS parameter sets ####
LHS.data <- get(load(file = "LHS parameter sets.Rdata"))
reps <- 50 # number of simulations per parameter set

params.to.run <- seq(1, 150)


##### Set seed #####
set.seed(67459)


#### LOAD FUNCTIONS ####

#### this version of Gillespie algorithm adapted from Ben Bolker, https://rpubs.com/bbolker/SIRgillespie ####

### Function for computing the event rates
ratefun <- function(y,p,t) {
  if(t>52){ # if after first year, add vaccination
    with(as.list(c(y,p)),{
      c(inf=((c.rate)*((beta*net.den*S*Ip) + (c.beta*beta*net.den*S*Ir))),  ## progressive infection + regressive infection
        infv=((c.rate)*((beta*net.den*V*(1-ve)*Ip) + (c.beta*beta*net.den*V*(1-ve)*Ir))), ## infection of vaccinated individuals
        recover=gamma*Ir,  # recovery rate of regressive individuals
        death=mu*Ip,        # infection-induced mortality
        respawn=nu*Dt,     # "respawning" (i.e. birth) of dead individuals
        vax=tau*S/N          # vaccination rate
      )    
    })
  }else if(t<=52){ # don't include vaccination if within the first year
    with(as.list(c(y,p)),{
      c(inf=((c.rate)*((beta*net.den*S*Ip) + (c.beta*beta*net.den*S*Ir))),  ## progressive infection + regressive infection
        infv=((c.rate)*((beta*net.den*V*(1-ve)*Ip) + (c.beta*beta*net.den*V*(1-ve)*Ir))), ## infection of vaccinated individuals
        recover=gamma*Ir,  # recovery rate of regressive individuals
        death=mu*Ip,        # infection-induced mortality
        respawn=nu*Dt,     # "respawning" (i.e. birth) of dead individuals
        vax=0              # vaccination rate 
      )    
    })
  }
}

### epi transitions function
epifun <- function(y,w) { # w indicates the switch element that is observed
  switch(w,
         y + c(-1,1,0,0,0,0,0,0,0),     ## 1 = persistent/progressive infection: S-1, Ip+1, Ir+0, R+0, V+0, Dt+0, D+0, Ir_f+0
         y + c(0,1,0,0,-1,0,0,0,0),     ## 2 = vax becomes persistent/progressive infection: S+0, Ip+1, Ir+0, R+0, V-1, Dt+0, D+0, Ir_f+0
         y + c(0,0,-1,1,0,0,0,0,0),     ## 3 = regressive infection recovery: S+0, Ip+0, Ir-1, R+1, V+0, Dt+0, D+0, Ir_f+0
         y + c(0,-1,0,0,0,1,1,0,0),     ## 4 = disease-induced mortality: S+0, Ip-1, Ir+0, R+0, V+0, Dt+1, D+1, Ir_f+0 (both tracking/permanent and temporary death categories)
         y + c(1,0,0,0,0,-1,0,0,0),     ## 5 = respawning: S+1, Ip+0, Ir+0, R+0, V+0, Dt-1, D+0, Ir_f+0 (don't remove from permanent death)
         y + c(-1,0,0,0,1,0,0,0,1),     ## 6 = vaccination: S-1, Ip+0, Ir+0, R+0, V+1, Dt+0, D+0, Ir_f+0, V_f+1
         y + c(-1,0,1,0,0,0,0,1,0),     ## 7 = regressive infection: S-1, Ip+0, Ir+1, R+0, V+0, Dt+0, D+0, Ir_f+1
         y + c(-1,0,0,1,0,0,0,0,0),     ## 8 = abortive infection: S-1, Ip+0, Ir+0, R+1, V+0, Dt+0, D+0, Ir_f+0
         y + c(0,0,1,0,-1,0,0,1,0),     ## 9 = vax becomes regressive infection: S+0, Ip+0, Ir+1, R+0, V-1, Dt+0, D+0, Ir_f+1
         y + c(0,0,0,1,-1,0,0,0,0)      ## 10 = vax becomes abortive infection: S+0, Ip+0, Ir+0, R+1, V-1, Dt+0, D+0, Ir_f+0
  )
  
} 



# A wrapper function to run the simulation with specified parameters/starting values:
run <- function(p=c(beta=2, c.beta = 1, net.den = 0.1, c.rate = 0.2, gamma = 0.5, mu=0.1, nu=1, tau = 0.8, ve = 0.8, N=100),
                I0=1,
                tmax=52*2.5,
                reps = 10,
                outcome.props = c(0.25, 0.25, 0.5)) {
  
  SIRdata <- list()
  for(z in 1:reps){
   repeat{
     # set up results matrix
    epimat <- matrix(NA, nrow=1e5, ncol=11,
                     dimnames=list(NULL,c("S", "Ip", "Ir", "R", "V", "Dt", "D", "Ir_f", "V_f", "t", "it")))
    
    # set up initial conditions
    it <- 1
    t <- 0
    y <- c(S=unname(p["N"])-I0,Ip=I0, Ir=0, R=0, V=0, Dt=0, D=0, Ir_f=0, V_f=0) # initial population conditions for use in tracking number in each compartment over time
    epimat[1,] <- c(y, t, it)
    
    while ((y["Ip"] + y["Ir"])>0 & t<=tmax) {
      #### core of the Gillespie algorithm ####
      #update rates
      r <- ratefun(y,p,t)
      
      # determine time to next event; tau ~ Exp(sum(lamba_i))
      t <- t+rexp(1,rate=sum(r)) 
      
      # determine "type" of event that occurs
      w <- sample(length(r),size=1,prob=r) # probability the event is type i; p_i = lambda_i/sum(lambda_i)
      
      
      # if infection occurs, random draw for if infected individual becomes progressively, regressively, or abortively infected (1, 7, 8, respectively)
      if(w==1){ 
        w <- sample(c(1,7,8), size=1, prob=outcome.props)
      } 
      
      # if infection of a vaccinated individual occurs, random draw for type of infection:
      if(w==2){
        w <- sample(c(2,9,10), size=1, prob=outcome.props)
      }
      
      
      #########################################
      
      # update iteration number (ensures we keep the initial conditions in the results matrix)
      it <- it+1
      
      # update compartments with results of new event
      y <- epifun(y, w)
      epimat[it,] <- c(y, t, it) #epimat records numbers in compartments over time

    }
    if((tail(na.omit(epimat[,7]),1) + tail(na.omit(epimat[,8]),1)) > 1) break # if final D + Ir_f is greater than 1, keep; otherwise repeat
    }
    
    SIRdata[[z]] <- as.data.frame(epimat[!is.na(epimat[,1]),])
  }
  return(SIRdata)
}


# function for processing simulation results
# pulls results from end of outbreak or 2-year mark, whichever came first
process <- function(temp.data, c.beta){
  # duration of outbreak = time point at which there are 0 individuals with status = 1
  if(c.beta > 0){ 
    # if regressives are infectious, outbreak ends when Ip and Ir = 0, which will be last row in temp.data
    dur.time <- tail(temp.data$t, 1) # time at which Ip and Ir = 0
  }else if(c.beta == 0){
    # if regressives are NOT infectious, outbreak ends when Ip = 0
    # regressives may still exist after end of outbreak, so have to take first time in which Ip==0
    dur.time <- temp.data$t[which(temp.data$Ip==0)[1]]
  }
  
  # pull the last observation within the two year observation period
  if(max(temp.data$t)>=(104)){
    # pull last observation at or before the two year mark
    maxless <- max(temp.data$t[temp.data$t <= 104])
    last.time <- temp.data[which(temp.data$t == maxless), ]  
  }else if(max(temp.data$t)<(104)){
    last.time <- temp.data[nrow(temp.data),]
  }
  
  
  # total progressive infections (number dead at end of simulation, since only Ip can become D)
  total.prog <- last.time$D + last.time$Ip # previous Ip's plus not-yet-dead Ip's
  
  
  # total latent/regressive infections (final number Ir_f)
  total.lr <- last.time$Ir_f # Ir_f includes both recovered and not-yet-recovered Ir's
  
  
  # total abortive infections ("recovered" individuals minus latent/regressives)
  total.ab <- last.time$R-total.lr
  
  
  # total vaccinated individuals (V_f counts all EVER vaccinated)
  total.vax <- last.time$V_f
  
  
  # number of "failed" epidemics (similar to number that initiated with an isolate)
  # setting as NA because not really comparable measure for this homogeneous mixing
  num.failed <- NA
  
  
  
  
  
  temp.results <- data.frame(sim.num = NA,
                             dur.time = dur.time,
                             total.prog = total.prog,
                             total.lr = total.lr,
                             total.ab = total.ab,
                             total.vax = total.vax,
                             num.failed = num.failed
  )
  return(temp.results)
}


# add final observation at maxtime for plotting purposes
new.final <- function(temp.data, t.max){
  if(tail(temp.data$t,1)<max(t.max)){
    new.end <- tail(temp.data,1)
    new.end$t <- max(t.max)
    temp.data <- rbind(temp.data, new.end)
  }else{
    temp.data <- temp.data
  }
  return(temp.data)
}



##### LOOP THROUGH PARAMETER SETS ####
for(k in 1:length(params.to.run)){
  # set individual parameters
  param.set.num <- params.to.run[k]
  print(paste("set ", param.set.num, sep = ""))
  
  
  # pull parameters
  params <- LHS.data[LHS.data$set.id==param.set.num,]
  
  
  if(params$prop.outcomes==1){
    outcome.props <- c(0.25, 0.25, 0.50)
  } 
  
  if(params$progr.dur==1){
    progr.dur <- 18
  } else if(params$progr.dur==2){
    progr.dur <- 26
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
  
  
  ##### Reps of simulations ######
  
  ex1 <- run(p=c(beta=params$beta, c.beta = c.beta, net.den = params$net.den, c.rate = params$c.rate, 
                 gamma = 1/(progr.dur*regr.dur.c), mu=1/progr.dur, nu=params$terr.repop, 
                 tau = params$vax.rate, ve = params$vax.eff, N=params$pop.size),
             I0=1,
             tmax=52*2.5,
             reps = reps,
             outcome.props = outcome.props
             )

  
  
  
  
  ##### Store primary outcomes of interest ######
  # store 1. duration of outbreak, 2. total progressive infections, 3. total latent/regressive infections, 4. total abortive infections, 5, total vaccinated, 6. number of "failed" epidemics (NA for homogmix)
  full.set.results1 <- lapply(ex1, process, c.beta = c.beta)
  full.set.results2 <- as.data.frame(do.call(rbind, full.set.results1))
  full.set.results2$sim.num <- seq(1:nrow(full.set.results2))
  
  full.set.results2$param.set <- paste("set_", param.set.num, sep = "")
  
  # save "full.set.results2" and "ex1"
  ex1.name <- paste("Ho_Simulation_Results/ex1 data_Ho paramset ", param.set.num, "_sims ", min(full.set.results2$sim.num), "-", max(full.set.results2$sim.num), ".Rdata", sep = "")
  save(ex1, file = ex1.name)
  
  full.name <- paste("Ho_Simulation_Results/full set results_Ho paramset ", param.set.num, ".Rdata", sep = "")
  save(full.set.results2, file = full.name)
  
  
  #### process results for plotting ######
  
  # determine max number of individuals to plot
  y.max <- NULL
  for(i in 1:reps){
    temp.ex <- ex1[[i]]
    y.max <- c(y.max, max(c(max(temp.ex$Ip), max(temp.ex$R), max(temp.ex$D))))
  }
  # max(y.max)
  
  # determine max time
  t.max <- NULL
  for(i in 1:reps){
    temp.ex <- ex1[[i]]
    t.max <- c(t.max, max(temp.ex$t))
  }
  # max(t.max)
  
  
  # add final observation at maxtime for plotting purposes
  ex1_p <- lapply(ex1, new.final, t.max = t.max)
  
  #### get rolling averages for compartments based on sliding windows ####
  all.wind.avg <- NULL
  wind.width <- 5
  start.t <- 0
  wind.overlap <- 5/2
  
  maxneg <- function(x){
    max(x[x<0])
  }
  
  while((start.t+wind.width) < ceiling(max(t.max))){
    
    # update time window
    temp.wind <- c(start.t, start.t+wind.width)
    # update rep-based storage
    temp.wind.keep <- NULL
    
    # loop through each dataset in the full reps data, pulling all observations within current time window
    for(j in 1:reps){
      
      temp.ex <- ex1[[j]]
      
      # extract data within the current time window
      temp.ex2 <- temp.ex[temp.ex$t>=temp.wind[1] & temp.ex$t < temp.wind[2],]
      
      # if only 1 observation in this window, keep and move on
      if(nrow(temp.ex2)==1){
        to.keep <- temp.ex2[,c(1:10)]
        to.keep$wind.min <- temp.wind[1]
        to.keep$wind.max <- temp.wind[2]
        to.keep$wind.mean <- mean(temp.wind)
      }
      # if more than one observation per dataset in time window, take the mean of these points so have one value per dataset
      if(nrow(temp.ex2)>1){
        to.keep <- data.frame(t(colMeans(temp.ex2[,c(1:10)])))
        to.keep$wind.min <- temp.wind[1]
        to.keep$wind.max <- temp.wind[2]
        to.keep$wind.mean <- mean(temp.wind)
      }
      # if no data in the current time window, take the mean of the points on either side of the current window
      if(nrow(temp.ex2)==0){
        # get lower bound
        min.point <- temp.ex[which(temp.ex$t-temp.wind[1]==maxneg(temp.ex$t-temp.wind[1])),]
  
        
        # get upper bound
        max.point <- temp.ex
        max.point$max.diff <- max.point$t-temp.wind[2]
        max.point <- max.point[max.point$max.diff>0,]
        max.point <- max.point[which.min(max.point$max.diff),]
        max.point <- max.point[,-c(12)] # drop the "max.diff" column
        
        
        temp.ex2 <- rbind(min.point, max.point)
        # if epidemic is still going in this window, keep average of observations at upper and lower bounds
        if(nrow(temp.ex2)>1){
          to.keep <- data.frame(t(colMeans(temp.ex2[,c(1:10)])))
          to.keep$wind.min <- temp.wind[1]
          to.keep$wind.max <- temp.wind[2]
          to.keep$wind.mean <- mean(temp.wind)
        }
        # if epidemic has ended, keep the final observation from the epidemic
        if(nrow(temp.ex2)==1){
          to.keep <- tail(temp.ex[,c(1:10)], 1)
          to.keep$wind.min <- temp.wind[1]
          to.keep$wind.max <- temp.wind[2]
          to.keep$wind.mean <- mean(temp.wind)
        }
      }
      
      # save data from each rep in this time window
      temp.wind.keep <- rbind(temp.wind.keep, to.keep)
    }
    
    # take average values for each compartment in this time window and save
    wind.avg <- data.frame(t(colMeans(temp.wind.keep)))
    all.wind.avg <- rbind(all.wind.avg, wind.avg)
    
    # update start of time window
    start.t <- start.t+wind.overlap
  }
  
  
  
  #### plot all results ####
  
  # Save plotting data
  ex1p.name <- paste("Ho_Simulation_Results/ex1p data_Ho paramset ", param.set.num, "_sims ", min(full.set.results2$sim.num), "-", max(full.set.results2$sim.num), ".Rdata", sep = "")
  save(ex1_p, file = ex1p.name)
  
  # Save plot 
  eplot.name <- paste("Ho_Simulation_Results/Simulation_Figures/epi plot_Ho paramset ", param.set.num, "_sims ", min(full.set.results2$sim.num), "-", max(full.set.results2$sim.num), ".jpg", sep = "")
  jpeg(eplot.name)
  
  # par(mfrow=c(1,1))
  with(ex1_p[[1]],
       {plot(t,R,col=alpha("purple", 0.4),lty=2, type = "l", ylim=c(0,max(y.max)), xlim = c(0,max(t.max)),
             xlab = "time",
             ylab = "number of individuals")
         lines(t,Ip,col=alpha("red", 0.4),lty=2)
         lines(t,D,col = alpha("black", 0.4), lty=2)  
       })
  for (i in 2:reps){
    with(ex1_p[[i]],
         {lines(t,R,col=alpha("purple", 0.4),lty=2)
           lines(t,Ip,col=alpha("red", 0.4),lty=2)
           lines(t,D,col = alpha("black", 0.4), lty=2)  
         })
  }
  
  
  lines(all.wind.avg$wind.min, all.wind.avg$Ip,col="red",lty=1, lwd = 2)
  lines(all.wind.avg$wind.min, all.wind.avg$R,col="purple",lty=1, lwd = 2)
  lines(all.wind.avg$wind.min, all.wind.avg$D,col="black",lty=1, lwd = 2)
  
  dev.off()
}
