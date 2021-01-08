## CuzEd_sims.R
# 
#========================================================	
# ---
### title: Cuzick and Edward's test of FeLV simulations
# author: Marie Gilbertson
# date: "12/13/2020"
#---
# ##  Preamble	
# 
# What this code does:
# 1. Runs Cuzick and Edward's test for spatial clustering across "feasible" simulations (FIV_geo & SO_geo only)

##### Clear Environment #####
remove(list=ls())


#### load libraries ####
library(smacpod) # for Cuzick and Edward's test
library(spatstat) # to make ppp for C/E test
library(sp)
library(rgdal) #for lat/long to UTM
library(igraph)
library(plyr)
library(dplyr)
library(ggplot2)
library(forcats) # for releveling in plotting
library(stringr)
library(reshape2)
library(ggpubr)



#### set seed ####
set.seed(2987)


#### functions ####

# function for calculating expected test statistic for Cuzick-Edward's
exp_Tq <- function(pop.size,
                    neighbors = c(3, 5, 7), # knowing that 3, 5, and 7 were significant neighbor levels in empirical outbreak
                    cases){
  all_Etq <- c()
  
  for(k in 1:length(neighbors)){
    temp.n <- neighbors[k]
    
    p <- (cases/pop.size)*( (cases-1)/(pop.size-1) )
    
    E_tq <- p*temp.n*pop.size
    
    all_Etq <- c(all_Etq, E_tq)
  }
  return(all_Etq)
}


## Classify cases
#0 means susceptible; 1 means infectious; 2 means regressively/latently infected; 3 means recovered/immune; 4 means dead from infection; 5 means "recovered" latent/regressive individuals; 6 means vaccinated
classify_case <- function(x){
  if(x %in% c("0", "3", "6")){
    type <- "Control"
  }else if(x %in% c("1", "2", "4", "5")){
    type <- "Case"
  }
}


#### function for classifying p-values for plotting ####
p.classify <- function(x){
  if(x=="a_0.05"){
    classify <- "p ≤ 0.05"
  }else if(x=="a_0.1"){
    classify <- "0.05 < p ≤ 0.1"
  }else if(x=="a_0.5"){
    classify <- "0.1 < p ≤ 0.5"
  }else if(x=="a_1.0"){
    classify <- "0.5 < p ≤ 1.0"
  }
}

# Cuzick-Edward's test function (for neighbor levels statistically significant in empirical outbreak)
ce.test <- function(sim.data, neighbors = c(3, 5, 7)){
  
  # convert lat/long to UTM for Cuzick/Edward's test
  spgeo <- SpatialPoints(cbind(sim.data$longitude, sim.data$latitude), proj4string=CRS("+proj=longlat +datum=WGS84"))
  sputm <- data.frame(spTransform(spgeo, CRS("+proj=utm +zone=17N +datum=WGS84")))
  colnames(sputm) <- c("UTM_East", "UTM_North")
  
  sim.data2 <- cbind(sim.data, sputm)
  
  # run the Cuzick/Edward's test
  tele.window <- owin(xrange=c(min(sim.data2$UTM_East), max(sim.data2$UTM_East)), yrange = c(min(sim.data2$UTM_North), max(sim.data2$UTM_North)))
  sim.data.pp <- ppp(x = sim.data2$UTM_East, y = sim.data2$UTM_North, window = tele.window, marks = as.factor(sim.data2$classbin))
  invisible(capture.output(runtest <- qnn.test(sim.data.pp, q = neighbors, nsim = 999, case = 2, longlat = F)))
  
  return(runtest$qsum)
}


# function to prepare data for plotting at each neighborhood level analyzed
plot.data <- function(results.data, 
                      neighbor.level = c("p3", "p5", "p7")){
  
  plot.list <- vector(mode = "list", length = length(neighbor.level))
  
  for(z in 1:length(neighbor.level)){
    temp.results.data <- results.data[,c(1, 2, which(colnames(results.data)==neighbor.level[z]))]
    colnames(temp.results.data) <- c("param.set", "sim.num", "p")
    
    plot.results <- data.frame(matrix(nrow = length(unique(results.data$param.set)), ncol = 5))
    colnames(plot.results) <- c("param.set", "a_0.05", "a_0.1", "a_0.5", "a_1.0")
    
    plot.results$param.set <- unique(temp.results.data$param.set)
    
    plot.results$a_0.05 <- ddply(temp.results.data, .(param.set), function(x) length(which(x$p <= 0.05)))[,2]
    plot.results$a_0.1 <- ddply(temp.results.data, .(param.set), function(x) length(which(x$p > 0.05 & x$p <= 0.1)))[,2]
    plot.results$a_0.5 <- ddply(temp.results.data, .(param.set), function(x) length(which(x$p > 0.1 & x$p <= 0.5)))[,2]
    plot.results$a_1.0 <- ddply(temp.results.data, .(param.set), function(x) length(which(x$p > 0.5 & x$p <= 1)))[,2]
    
    # check if each row adds up to the number of simulations (should give FALSE if correct)
    if(any(apply(plot.results, 1, function(x) sum(x[2:5]))!=length(unique(results.data$sim.num)))){
      print("Warning! Row math is wrong!")
    }
    
    # convert from wide to long format
    plot.results_long <- melt(plot.results, id.vars = "param.set")
    
    plot.results_long$p.class <- apply(plot.results_long, 1, function(x) p.classify(x[2]))
    
    plot.results_long$p.class <- factor(plot.results_long$p.class, levels = c("0.5 < p ≤ 1.0",  "0.1 < p ≤ 0.5", "0.05 < p ≤ 0.1", "p ≤ 0.05"))
    
    plot.results_long$n.level <- neighbor.level[z]
    
    
    plot.list[[z]] <- plot.results_long
    names(plot.list)[z] <- neighbor.level[z]
  }
  
  return(plot.list)
}


# function to run Cuzick/Edward's test across simulation results
ce.loop <- function(param.sets = param.sets,
                        sim.nums = sim.nums,
                        type = type,
                        use.respawns = T # NOTE: has option to ignore respawns if concerned about multiple individuals in the same location
                    ){
  #### start of C/E loop ####
  
  full.results <- NULL
  
  for(i in 1:length(param.sets)){
    
    param.set.num <- param.sets[i]
    print(param.set.num)
    
    
    temp.results <- data.frame(matrix(ncol = 11, nrow = length(sim.nums)))
    colnames(temp.results) <- c("param.set", "sim.num", "Tq3", "Tq5", "Tq7", "p3", "p5", "p7", "E3", "E5", "E7")
    
    for(j in 1:length(sim.nums)){
      #### load and reformat data ####
      sim.num <- sim.nums[j]
      
      
      g.name <- paste("FIV_Simulation_Results/g network_FIV paramset ", param.set.num, "_sim ",sim.num, ".Rdata", sep = "")
      g <- get(load(g.name))
      
      if(type=="FIV"){
        mnew.name <- paste("FIV_Simulation_Results/mnew data_FIV paramset ", param.set.num, "_sim ", sim.num, ".Rdata", sep = "")
        m.new <- get(load(mnew.name))
      }else if(type=="SO"){
        mnew.name <- paste("SO_Simulation_Results/mnew data_SO paramset ", param.set.num, "_sim ", sim.num, ".Rdata", sep = "")
        m.new <- get(load(mnew.name))
      }
      
      # extract lat and long coordinates from g
      ids <- vertex_attr(g, "vertex.names", index = V(g))
      lat <- vertex_attr(g, "latitude", index = V(g))
      long <- vertex_attr(g, "longitude", index = V(g))
      ss.data <- data.frame(id = ids,
                            latitude = lat,
                            longitude = long
      )
      ss.data$id <- as.character(ss.data$id)
      
      if(use.respawns==TRUE){
        # pull lat and long coordinates for the respawns
        
        respawns <- data.frame(id = colnames(m.new)[!colnames(m.new) %in% ss.data$id])
        respawns$id <- as.character(respawns$id)
        
        if(nrow(respawns)>0){
          o.id <- matrix(unlist(strsplit(respawns$id, ".", fixed = T)), ncol = 2, byrow = T)
          respawns$o.id <- o.id[,1]
          
          
          respawns2 <- left_join(respawns, ss.data, by = c("o.id" = "id"))
          
          respawns <- respawns2[,c("id", "latitude", "longitude")]
          
          ss.data <- rbind(ss.data, respawns)
        }
        
        # merge lat and long coordinates per individual with "outcome" from m.new
        outcome <- data.frame(id = paste(colnames(m.new)),
                              outcome = paste(tail(m.new, 1))
        )
        outcome$id <- as.character(outcome$id)
        
        if(nrow(outcome)!=nrow(ss.data)){
          print("Warning! Data mismatch!")
        }
        ss.data <- left_join(ss.data, outcome, by = "id")
        
      }else if(use.respawns==FALSE){
        ### if not including respawns, just take the "original" individuals
        
        # merge lat and long coordinates per individual with "outcome" from m.new
        outcome <- data.frame(id = paste(colnames(m.new)[1:nrow(ss.data)]),
                              outcome = paste(tail(m.new[,1:nrow(ss.data)], 1))
        )
        outcome$id <- as.character(outcome$id)
        
        if(nrow(outcome)!=nrow(ss.data)){
          print("Warning! Data mismatch!")
        }
        ss.data <- left_join(ss.data, outcome, by = "id")
        
      }
      
      # classify cases vs controls
      ss.data$class <- apply(as.matrix(ss.data[,c("outcome")]), 1, function(x) classify_case(x))
      ss.data$classbin <- ifelse(ss.data$class=="Case", 1, 0)
      
      
      #### run the C/E test ####
      suppressWarnings(ss.ce <- ce.test(sim.data = ss.data,
                       neighbors = c(3, 5, 7)))
      
      
      # calculate expected Tq
      sim_Etq <- exp_Tq(pop.size = nrow(ss.data),
                        neighbors = c(3, 5, 7),
                        cases = sum(ss.data$classbin))
      
      # store results
      temp.results$param.set[j] <- param.set.num
      temp.results$sim.num[j] <- sim.num
      temp.results[j, 3:5] <- ss.ce$Tq
      temp.results[j, 6:8] <- ss.ce$pvalue
      temp.results[j, 9:11] <- sim_Etq
    }
    
    full.results <- rbind(full.results, temp.results)
  }
  
  
  
  #### EXTRACT RESULTS FOR PLOTTING ####
  plotting <- plot.data(results.data = full.results,
                        neighbor.level = c("p3", "p5", "p7"))
  
    
  

  export.results <- plotting
  export.results[[length(plotting)+1]] <- full.results
  names(export.results)[[length(plotting)+1]] <- "full.results"
  
  return(export.results)
}



#### run C/E looping function ####

type <- "FIV" 

# load sets to evaluate
if(type=="FIV"){
  ideal.sub_set <- get(load("FIV_ideal_sets.Rdata"))
}else if(type=="SO"){
  ideal.sub_set <- get(load("SO_ideal_sets.Rdata"))
}

ideal.sets_full <- unique(ideal.sub_set$param.set)
ideal.sets <- as.numeric(str_extract(ideal.sets_full, "[0-9]+$"))


fiv.results <- ce.loop(param.sets = ideal.sets,
                           sim.nums = seq(1, 50),
                           type = "FIV_geo",
                           use.respawns = TRUE
)



type <- "SO" 

# load sets to evaluate
if(type=="FIV"){
  ideal.sub_set <- get(load("FIV_ideal_sets.Rdata"))
}else if(type=="SO"){
  ideal.sub_set <- get(load("SO_ideal_sets.Rdata"))
}

ideal.sets_full <- unique(ideal.sub_set$param.set)
ideal.sets <- as.numeric(str_extract(ideal.sets_full, "[0-9]+$"))


so.results <- ce.loop(param.sets = ideal.sets,
                          sim.nums = seq(1, 50),
                          type = "SO_geo",
                          use.respawns = TRUE
)



#### C/E STATISTIC PLOT ####

fiv.full <- fiv.results$full.results
fiv.full$model.type <- "FIV-based"

so.full <- so.results$full.results
so.full$model.type <- "Overlap-based"

ce.results_all <- rbind(fiv.full, so.full)


ce.results_all$oe_3 <- ce.results_all$Tq3/ce.results_all$E3
ce.results_all$oe_5 <- ce.results_all$Tq5/ce.results_all$E5
ce.results_all$oe_7 <- ce.results_all$Tq7/ce.results_all$E7

ce.results_3 <- ce.results_all[ce.results_all$p3<=0.1,]
ce.results_3$bin <- ifelse(ce.results_3$p3<=0.05,  "p ≤ 0.05", "0.05 < p ≤ 0.1")
ce.results_5 <- ce.results_all[ce.results_all$p5<=0.1,]
ce.results_5$bin <- ifelse(ce.results_5$p5<=0.05,  "p ≤ 0.05", "0.05 < p ≤ 0.1")
ce.results_7 <- ce.results_all[ce.results_all$p7<=0.1,]
ce.results_7$bin <- ifelse(ce.results_7$p7<=0.05,  "p ≤ 0.05", "0.05 < p ≤ 0.1")


#### Box plots per model type ####

g3_3 <- ggplot(ce.results_3, aes(x = model.type, y = oe_3)) +
  geom_jitter(aes(colour = bin), width = 0.1) +
  geom_boxplot(colour = "black", width = 0.3, outlier.shape = NA, fill = NA) +
  geom_hline(yintercept = 1.515152, colour = "red") +
  scale_colour_manual(values = c("0.05 < p ≤ 0.1" = "#41b6c4",
                                 "p ≤ 0.05" = "#225ea8")) +
  scale_y_log10(limits = c(1, 36)) +
  xlab("Model Type") +
  ylab("Observed/Expected") +
  theme(legend.title = element_blank(), axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12), legend.text=element_text(size=12))
g3_3



g3_5 <- ggplot(ce.results_5, aes(x = model.type, y = oe_5)) +
  geom_jitter(aes(colour = bin), width = 0.1) +
  geom_boxplot(colour = "black", width = 0.3, outlier.shape = NA, fill = NA) +
  geom_hline(yintercept = 1.454545, colour = "red") +
  scale_colour_manual(values = c("0.05 < p ≤ 0.1" = "#41b6c4",
                                 "p ≤ 0.05" = "#225ea8")) +
  scale_y_log10(limits = c(1, 36)) +
  # xlab("Model Type") +
  # ylab("Observed/Expected") +
  theme(legend.title = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12), legend.text=element_text(size=12)) 
g3_5


g3_7 <- ggplot(ce.results_7, aes(x = model.type, y = oe_7)) +
  geom_jitter(aes(colour = bin), width = 0.1) +
  geom_boxplot(colour = "black", width = 0.3, outlier.shape = NA, fill = NA) +
  geom_hline(yintercept = 1.396104, colour = "red") +
  scale_colour_manual(values = c("0.05 < p ≤ 0.1" = "#41b6c4",
                                 "p ≤ 0.05" = "#225ea8")) +
  scale_y_log10(limits = c(1, 36)) +
  # xlab("Model Type") +
  # ylab("Observed/Expected") +
  theme(legend.title = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12), legend.text=element_text(size=12)) 
g3_7


g3 <- ggarrange(g3_3, g3_5, g3_7, labels = c("A", "B", "C"), ncol = 3, common.legend = T) 
g3
