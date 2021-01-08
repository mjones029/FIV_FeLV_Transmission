## SatScan_sims
# 
#========================================================	
# ---
### title: SatScan local cluster analysis of simulation datas
# author: Marie Gilbertson
# date: "03/02/2020"
#---
# ##  Preamble	
# 
# What this code does:
# 1. Runs SatScan cluster analysis of simulated FeLV outbreak 
# based on: https://www.satscan.org/rsatscan/rsatscan.html

##### Clear Environment #####
remove(list=ls())


#### load libraries ####
detach("package:rsatscan", unload=TRUE) # unload if already attached
library(rsatscan)
library(igraph)
library(dplyr)
library(ggplot2)
library(forcats) # for releveling in plotting



#### set seed ####
set.seed(2987)


#### update set of files to read in ####

# parameter sets to load
param.sets <- seq(1, 150)

# sims to load
sim.nums <- seq(1, 50)

# model type
type <- "FIV" # options are: "FIV", "SO"


## Function for creating datasets of appropriate format for SatScan analysis
#0 means susceptible; 1 means infectious; 2 means regressively/latently infected; 3 means recovered/immune; 4 means dead from infection; 5 means "recovered" latent/regressive individuals; 6 means vaccinated
classify_case <- function(x){
  if(x %in% c("0", "3", "6")){
    type <- "Control"
  }else if(x %in% c("1", "2", "4", "5")){
    type <- "Case"
  }
}

#### GENERATE FILES FOR SATSCAN ####
for(i in 1:length(param.sets)){
  print(i)
  param.set.num <- param.sets[i]
  
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
    
    ## Create datasets of appropriate format for SatScan analysis
    #0 means susceptible; 1 means infectious; 2 means regressively/latently infected; 3 means recovered/immune; 4 means dead from infection; 5 means "recovered" latent/regressive individuals; 6 means vaccinated
    
    ss.data$class <- apply(as.matrix(ss.data[,c("outcome")]), 1, function(x) classify_case(x))
    cases <- ss.data[ss.data$class=="Case",]
    cases$cases <- rep(1, nrow(cases))
    cases <- cases[,c("id", "cases")]
    
    controls <- ss.data[ss.data$class=="Control",]
    controls$cases <- rep(1, nrow(controls))
    controls <- controls[,c("id", "cases")]
    
    geo <- ss.data[,c("id", "latitude", "longitude")]
    
    
    #### reset parameter file ####
    
    invisible(ss.options(reset=TRUE))
    
    
    
    #### write the SatScan files ####
    # write parameter file, case file, geometry file
    
    ### NOTE: UPDATE TD AS NEEDED FOR PERSONAL FILE NAMING SYSTEM
    td = "SaTScan Files"
    write.cas(cases, td, paste(type, "_paramset ", param.set.num, "_sim ", sim.num, sep = ""), userownames=F)
    write.ctl(controls, td, paste(type, "_paramset ", param.set.num, "_sim ", sim.num, sep = ""), userownames=F)
    write.geo(geo, td, paste(type, "_paramset ", param.set.num, "_sim ", sim.num, sep = ""), userownames=F)
    
    

    #### update SatScan parameters ####
    ss.options(list(CaseFile=paste(td, "/", type, "_paramset ", param.set.num, "_sim ", sim.num, ".cas", sep = ""), ControlFile=paste(td, "/", type, "_paramset ", param.set.num, "_sim ", sim.num,".ctl", sep = "")))
    ss.options(list(PrecisionCaseTimes=0))
    ss.options(list(CoordinatesFile=paste(td, "/", type, "_paramset ", param.set.num, "_sim ", sim.num, ".geo", sep = ""), 
                    CoordinatesType=1, AnalysisType=1, ModelType=1, ScanAreas=1))
    ss.options(list(TimeAggregationUnits=0, NonCompactnessPenalty=0, MaxSpatialSizeInPopulationAtRisk=50,
                    SpatialWindowShapeType=0, IncludePurelySpatial="y"))
    ss.options(list(ReportGiniClusters="n", LogRunToHistoryFile="n"))
    
    
    
    
    
    
    #### write out parameters ####
    write.ss.prm(td, paste(type, "_paramset ", param.set.num, "_sim ", sim.num, sep = ""))
    

  }
}

#### RUN SATSCAN EXTERNALLY ####
# in MacOS command line, can run the following in a bash shell script
#### for j in {1..150}; do
#### for i in {1..50}; do /Volumes/Macintosh\ HD/Applications/SaTScan/satscan FIV_paramset\ ${j}_sim\ ${i}.prm; done
#### done
# (if in shell script, execute the shell with "bash" command)
# j loop for parameter sets; i loop for sims
# once complete, proceed to extracting relevant results



#### READ IN, EXTRACT, AND SAVE RELEVANT RESULTS ####
remove(list=ls())


library(stringr)
library(igraph)
library(dplyr)
library(ggplot2)
library(forcats) # for releveling in plotting

#### set seed ####
set.seed(2987)

type <- "FIV" # options are: "FIV", "SO"

# load sets to evaluate
ideal.sub_set <- get(load("FIV_ideal_sets.Rdata"))
# ideal.sub_set <- get(load("SO_ideal_sets.Rdata"))

sets <- unique(ideal.sub_set$param.set)
param.sets <- as.numeric(str_extract(sets, "[0-9]+$"))


# sims to load
sim.nums <- seq(1, 50)


# loop through SaTScan results
full.satscan.results <- NULL

for(i in 1:length(param.sets)){
  param.set.num <- param.sets[i]
  temp.set.results <- NULL

  
  for(j in 1:length(sim.nums)){
    #### load and reformat data ####
    sim.num <- sim.nums[j]
    
    ### NOTE: UPDATE TD AS NEEDED
    td = "SaTScan Files"
    
    f <- readLines(paste(td, "/", type, "_paramset ", param.set.num, "_sim ", sim.num, ".txt", sep = ""), n = 100)
    
    # NOTE: commented are examples of data in original text file
    # Coordinates / radius..: (26.191001 N, 81.264783 W) / 14.66 km"
    crline <- grep("Coordinates / radius", f, value = T)[1]
    crline2 <- sub(".*: ", "", crline)
    crval1 <- sub(" /.*", "", crline2)
    crval2 <- sub(".*/ ", "", crline2)
    
    # Population............: 17"
    popline <- grep("Population", f, value = T)[1]
    popval <- as.numeric(str_extract(popline,".[0-9]+$"))
    
    # Number of cases.......: 12"  
    cline <- grep("Number of cases", f, value = T)[1]
    cval <- as.numeric(str_extract(cline,".[0-9]+$"))
    
    # Expected cases........: 5.29" 
    eline <- grep("Expected cases", f, value = T)[1]
    eval <- as.numeric(str_extract(eline,"[0-9]+\\.[0-9]+$"))
    
    # Relative risk.........: 2.88"
    rline <- grep("Relative risk", f, value = T)[1]
    rval <- as.numeric(str_extract(rline,"[0-9]+\\.[0-9]+$"))
    
    # Percent cases in area.: 70.6"
    pcline <- grep("Percent cases in area", f, value = T)[2] # take only the second instance (the first is from the summary data)
    pcval <- as.numeric(str_extract(pcline,"[0-9]+\\.[0-9]+$"))
    
    # Log likelihood ratio..: 6.659883"
    lline <- grep("Log likelihood ratio", f, value = T)[1]
    lval <- as.numeric(str_extract(lline,"[0-9]+\\.[0-9]+$"))
    
    # takes the line for p-value and extracts the actual value
    pline <- grep("P-value",f,value=TRUE)[1]
    pval <- as.numeric(str_extract(pline,"[0-9]+\\.[0-9]+$"))
    
    
    
    # assemble into dataframe
    temp.results <- data.frame(coordinates = crval1,
                               radius = crval2,
                               population = popval,
                               no.cases = cval,
                               exp.cases = eval,
                               rel.risk = rval,
                               perc.cases = pcval,
                               ll.ratio = lval,
                               p.value = pval
                               )
    temp.results$model.type <- type
    temp.results$param.set <- param.set.num
    temp.results$sim.num <- sim.num
    
    temp.set.results <- rbind(temp.set.results, temp.results)

  }
  
  full.satscan.results <- rbind(full.satscan.results, temp.set.results)
}

# classify SaTScan results by statistically significant clustering
cl.results <- full.satscan.results
cl.results$p.value[is.na(cl.results$p.value)] <- 10

p.classify <- function(x){
  if(x<=0.05){
    classify <- "p ≤ 0.05"
  }else if(x>0.05 & x<=0.1){
    classify <- "0.05 < p ≤ 0.1"
  }else if(x>0.1 & x<=1){
    classify <- "p > 0.1"
  }else if(x==10){
    classify <- "No Clusters"
  }
}
cl.results$p.class <- apply(as.matrix(cl.results[,c("p.value")]), 1, function(x) p.classify(x))

fiv.cl.results <- cl.results




#### repeat extraction for SO model type ####
type <- "SO" # options are: "FIV", "SO"

# load sets to evaluate
# ideal.sub_set <- get(load("FIV_ideal_sets.Rdata"))
ideal.sub_set <- get(load("SO_ideal_sets.Rdata"))

sets <- unique(ideal.sub_set$param.set)
param.sets <- as.numeric(str_extract(sets, "[0-9]+$"))


# sims to load
sim.nums <- seq(1, 50)


# loop through SaTScan results
full.satscan.results <- NULL

for(i in 1:length(param.sets)){
  param.set.num <- param.sets[i]
  temp.set.results <- NULL
  
  
  for(j in 1:length(sim.nums)){
    #### load and reformat data ####
    sim.num <- sim.nums[j]
    
    ### NOTE: UPDATE TD AS NEEDED FOR PERSONAL FILE NAMING SYSTEM
    td = "SaTScan Files"
    
    f <- readLines(paste(td, "/", type, "_paramset ", param.set.num, "_sim ", sim.num, ".txt", sep = ""), n = 100)
    
    # again, commented are examples from SatScan results text file
    # Coordinates / radius..: (26.191001 N, 81.264783 W) / 14.66 km"
    crline <- grep("Coordinates / radius", f, value = T)[1]
    crline2 <- sub(".*: ", "", crline)
    crval1 <- sub(" /.*", "", crline2)
    crval2 <- sub(".*/ ", "", crline2)
    
    # Population............: 17"
    popline <- grep("Population", f, value = T)[1]
    popval <- as.numeric(str_extract(popline,".[0-9]+$"))
    
    # Number of cases.......: 12"  
    cline <- grep("Number of cases", f, value = T)[1]
    cval <- as.numeric(str_extract(cline,".[0-9]+$"))
    
    # Expected cases........: 5.29" 
    eline <- grep("Expected cases", f, value = T)[1]
    eval <- as.numeric(str_extract(eline,"[0-9]+\\.[0-9]+$"))
    
    # Relative risk.........: 2.88"
    rline <- grep("Relative risk", f, value = T)[1]
    rval <- as.numeric(str_extract(rline,"[0-9]+\\.[0-9]+$"))
    
    # Percent cases in area.: 70.6"
    pcline <- grep("Percent cases in area", f, value = T)[2] # take only the second instance (the first is from the summary data)
    pcval <- as.numeric(str_extract(pcline,"[0-9]+\\.[0-9]+$"))
    
    # Log likelihood ratio..: 6.659883"
    lline <- grep("Log likelihood ratio", f, value = T)[1]
    lval <- as.numeric(str_extract(lline,"[0-9]+\\.[0-9]+$"))
    
    # takes the line for p-value and extracts the actual value
    pline <- grep("P-value",f,value=TRUE)[1]
    pval <- as.numeric(str_extract(pline,"[0-9]+\\.[0-9]+$"))
    
    
    
    # assemble into dataframe
    temp.results <- data.frame(coordinates = crval1,
                               radius = crval2,
                               population = popval,
                               no.cases = cval,
                               exp.cases = eval,
                               rel.risk = rval,
                               perc.cases = pcval,
                               ll.ratio = lval,
                               p.value = pval
    )
    temp.results$model.type <- type
    temp.results$param.set <- param.set.num
    temp.results$sim.num <- sim.num
    
    temp.set.results <- rbind(temp.set.results, temp.results)
    
  }
  
  full.satscan.results <- rbind(full.satscan.results, temp.set.results)
}

# classify SaTScan results by statistically significant clustering
cl.results <- full.satscan.results
cl.results$p.value[is.na(cl.results$p.value)] <- 10

cl.results$p.class <- apply(as.matrix(cl.results[,c("p.value")]), 1, function(x) p.classify(x))




cl.results <- rbind(fiv.cl.results, cl.results)


#### PLOT SATSCAN RESULTS ####

# supplementary labels for the following plots
supp.labs <- c("FIV-based", "Overlap-based")
names(supp.labs) <- unique(cl.results$model.type)


#### cluster size plots ####

# convert "radius" to numeric (without "km")
cl.results$radius_num <- as.numeric(gsub(" km", "", cl.results$radius))
cl.results$radius_num[is.na(cl.results$radius_num)] <- 0


c2 <- cl.results %>%
  mutate(p.class = fct_relevel(p.class, "No Clusters", "p > 0.1", "0.05 < p ≤ 0.1", "p ≤ 0.05")) %>%
  ggplot(aes(x = model.type, y = radius_num)) +
  geom_violin(data = cl.results[cl.results$p.value <= 0.1,], aes(x = model.type, y = radius_num), colour = "black", scale = "count") +
  geom_boxplot(data = cl.results[cl.results$p.value <= 0.1,], aes(x = model.type, y = radius_num), colour = "black", width = 0.1, outlier.shape = NA, fill = "#d9d9d9") +
  scale_colour_manual(values = c("No Clusters" = "#ffffcc",
                                 "p > 0.1" = "#a1dab4",
                                 "0.05 < p ≤ 0.1" = "#41b6c4",
                                 "p ≤ 0.05" = "#225ea8"))+
  geom_jitter(data = cl.results[cl.results$p.value <= 0.1,], aes(x = model.type, y = radius_num, colour=as.factor(p.class)), width = 0.1) +
  scale_x_discrete(labels = c("FIV-based", "Overlap-based")) +
  ylab("Cluster Radius (km)") +
  theme(legend.title = element_blank(), axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12), legend.text=element_text(size=12)) +
  geom_hline(yintercept = 19.87, colour = "red")


c2


#### O/E plots ####
sub.clresults <- cl.results[cl.results$p.value <= 0.1,]
sub.clresults$oe <- sub.clresults$no.cases/sub.clresults$exp.cases

oe2 <- ggplot(data = sub.clresults, aes(x = model.type, y = oe)) +
  geom_boxplot(colour = "black", width = 0.3, outlier.shape = NA) +
  
  geom_jitter(aes(colour=as.factor(p.class)), width = 0.1) +
  scale_colour_manual(values = c("0.05 < p ≤ 0.1" = "#41b6c4",
                                 "p ≤ 0.05" = "#225ea8")) +
  scale_x_discrete(labels = c("FIV-based", "Overlap-based")) +
  ylab("Observed/Expected Cases") +
  theme(legend.title = element_blank(), axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12), legend.text=element_text(size=12)) +
  geom_hline(yintercept = 2.21, colour = "red")


oe2



#### joint cluster size and O/E plots ####
library(ggpubr)

# make c2 and oe2 plots above
b <- ggarrange(c2, oe2, labels = c("A", "B"), ncol = 2, common.legend = T, widths = c(1.5, 1)) 
b

