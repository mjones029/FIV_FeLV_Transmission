## RF_cbeta_figure
# 
#========================================================	
# ---
### title: Random forest figure for c.beta variable (Figure S10)
# author: Marie Gilbertson
# date: "11/16/2020"
#---
# ##  Preamble	
# 
# What this code does:
# 1. Variable  importance of "C" beta modification variable
# Specifically, generates supplementary Figure S10
# Adapted from random forest analysis in "Analyze_simulation_results.R"

##### Clear Environment #####
remove(list=ls())


#### load libraries ####
library(plyr)
library(ggplot2)
library(dplyr)
library(rpart) 
library(caTools)
library(caret)
library(rpart.plot)
library(randomForest)
library(ROSE)
library(DMwR)




#### LOOP PARAMETERS ####
model.types <- c("FIV", "rand", "SO", "Ho")
model.names <- c("FIV", "Random", "Overlap-based", "Homogeneous")

# parameter sets to load
sets <- seq(1, 150)

LHS.data <- get(load(file = "LHS parameter sets.Rdata"))
LHS.data$param.set <- paste("set", LHS.data$set.id, sep = "_")


# update c.beta to actual value for plotting purposes
for(i in 1:nrow(LHS.data)){
  if(LHS.data$c.beta[i]==1){
    LHS.data$c.beta[i] <- 0
  }else if(LHS.data$c.beta[i]==2){
    LHS.data$c.beta[i] <- 0.1
  }else if(LHS.data$c.beta[i]==3){
    LHS.data$c.beta[i] <- 0.5
  }else if(LHS.data$c.beta[i]==4){
    LHS.data$c.beta[i] <- 1
  }
}


#### SET IDEAL BOUNDS ####

### set ideal metrics for medians ###
dur.range <- c(78, 117)
prog.range <- c(5, 20)
lr.range <- c(0, 200)
ab.range <- c(5, 200)


#### Load function for updating "ideal" metrics ####
ideal.data <- function(dur.time, total.prog, total.lr, total.ab,
                       dur.range. = dur.range,
                       prog.range. = prog.range,
                       lr.range. = lr.range,
                       ab.range. = ab.range
){
  
  if(
    (dur.time >= dur.range[1] & dur.time <= dur.range[2]) &
    (total.prog >= prog.range[1] & total.prog <= prog.range[2]) &
    (total.lr >= lr.range[1] & total.lr <= lr.range[2]) &
    (total.ab >= ab.range[1] & total.ab <= ab.range[2])
  ){
    ideal <- 1
  }else{
    ideal <- 0
  }
  return(ideal)
}






#### START OF MAIN LOOP ####
par(mfrow = c(2, 2))

for(z in 1:length(model.types)){

  #### Set seed ####
  set.seed(9346)
  
  #### SET WHICH DATA TO LOAD ####
  
  # model type to load
  model.type <- model.types[z]
  print(model.type)
  
  ######### LOAD DATA ############
  all.data <- NULL
  
  if(model.type=="FIV"){
    insert1 <- insert2 <- ""
  }else if(model.type=="rand"){
    insert1 <- "Rand_"
    insert2 <- "rand "
  }else if(model.type=="SO"){
    insert1 <- "SO_"
    insert2 <- "SO "
  }else if(model.type=="Ho"){
    insert1 <- "Ho_"
    insert2 <- "Ho "
  }
  
  for(i in 1:length(sets)){
    
    file.name <- paste(insert1, "Simulation_Results/full set results_", insert2, "paramset ", sets[i], ".Rdata", sep = "")
    temp.data <- get(load(file.name))
    
    all.data <- rbind(all.data, temp.data)
  }
  
  all.data$dur.time[is.na(all.data$dur.time)] <- 2.5*52
  all.data$dur.time[all.data$dur.time>2.5*52] <- 2.5*52
  
  
  
  
  
  #### WHICH SET MEDIANS MATCH OBSERVED DYNAMICS ####
  
  median.results <- ddply(all.data, .(param.set), summarize, dur.time = median(dur.time),
                          total.prog = median(total.prog), total.lr = median(total.lr), total.ab = median(total.ab))
  
  
  
  median.results$ideal <-  apply(median.results[,c("dur.time", "total.prog", "total.lr", "total.ab")], 1, 
                                 function(y) ideal.data(y[1], y[2], y[3], y[4]))
  
  ideal.sets <- subset(median.results, median.results$ideal==1) 
  ideal.sets$param.set
  
  ideal.sub_set <- all.data[all.data$param.set %in% ideal.sets$param.set,]
  

  ##### RANDOM FOREST FOR PARAMETER IMPORTANCE ####
  type.sets <- LHS.data[sets,]
  type.sets$ideal <- NA
  
  if(model.type=="FIV"){
    type.sets$param.set <- type.sets$param.set
  }else if(model.type=="rand"){
    type.sets$param.set <- paste("rand_", type.sets$param.set, sep = "")
  }else if(model.type=="SO"){
    type.sets$param.set <- paste("SO_", type.sets$param.set, sep = "")
    type.sets$net.den <- NA
  }else if(model.type=="Ho"){
    type.sets$param.set <- type.sets$param.set
  }
  
  type.sets$ideal <- ifelse(type.sets$param.set %in% ideal.sets$param.set, 1, 0)
  type.sets$mod.type <- model.type
  
  
  all.sets <- type.sets
  
  # Extract outcome and predictors of interest.
  if(model.type=="FIV"|model.type=="rand"){
    model.data <- all.sets[,c("ideal", "beta", "progr.dur", "regr.dur.c", "c.beta", "pop.size", "c.rate",
                              "terr.repop", "net.den", "vax.rate", "vax.eff")]
  }else if(model.type=="SO"|model.type=="Ho"){ # don't include net.den for SO or Ho
    model.data <- all.sets[,c("ideal", "beta", "progr.dur", "regr.dur.c", "c.beta", "pop.size", "c.rate",
                              "terr.repop", "vax.rate", "vax.eff")]
  }
  
  model.data$progr.dur <- as.factor(model.data$progr.dur)
  model.data$regr.dur.c <- as.factor(model.data$regr.dur.c)
  model.data$c.beta <- as.factor(model.data$c.beta)
  
  model.data$ideal <- ifelse(model.data$ideal==1, "Y", "N")
  model.data$ideal <- as.factor(model.data$ideal)
  
  #### repeat random forest until area under the curve is greater than 0.8 ####
  trial <- 1
  repeat{
    # Split data into training and testing data sets.
    split  <- sample.split(model.data$ideal,SplitRatio = 0.8)
    
    #use the split vector to select the rows to be included in each set
    d.train <- model.data[split==T,]
    d.test <- model.data[split==F,]
    
    
    # Build a random forest
    rf <- randomForest(ideal~., data=d.train, type="classification", importance=T, ntree=1000)
    rf
    
    # Test model performance on the testing data.
    pred<-predict(rf,d.test) 
    conf <- confusionMatrix(table(pred,d.test$ideal),positive="Y")
    conf
    
    
    # Make an ROC plot.
    # roc.curve(d.test$ideal,pred)
    
    
    # Data imbalance can create problems.
    tab <- table(d.train$ideal)
    tab
    
    
    
    # We can try down-sampling to balance the data.
    d.down <- ovun.sample(ideal ~ ., data = d.train, method = "under",N=min(tab)*2)$data
    table(d.down$ideal)
    
    
    # Now we can try up-sampling, another strategy to balance the data.
    d.up <- ovun.sample(ideal ~ ., data = d.train, method = "over",N=max(tab)*2)$data
    table(d.up$ideal)
    
    
    # Another strategy is to use both up and down sampling.
    d.both <- ovun.sample(ideal ~ ., data = d.train, method = "both")$data
    table(d.both$ideal)
    
    # try SMOTE sampling.
    d.smote <- SMOTE(ideal ~ ., data = d.train)
    table(d.smote$ideal)
    
    
    
    # Now, let's run a random forest on each of the balanced datasets. We always test our data on unbalanced data.
    rf.down <- randomForest(ideal~.,data=d.down,type="classification",importance=T,ntree=1000)
    rf.up <- randomForest(ideal~.,data=d.up,type="classification",importance=T,ntree=1000)
    rf.both <- randomForest(ideal~.,data=d.both,type="classification",importance=T,ntree=1000)
    rf.smote <- randomForest(ideal~.,data=d.smote,type="classification",importance=T,ntree=1000)
    
    # the original random forest on the unbalanced data
    or <- confusionMatrix(predict(rf,newdata=d.test),
                          d.test$ideal,positive="Y")
    
    # the down sampled data
    ds <- confusionMatrix(predict(rf.down,newdata=d.test) ,
                          d.test$ideal,positive="Y")
    
    # the up sampled data
    us <- confusionMatrix(predict(rf.up,newdata=d.test),
                          d.test$ideal,positive="Y")
    
    # the up-down sampled data
    ud <- confusionMatrix(predict(rf.both,newdata=d.test),
                          d.test$ideal,positive="Y")
    
    # the SMOTE data
    sm <- confusionMatrix(predict(rf.smote,newdata=d.test),
                          d.test$ideal,positive="Y")
    
    
    # set final data based on results of balancing
    best.balance <- which.max(c(or$byClass["Balanced Accuracy"],
                                ds$byClass["Balanced Accuracy"],
                                us$byClass["Balanced Accuracy"],
                                ud$byClass["Balanced Accuracy"],
                                sm$byClass["Balanced Accuracy"]))
    
    
    
    if(best.balance==1){
      d.final <- d.train
    }else if(best.balance==2){
      d.final <- d.down
    }else if(best.balance==3){
      d.final <- d.up
    }else if(best.balance==4){
      d.final <- d.both
    }else if(best.balance==5){
      d.final <- d.smote
    }
    
    
    
    
    
    # Let's find the hyper-parameters of the random forest algorithm that result in the best-performing model. Start with creating the grid of parameters to try.
    mtry <- seq(1, ncol(d.train) * 0.8, 1)
    nodesize <- seq(1, 8, 1)
    
    # Create a data frame containing all combinations 
    hyper_grid2 <- expand.grid(mtry = mtry, nodesize = nodesize, oob_err = 0,bal.acc=0)
    
    
    # Loop through the different hyper-parameters using our "final" balanced data, and store the performance of each model.
    set.seed(1000)
    for (i in 1:nrow(hyper_grid2)) {
      # Train a Random Forest model
      model <- randomForest(ideal~.,data=d.final, 
                            mtry = hyper_grid2$mtry[i],importance=F,
                            nodesize = hyper_grid2$nodesize[i])
      
      # Store OOB error for the model                      
      hyper_grid2$oob_err[i] <- model$err.rate[nrow(model$err.rate), "OOB"]
      p <- predict(model,newdata=d.test) #Predictions on Test Set for each Tree
      ba <- confusionMatrix(p,d.test$ideal,positive="Y")$byClass["Balanced Accuracy"]#record the balanced accuracy
      hyper_grid2$bal.acc[i] <- ba 
      
    }
    
    # Find the settings that worked the best.
    opt_i <- which.max(hyper_grid2$bal.acc)
    optimal <- hyper_grid2[opt_i,]
    
    
    
    # Run random fit with our optimized hyper-parameters, and check fit with the testing data.
    rf.final <- randomForest(ideal~.,data=d.final, 
                             mtry = optimal$mtry,
                             nodesize = optimal$nodesize,
                             importance=T)
    rf.final
    
    
    p <- predict(rf.final,newdata=d.test) #Predictions on Test Set for each Tree
    confusionMatrix(p,d.test$ideal,positive="Y")
    
    # plot ROC 
    # roc.curve(d.test$ideal,p)
    auc <- roc.curve(d.test$ideal,p, plotit = F) 

    print(paste(trial, auc$auc, sep = "; "))

    trial <- trial + 1

    if(auc$auc >= 0.8) break
  }
  
  # Need partial dependence plots to tell us about the strength and direction of effects.
  #extract the importance values
  imp <- importance(rf.final,type=1)
  
  #order the variables from most to least important.
  impvar <- rownames(imp)[order(imp[, 1], decreasing=TRUE)]

  partialPlot(rf.final, d.final, "c.beta", xlab="C Value",
              main=paste(model.names[z]),which.class="Y",
              cex.lab=1.5, cex.axis = 1.5, cex.names = 1.5, cex.main = 2)

}

#### END OF MAIN LOOP ####


