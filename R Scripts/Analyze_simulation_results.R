## Analyze_simulation_results
# 
#========================================================	
# ---
### title: Analyze simulation results
# author: Marie Gilbertson
# date: "09/17/2020"
#---
# ##  Preamble	
# 
# What this code does:
# 1. Analysis steps for evaluating model performance (i.e. how well do models predict FeLV transmission dynamics?)

##### Clear Environment #####
remove(list=ls())


#### load libraries ####
library(plyr)
library(ggplot2)
library(dplyr)


#### Set seed ####
set.seed(9346)


#### SET WHICH DATA TO LOAD ####

# model type to load
model.type <- "FIV" # options are: FIV, rand, SO, Ho

# parameter sets to load
sets <- seq(1, 150)

LHS.data <- get(load(file = "LHS parameter sets.Rdata"))
LHS.data$param.set <- paste("set", LHS.data$set.id, sep = "_")


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



######### LOAD DATA ############
all.data <- NULL

# manages different file names
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

# if an outbreak duration exceeded max time, assign max time
all.data$dur.time[is.na(all.data$dur.time)] <- 2.5*52
all.data$dur.time[all.data$dur.time>2.5*52] <- 2.5*52




#### WHICH PARAMETER SET MEDIANS MATCH OBSERVED DYNAMICS ####

median.results <- ddply(all.data, .(param.set), summarize, dur.time = median(dur.time),
      total.prog = median(total.prog), total.lr = median(total.lr), total.ab = median(total.ab))



median.results$ideal <-  apply(median.results[,c("dur.time", "total.prog", "total.lr", "total.ab")], 1, 
                             function(y) ideal.data(y[1], y[2], y[3], y[4]))

ideal.sets <- subset(median.results, median.results$ideal==1) 
ideal.sets$param.set

# save "feasible" data sets
ideal.sub_set <- all.data[all.data$param.set %in% ideal.sets$param.set,]

# save(ideal.sub_set, file = "FIV_ideal_sets.Rdata")
# save(ideal.sub_set, file = "SO_ideal_sets.Rdata")
# save(ideal.sub_set, file = "Rand_ideal_sets.Rdata")

g1 <- ggplot(ideal.sub_set, aes(x = as.factor(param.set), y = total.prog)) + geom_boxplot() +
  ylab("Total number of progressive infections") +
  xlab("Parameter set") 


g2 <- ggplot(ideal.sub_set, aes(x = as.factor(param.set), y = dur.time)) + geom_boxplot() +
  ylab("Duration (Weeks)") +
  xlab("Parameter set") 




##### RANDOM FOREST FOR PARAMETER IMPORTANCE ####
# Load libraries
library(rpart) 
library(caTools)
library(caret)
library(rpart.plot)
library(randomForest)
library(ROSE)
library(DMwR)

# Key question: within each model type, what parameters are most important for generating "ideal" predictions?


type.sets <- LHS.data[sets,]
type.sets$ideal <- NA

if(model.type=="FIV"){
  type.sets$param.set <- type.sets$param.set
}else if(model.type=="rand"){
  type.sets$param.set <- paste("rand_", type.sets$param.set, sep = "")
}else if(model.type=="SO"){
  type.sets$param.set <- paste("SO_", type.sets$param.set, sep = "")
}else if(model.type=="Ho"){
  type.sets$param.set <- type.sets$param.set
}

# assign outcome ("feasible" = 1; or not = 0)
type.sets$ideal <- ifelse(type.sets$param.set %in% ideal.sets$param.set, 1, 0)
type.sets$mod.type <- model.type


all.sets <- type.sets

#### Extract outcome and predictors of interest ####
if(model.type=="FIV"|model.type=="rand"){
  model.data <- all.sets[,c("ideal", "beta", "progr.dur", "regr.dur.c", "c.beta", "pop.size", "c.rate",
                          "terr.repop", "net.den", "vax.rate", "vax.eff")]
}else if(model.type=="SO"|model.type=="Ho"){ # don't include net.den for SO or Ho
  model.data <- all.sets[,c("ideal", "beta", "progr.dur", "regr.dur.c", "c.beta", "pop.size", "c.rate",
                            "terr.repop", "vax.rate", "vax.eff")]
}


# update c.beta to actual value for plotting purposes
for(i in 1:nrow(model.data)){
  if(model.data$c.beta[i]==1){
    model.data$c.beta[i] <- 0
  }else if(model.data$c.beta[i]==2){
    model.data$c.beta[i] <- 0.1
  }else if(model.data$c.beta[i]==3){
    model.data$c.beta[i] <- 0.5
  }else if(model.data$c.beta[i]==4){
    model.data$c.beta[i] <- 1
  }
}

# update progr.dur to actual value for plotting purposes
for(i in 1:nrow(model.data)){
  if(model.data$progr.dur[i]==1){
    model.data$progr.dur[i] <- "1/18"
  }else if(model.data$progr.dur[i]==2){
    model.data$progr.dur[i] <- "1/26"
  }
}


# update regr.dur.c to actual value for plotting purposes
for(i in 1:nrow(model.data)){
  if(model.data$regr.dur.c[i]==1){
    model.data$regr.dur.c[i] <- 0.5
  }else if(model.data$regr.dur.c[i]==2){
    model.data$regr.dur.c[i] <- 1
  }
}

# update class of factor variables
model.data$progr.dur <- as.factor(model.data$progr.dur)
model.data$regr.dur.c <- as.factor(model.data$regr.dur.c)
model.data$c.beta <- as.factor(model.data$c.beta)

model.data$ideal <- ifelse(model.data$ideal==1, "Y", "N")
model.data$ideal <- as.factor(model.data$ideal)

# for manuscript plotting 
colnames(model.data) <- c("ideal", "Beta", "Mu", "K", "C", "Pop_size", "Omega", "Nu", "Net_dens",
                          "Tau", "ve")





#### Split data into training and testing data sets ####
split  <- sample.split(model.data$ideal,SplitRatio = 0.8)

#use the split vector to select the rows to be included in each set
d.train <- model.data[split==T,]
d.test <- model.data[split==F,]


#### Build a random forest ####
rf <- randomForest(ideal~., data=d.train, type="classification", importance=T, ntree=1000)
rf

# Test model performance on the testing data.
pred<-predict(rf,d.test) 
conf <- confusionMatrix(table(pred,d.test$ideal),positive="Y")
conf


# Make an ROC plot.
roc.curve(d.test$ideal,pred)


#### balance data ####
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


#### optimize hyper-parameters #####
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


#### run final random forest ####
# Run random forest with our optimized hyper-parameters, and check fit with the testing data.
rf.final <- randomForest(ideal~.,data=d.final, 
                         mtry = optimal$mtry,
                         nodesize = optimal$nodesize,
                         importance=T)
rf.final


p <- predict(rf.final,newdata=d.test) #Predictions on Test Set for each Tree
confusionMatrix(p,d.test$ideal,positive="Y")

# plot ROC 
roc.curve(d.test$ideal,p)


#### variable importance ####
# Look at variable importance plots.
varImpPlot(rf.final)



# Need partial dependence plots to tell us about the strength and direction of effects.
#extract the importance values
imp <- importance(rf.final,type=1)

#order the variables from most to least important.
impvar <- rownames(imp)[order(imp[, 1], decreasing=TRUE)]

par(mfrow=c(3, 4))

#plot each partial dependence plot using a for loop.
for (i in seq_along(impvar)) {
  partialPlot(rf.final, d.final, impvar[i], xlab=impvar[i],
              main=paste( impvar[i]),which.class="Y")
}






#### GLMM FOR MODEL-TYPE IMPORTANCE ####
## Key question: are FIV-based networks best at generating "ideal" predictions?
# logistic mixed effects model 
# "ideal" (i.e. "feasible") as outcome, model type as explanatory variable, random intercept for parameter set


#### Load results for ALL simulations ####
# for all model types, not just a single one
all.data <- NULL

all.model.types <- c("FIV", "rand", "SO", "Ho")

# parameter sets to load
sets <- seq(1, 150)

LHS.data <- get(load(file = "LHS parameter sets.Rdata"))
LHS.data$param.set <- paste("set", LHS.data$set.id, sep = "_")

for(i in 1:length(all.model.types)){
  
  model.type <- all.model.types[i]
  
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
  
  
  set.data <- NULL
  
  for(j in 1:length(sets)){
    file.name <- paste(insert1, "Simulation_Results/full set results_", insert2, "paramset ", sets[j], ".Rdata", sep = "")
    temp.data <- get(load(file.name))
    
    set.data <- rbind(set.data, temp.data)
  }
  set.data$model.type <- model.type
  
  
  if(model.type=="Ho"){
    set.data$param.set <- paste("Ho_", set.data$param.set, sep = "")
  }
  
  all.data <- rbind(all.data, set.data)
}

# if outbreak duration was greater than max time, assign max time
all.data$dur.time[is.na(all.data$dur.time)] <- 52*2.5
all.data$dur.time[all.data$dur.time>2.5*52] <- 2.5*52


#### now (re-)determine which set medians are within the "ideal" range ####

median.results <- ddply(all.data, .(param.set), summarize, dur.time = median(dur.time),
                        total.prog = median(total.prog), total.lr = median(total.lr), total.ab = median(total.ab))



median.results$ideal <-  apply(median.results[,c("dur.time", "total.prog", "total.lr", "total.ab")], 1, 
                               function(y) ideal.data(y[1], y[2], y[3], y[4]))

ideal.sets <- subset(median.results, median.results$ideal==1) 



### merge parameters, model types, and "ideal" status ####

all.model.sets <- NULL

for(i in 1:length(all.model.types)){
  
  model.sets <- LHS.data[sets,]
  model.type <- all.model.types[i]

  if(model.type=="FIV"){
    model.sets$param.set <- model.sets$param.set
  }else if(model.type=="rand"){
    model.sets$param.set <- paste("rand_", model.sets$param.set, sep = "")
  }else if(model.type=="SO"){
    model.sets$param.set <- paste("SO_", model.sets$param.set, sep = "")
  }else if(model.type=="Ho"){
    model.sets$param.set <- paste("Ho_", model.sets$param.set, sep = "")
  }
  model.sets$model.type <- model.type
  all.model.sets <- rbind(all.model.sets, model.sets)
}

all.model.sets$ideal <- NA
all.model.sets$ideal <- ifelse(all.model.sets$param.set %in% ideal.sets$param.set, 1, 0)

sum(all.model.sets$ideal)


#### convert data classes as needed ####

# set.id to factor
all.model.sets$set.id <- as.factor(all.model.sets$set.id)

# model.type to factor
all.model.sets$model.type <- as.factor(all.model.sets$model.type)

# set ideal to factor
all.model.sets$ideal <- as.factor(all.model.sets$ideal)



#### fit GLMM ####

library(lme4)
library(lmerTest)

all.model.sets <- within(all.model.sets, model.type <- relevel(model.type, ref = "Ho"))
mod1 <- glmer(ideal ~ model.type + (1|set.id), family = binomial, data = all.model.sets)
summary(mod1)
plot(mod1)


