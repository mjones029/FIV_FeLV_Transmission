## Phyloscanner_to_network.R
# 
#========================================================	
# ---
### title: Phyloscanner to network
# author: Marie Gilbertson
# date: "05/05/2020"
#---
# ##  Preamble	
# 
# What this code does:
# 1. Read in pairwise relationship output from Phyloscanner and generate network object for analysis
# 2. Remove edges (inferred transmission events) between individuals that were not alive at the same time.

##### Clear Environment #####
remove(list=ls())

##### Set seed #####
set.seed(89490)


#### load libraries ####
library(igraph)
library(lubridate)


#### CREATE RAW TRANSMISSION NETWORK ####

#### load results and generate dyad support data ####

# function for comparing support for any transmission relationship to support for no transmission on a per-dyad basis
generate_dyad_support <- function(phylo.data = phylo){
  dyads <- phylo$dyad <- paste(phylo$host.1, phylo$host.2, sep = "_")
  dyads <- unique(dyads)
  
  support.data <- data.frame(dyad=dyads,
                             host.1 = as.character(NA),
                             host.2 = as.character(NA),
                             any.ancestry = as.numeric(-1),
                             no.trans = as.numeric(-1),
                             stringsAsFactors = FALSE)
  
  for(i in 1:length(dyads)){
    temp.dyad <- dyads[i]
    temp.data <- phylo[phylo$dyad==temp.dyad,]
    
    any.ancestry <- temp.data[temp.data$ancestry!="noAncestry",]
    
    support.data$any.ancestry[i] <- sum(any.ancestry$fraction.R)
    support.data$no.trans[i] <- sum(temp.data[temp.data$ancestry=="noAncestry","fraction.R"])
    
    support.data$host.1[i] <- paste(temp.data$host.1[1])
    support.data$host.2[i] <- paste(temp.data$host.2[1])
  }
  
  
  support.data$a_n.ratio <- support.data$any.ancestry/support.data$no.trans
  return(support.data)
}



# load phyloscanner results
phylo <- read.csv("Phyloscanner_networks/FL FIV_hostRelationshipSummary.csv", header = T)
phylo$fraction.R <- phylo$ancestry.tree.count/phylo$both.exist


# use function to generate data about dyad support
support.data_all <- generate_dyad_support(phylo)


# store tree support for transmission network figure
# save(support.data_all, file = "Phyloscanner_networks/FL_treesupport_all.Rdata")





#### turn phyloscanner output into a network object ####

# generate edgelist
el <- support.data_all[,c("host.1", "host.2")] 
ph.g <- graph_from_data_frame(el, directed = F)
is.weighted(ph.g)
ph.g

# check for mutualisms
ph.g <- simplify(ph.g) 
ph.g


#### REMOVE TRANSMISSION EVENTS BETWEEN NON-OVERLAPPING INDIVIDUALS ####
# remove those transmission events that are suggested between individuals that were not alive at the same time
# based on overlap of birth and death/removal dates

#### load attribute data ####
cats.fiv <- get(load("Attribute Data/FL FIV covariates.Rdata"))

# update single "kitten" to be in same category as other non-adults (i.e. yearlings)
cats.fiv$Age.categorical[cats.fiv$Age.categorical=="Kitten"] <- "Yearling"



# convert network to adjacency matrix 
tt.matrix <- as_adj(ph.g)
tt.matrix <- as.matrix(tt.matrix)

# reorder individuals to match ordering in adjacency matrix
cats.fiv2 <- cats.fiv[cats.fiv$SVRG.ID %in% rownames(tt.matrix),]
reorder_idx <- match(rownames(tt.matrix), cats.fiv2$SVRG.ID)
cats.fiv3 <- cats.fiv2[reorder_idx,]

# verify that ordering of individuals is consistent with attribute data objects
if(identical(cats.fiv3$SVRG.ID, rownames(tt.matrix))){
  cats.fiv <- cats.fiv3
}else{
  print("Not identical. Try again.")
}




#### check if transmission connections have overlapping birth/death dates ####

# loop through transmission matrix
# for every "1" in transmission matrix, check if row and column-name individuals had overlapping birth/death dates
# if so, new transmission matrix gets a "1"
# in all other cases, new transmission matrix gets a "0"

new.tt <- matrix(NA, nrow = nrow(tt.matrix), ncol = ncol(tt.matrix))
colnames(new.tt) <- colnames(tt.matrix)
row.names(new.tt) <- row.names(tt.matrix)

for(i in 1:nrow(tt.matrix)){
  for(j in 1:ncol(tt.matrix)){
    print(paste(rownames(tt.matrix)[i],colnames(tt.matrix)[j], sep = "_"))
    
    if(tt.matrix[i,j]==1){
      ind1.id <- row.names(tt.matrix)[i]
      ind2.id <- colnames(tt.matrix)[j]
      
      ind1 <- cats.fiv[cats.fiv$SVRG.ID==ind1.id,]
      ind2 <- cats.fiv[cats.fiv$SVRG.ID==ind2.id,]
      
      int1 <- interval(start = ind1$Birth.R.Date, end = ind1$Remove.R.Date)
      int2 <- interval(start = ind2$Birth.R.Date, end = ind2$Remove.R.Date)
      
      if(int_overlaps(int1, int2)){
        new.tt[i,j] <- 1
      }else{
        new.tt[i,j] <- 0
      }
    }
    
    if(tt.matrix[i,j]==0){
      new.tt[i,j] <- 0
    }
    
  }
}

# check if transmission matrices are identical
identical(tt.matrix, new.tt)


# Function to convert adjacency/transmission matrices to edgelists and evaluate which/how many edges were dropped
am.to.el <- function(t.mat){
  
  el.results <- NULL
  
  for(i in 1:nrow(t.mat)){
    for(j in 1:ncol(t.mat)){
      if(t.mat[i,j]==1){
        temp.dyad <- data.frame(ind1 = row.names(t.mat)[i],
                                ind2 = colnames(t.mat)[j])
        
        el.results <- rbind(el.results, temp.dyad)
      }
      
    }
  }

  
  return(el.results)
}


old.el <- am.to.el(tt.matrix)
new.el <- am.to.el(new.tt)

nrow(old.el)
nrow(new.el)

old.net <- graph_from_edgelist(as.matrix(old.el), directed = F)
old.net <- simplify(old.net)

new.net <- graph_from_edgelist(as.matrix(new.el), directed = F)
new.net <- simplify(new.net)


old.net
new.net

# check for isolates, once birth/death dates are accounted for

old.ids <- get.vertex.attribute(old.net, "name")
new.ids <- get.vertex.attribute(new.net, "name")

isolate.ids <- old.ids[!(old.ids %in% new.ids)]

g <- new.net %>%
  add_vertices(length(isolate.ids), name = isolate.ids)

ph.g <- g



#### VIEW NETWORK AND SAVE ####

lay <- layout.fruchterman.reingold(ph.g)


plot.igraph(ph.g,
            vertex.label=NA,
            vertex.size = 5,
            layout=lay,
            edge.color="dark grey")



# save draft network
# save(ph.g, file = "Phyloscanner_networks/FL FIV transmission network.Rdata")




