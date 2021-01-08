## simulate_SOnet
# 
#========================================================	
# ---
### title: Simulate spatial overlap-based network
# author: Marie Gilbertson
# date: "04/03/2020"
#---
# ##  Preamble	
# 
# What this code does:
# 1. Function for simulating a spatial overlap-based network with a given size
# spatial overlap networks were used to generate the target degree distribution
# negative binomial distribution parameters describe the shape of the degree distribution from which to draw


## NOTE: requires igraph, ergm, and intergraph packages to run


simulate_sonet <- function(pop.size = pop.size,          # size of simulated population
                           model.mu = model.mu,          # "mu" parameter for negative binomial distribution
                           model.size = model.size       # "size" parameter for negative binomial distribution
                           ){
  
  test <- rnbinom(n= pop.size, mu = model.mu, size = model.size)

  targets <- as.vector(table(test))
  rand.degs <- sort(unique(test))
  target.dens <- sum(test)/choose(pop.size, 2)
  
  # make sure simulated network is reasonably well populated
  repeat{
  x <- network(pop.size, directed = FALSE)
  y <- san(x~degree(rand.degs), target.stats = targets) # tries to maintain mean degree

  if((network.density(y)/target.dens)>0.3) break
  }

  
  # convert to igraph object
  sim.net <- asIgraph(y)
  
  return(sim.net)
  
}