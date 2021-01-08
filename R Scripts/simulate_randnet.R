## simulate_randnet
# 
#========================================================	
# ---
### title: Simulate random network
# author: Marie Gilbertson
# date: "04/02/2020"
#---
# ##  Preamble	
# 
# What this code does:
# 1. Function for simulating a random network with a given density and size

## NOTE: requires igraph package to run
# this isn't really a necessary function, since just uses single function from igraph, but 
# using for now for consistency

simulate_randnet <- function(net.den = net.den,           # network density to constrain simulated network to
                             pop.size = pop.size          # size of simulated population
                             ){
  
  sim.net <- erdos.renyi.game(pop.size, net.den, direct=F)    #random graph generated
  
  return(sim.net)
  
}