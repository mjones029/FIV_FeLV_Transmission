# Function to extract key results of interest from network transmission simulations

extract_results <- function(p.results.a = p.results, 
                            m.new.a = m.new, 
                            c.beta.a = c.beta,
                            m.list.a = m.list,
                            sim.num.a = sim.num,
                            end.point = 2          # number of years after which to pull results
                            ){
  
  # store 1. duration of outbreak, 2. total progressive infections, 3. total latent/regressive infections, 4. total abortive infections, 5, total vaccinated, 6. number of "failed" epidemics (number that initiated with an isolate)
  
  # for duration of outbreak, only include regressives in determining end of outbreak if they were infectious (i.e. if c.beta > 0)
  if(c.beta.a>0){
    # duration of outbreak = time point at which there are 0 individuals with status = 1
    dur.data <- p.results.a
    dur.time <- (which(dur.data$prop.i==0 & dur.data$prop.lr==0)[1])-1 # first instance where proportion infectious = 0 (-1 because first entry is time=0)
    # (will give "NA" if there is no time at which prop.i = 0)
  }else if(c.beta.a==0){
    # duration of outbreak = time point at which there are 0 individuals with status = 1
    dur.data <- p.results.a
    dur.time <- (which(dur.data$prop.i==0)[1])-1 # first instance where proportion infectious = 0 (-1 because first entry is time=0)
    # (will give "NA" if there is no time at which prop.i = 0)
  }
  
  # pull results from "recording endpoint" (in this case, at the end of two years or 104 weeks)
  # however, if outbreak ended before week 104, just pull results from end of m.new
  if(nrow(m.new.a)>=(end.point*52+1)){
    last.time <- m.new.a[(end.point*52+1),] # add one because row 1 = time 0 in m.new
  }else if(nrow(m.new.a)<(end.point*52+1)){
    last.time <- m.new.a[nrow(m.new.a),]
  }
  # total progressive infections (status = 1 or 4 in last row of m.new)

  total.prog <- length(which(last.time==1 | last.time==4))
  
  
  # total latent/regressive infections (status = 2 or 5 in last row of m.new)
  total.lr <- length(which(last.time==2 | last.time==5))
  
  
  # total abortive infections (status = 3 in last row of m.new)
  total.ab <- length(which(last.time==3))
  
  
  # total vaccinated individuals (stored in simulation output to account for vaccinated individuals that are infected, i.e. vax fails)
  total.vax <- sum(m.list.a[[3]])
  
  
  # number of "failed" epidemics (number that initiated with an isolate)
  # (-1 is because initiate.dur gave the number of iterations before a non-isolate was selected
  # meaning that the number of times an isolate was selected = initiate.dur-1)
  num.failed <- m.list.a[[2]]-1
  
  
  
  
  temp.results <- data.frame(sim.num = sim.num.a,
                             dur.time = dur.time,
                             total.prog = total.prog,
                             total.lr = total.lr,
                             total.ab = total.ab,
                             total.vax = total.vax,
                             num.failed = num.failed
  )
  
  
  return(temp.results)
}