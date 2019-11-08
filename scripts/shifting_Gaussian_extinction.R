#
# script to run various simulations of character evolution and
# then estimate ancestral states using brownian motion
# plus a few options for visually summarizing the results
#

###########################
# Trait dependent extinction function = 'an upside-down' normal function
negnoroptimal.x=function (x, y0, y1, xmid, s2) {
  y0 - y1 * exp(-(x - xmid)^2/(2 * s2))
}

# setup simulation parameters
replicates <- 100
max_time <- 10
root_state <- 0.0
output_file <- "shifting_Gaussian_extinction"

# parameters values to simulate under
birth <- rep(0.6,7) 
death <- birth
beta <-  c(0.005,0.01,0.05,0.1,0.5,1,5)

peak =  seq(0.2,0.6,0.2)
sigma = 0.1
move =  seq(0,1.0,0.02)

# a vector to hold the times at which each regime ends
regimes <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)

results <- list()
iter <- 1

for (i in 1:length(birth)) {

    for (j in 1:length(peak)) {
    
        for (k in 1:length(sigma)) {
        
            for (l in 1:length(move)) {

                print(paste("Simulating under parameter set:", iter, "out of", length(birth) * length(peak) * length(sigma) * length(move)))

                # for each regime, speciation and extinction are each calculated as the sum of two functions: 
                # 1) a function of the character value x
                # 2) a function of the number of lineages (diversity dependence)

                # speciation as a function of the character value x
                lambda <- list(function(x) constant.x(x,birth[i]),
                               function(x) constant.x(x,birth[i]),
                               function(x) constant.x(x,birth[i]),
                               function(x) constant.x(x,birth[i]),
                               function(x) constant.x(x,birth[i]),
                               function(x) constant.x(x,birth[i]),
                               function(x) constant.x(x,birth[i]),
                               function(x) constant.x(x,birth[i]),
                               function(x) constant.x(x,birth[i]),
                               function(x) constant.x(x,birth[i]),
                               function(x) constant.x(x,birth[i]))

                # speciation as a function of the number of lineages d
                lambda_d <- list(function(d) constant.x(d, 0.0),
                                 function(d) constant.x(d, 0.0), 
                                 function(d) constant.x(d, 0.0), 
                                 function(d) constant.x(d, 0.0), 
                                 function(d) constant.x(d, 0.0), 
                                 function(d) constant.x(d, 0.0), 
                                 function(d) constant.x(d, 0.0), 
                                 function(d) constant.x(d, 0.0), 
                                 function(d) constant.x(d, 0.0), 
                                 function(d) constant.x(d, 0.0)) 

                # extinction as a function of the character value x
                mu <- list(function(x) negnoroptimal.x(x, y0=death[i], y1=peak[j], xmid=0, s2=sigma[k]),
                           function(x) negnoroptimal.x(x, y0=death[i], y1=peak[j], xmid=1*move[l], s2=sigma[k]),
                           function(x) negnoroptimal.x(x, y0=death[i], y1=peak[j], xmid=2*move[l], s2=sigma[k]),
                           function(x) negnoroptimal.x(x, y0=death[i], y1=peak[j], xmid=3*move[l], s2=sigma[k]),
                           function(x) negnoroptimal.x(x, y0=death[i], y1=peak[j], xmid=4*move[l], s2=sigma[k]),
                           function(x) negnoroptimal.x(x, y0=death[i], y1=peak[j], xmid=5*move[l], s2=sigma[k]),
                           function(x) negnoroptimal.x(x, y0=death[i], y1=peak[j], xmid=6*move[l], s2=sigma[k]),
                           function(x) negnoroptimal.x(x, y0=death[i], y1=peak[j], xmid=7*move[l], s2=sigma[k]),
                           function(x) negnoroptimal.x(x, y0=death[i], y1=peak[j], xmid=8*move[l], s2=sigma[k]),
                           function(x) negnoroptimal.x(x, y0=death[i], y1=peak[j], xmid=9*move[l], s2=sigma[k]),
                           function(x) negnoroptimal.x(x, y0=death[i], y1=peak[j], xmid=10*move[l], s2=sigma[k]))

                # extinction as a function of the number of lineages d
                mu_d <- list(function(d) constant.x(d, 0.0),
                             function(d) constant.x(d, 0.0),
                             function(d) constant.x(d, 0.0),
                             function(d) constant.x(d, 0.0),
                             function(d) constant.x(d, 0.0),
                             function(d) constant.x(d, 0.0),
                             function(d) constant.x(d, 0.0),
                             function(d) constant.x(d, 0.0),
                             function(d) constant.x(d, 0.0),
                             function(d) constant.x(d, 0.0))

                # character evolving with brownian motion, one function per regime
                char <- list(make.brownian.with.drift(0, beta[i]),
                             make.brownian.with.drift(0, beta[i]),
                             make.brownian.with.drift(0, beta[i]),
                             make.brownian.with.drift(0, beta[i]),
                             make.brownian.with.drift(0, beta[i]),
                             make.brownian.with.drift(0, beta[i]),
                             make.brownian.with.drift(0, beta[i]),
                             make.brownian.with.drift(0, beta[i]),
                             make.brownian.with.drift(0, beta[i]),
                             make.brownian.with.drift(0, beta[i]))

                # run the simulations:
                results[[iter]] <- simulate_trees(replicates, lambda, lambda_d, mu, mu_d, char, regimes, max_time, root_state, include_extinct=FALSE)
                iter <- iter + 1

                # plot some results
                #plot_simulations(replicates, simulations, est_type="lp")

                # save the results
                save(results, file=paste("C:\\Users\\User\\Desktop\\March\\outputs\\", output_file, ".RData", sep=""))
                #save(results, file=paste(output_file, ".RData", sep=""))
                
            }

        }

    }
}