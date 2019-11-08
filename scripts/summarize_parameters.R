#
# Script to create a dataframe summarizing the simulation results
# Skewness and kurtosis computed with library e1071.
# Tree imbalance (colless' and sackin's indices) computed with library apTreeshape.
# Gamma statistic of Pybus and Harvey and Faith's PD computed using ape.
#

summary_path = "C:\\Users\\User\\Desktop\\March\\results\\"
final_path = "C:\\Users\\User\\Desktop\\March\\final_results\\"

sim_types = c("uniform")
#sim_types = c("static_Gaussian_speciation")
#sim_types = c("static_skewed_speciation")
#sim_types = c("static_Gaussian_extinction")
#sim_types = c("shifting_Gaussian_speciation")
#sim_types = c("shifting_neg_skewed_speciation")
#sim_types = c("shifting_pos_skewed_speciation")
#sim_types = c("shifting_Gaussian_extinction")

combined_dataframe = data.frame()

# loop through results files and build dataframe 
for ( i in 1:length(sim_types) ) {

    file = paste(summary_path, sim_types[i], ".csv", sep="")
    d = read.table(file=file, sep=",", header=TRUE)

    # simulation parameters
    birth = vector()
    death = vector()
    beta = vector()
    peak_birth = vector()
    sigma_birth = vector()
    jump = vector()
    slide = vector()
    density_dep = vector()
    sim_time = vector()

    # for each simulation type add the correct
    # parameter values
    if (sim_types[i] == "uniform") {

# "uniform"
birth <- seq(0.1,0.6,0.05)
death <- rep(0,11)
beta <-  rep(0.01,11)
      
        peak_birth = rep(0, length(birth))
        sigma_birth = rep(0, length(birth))
        jump = rep(0, length(birth))
        slide = rep(0, length(birth))
        density_dep = rep(0, length(birth))

     } else if (sim_types[i] == "static_Gaussian_speciation") {
#    } else if (sim_types[i] == "static_skewed_speciation") {
#    } else if (sim_types[i] == "static_Gaussian_extinction") {

# "static_Gaussian_speciation"
birth_p <- rep(0,7) 
death_p <- rep(0,7)
beta_p <-  c(0.005,0.01,0.05,0.1,0.5,1,5)
peak_p =  seq(0.1,0.6,0.05)
sigma_p = 0.1

# "static_skewed_speciation"
#birth_p <- rep(0,7) 
#death_p <- rep(0,7) 
#beta_p <-  c(0.005,0.01,0.05,0.1,0.5,1,5)
#peak_p =  seq(0.1,0.6,0.05)
#sigma_p = 0.1

# "static_Gaussian_extinction"
#birth_p <- rep(0.6,7) 
#death_p <- birth
#beta_p <-  c(0.005,0.01,0.05,0.1,0.5,1,5)
#peak_p =  seq(0.1,0.6,0.05)
#sigma_p = 0.1

        for (j in 1:length(birth_p)) {

            for (k in 1:length(peak_p)) {

                for (l in 1:length(sigma_p)) {

                    birth = c(birth, birth_p[j])
                    death = c(death, death_p[j])
                    beta = c(beta, beta_p[j])
                    peak_birth = c(peak_birth, peak_p[k])
                    sigma_birth = c(sigma_birth, sigma_p[l])
                    jump = c(jump, 0)
                    slide = c(slide, 0)
                    density_dep = c(density_dep, 0)
                }
            }
        }

    } else if (sim_types[i] == "shifting_Gaussian_speciation") {
#    } else if (sim_types[i] == "shifting_neg_skewed_speciation") {
#    } else if (sim_types[i] == "shifting_pos_skewed_speciation") {
#    } else if (sim_types[i] == "shifting_Gaussian_extinction") {

# "shifting_Gaussian_speciation"
birth_p <- rep(0,7) 
death_p <- rep(0,7) 
beta_p <-  c(0.005,0.01,0.05,0.1,0.5,1,5)
peak_p =  seq(0.2,0.6,0.2)
sigma_p = 0.1
move_p =  seq(0,1.0,0.02)

# "shifting_neg_skewed_speciation"
#birth_p <- rep(0,7) 
#death_p <- rep(0,7) 
#beta_p <-  c(0.005,0.01,0.05,0.1,0.5,1,5)
#peak_p =  seq(0.2,0.6,0.2)
#sigma_p = 0.1
#move_p =  seq(0,1.0,0.02)

# "shifting_pos_skewed_speciation"
#birth_p <- rep(0,7) 
#death_p <- rep(0,7) 
#beta_p <-  c(0.005,0.01,0.05,0.1,0.5,1,5)
#peak_p =  seq(0.2,0.6,0.2)
#sigma_p = 0.1
#move_p =  seq(0,1.0,0.02)

# "shifting_Gaussian_extinction"
#birth_p <- rep(0.6,7) 
#death_p <- birth
#beta_p <-  c(0.005,0.01,0.05,0.1,0.5,1,5)
#peak_p =  seq(0.2,0.6,0.2)
#sigma_p = 0.1
#move_p =  seq(0,1.0,0.02)

        for (j in 1:length(birth_p)) {

            for (k in 1:length(peak_p)) {

                for (l in 1:length(sigma_p)) {

                    for (m in 1:length(move_p)) {
                        
                        birth = c(birth, birth_p[j])
                        death = c(death, death_p[j])
                        beta = c(beta, beta_p[j])
                        peak_birth = c(peak_birth, peak_p[k])
                        sigma_birth = c(sigma_birth, sigma_p[l])
                        jump = c(jump, 0)
                        slide = c(slide, move_p[m])
                        density_dep = c(density_dep, 0)
                    }
                }
            }
        }

    } 

    # add data to dataframe
    d["speciation"] = birth 
    d["extinction"] = death
    d["beta"] = beta
    d["peak_speciation"] = peak_birth
    d["sigma_speciation"] = sigma_birth
    d["jump_at_5_Ma"] = jump
    d["slide_every_1_Ma"] = slide
    d["density_dep_extinction"] = density_dep
    d["sim_time_Ma"] = rep(10, length(birth))
    d["root_char_state"] = rep(0, length(birth))
        
    if (length(combined_dataframe) == 0)
        combined_dataframe = d
    else
        combined_dataframe = rbind(d, combined_dataframe)

    print(paste("Writing file: ", final_path, sim_types[i], ".csv", sep=""))
    write.table(d, sep=",", row.names=FALSE, file=paste(final_path, sim_types[i], ".csv", sep=""))
}

print(paste("Writing file: ", final_path, "all_combined.csv", sep=""))
write.table(combined_dataframe, sep=",", row.names=FALSE, file=paste(final_path, "all_combined.csv", sep=""))

print("Done.")