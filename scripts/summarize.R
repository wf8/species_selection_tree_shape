#
# Script to create a dataframe summarizing the simulation results
# Skewness and kurtosis computed with library moments.
# Tree imbalance (colless' and sackin's indices) computed with library apTreeshape.
# Gamma statistic of Pybus and Harvey and Faith's PD computed using ape.
#

file_path = "C:\\Users\\User\\Desktop\\March\\outputs\\"
summary_path = "C:\\Users\\User\\Desktop\\March\\results\\"

sim_types = c("uniform")
#sim_types = c("static_Gaussian_speciation")
#sim_types = c("static_skewed_speciation")
#sim_types = c("static_Gaussian_extinction")
#sim_types = c("shifting_Gaussian_speciation")
#sim_types = c("shifting_neg_skewed_speciation")
#sim_types = c("shifting_pos_skewed_speciation")
#sim_types = c("shifting_Gaussian_extinction")

library(moments)
library(apTreeshape)

# function to return a vector of the number of nodes
# from each tip to the root
get_node_dist = function(tree) {
    node_dist = vector()
    n_tips = length( tree$tip.label )
    t = reorder(tree, "postorder")
    for (i in 1:n_tips) {
        next_node = i
        distance = 0
        for (j in 1:length(t$edge[,1])) {
            begin = t$edge[j, 1]
            end = t$edge[j, 2]
            if (end == next_node) {
                distance = distance + 1
                next_node = begin
            }
        }
        node_dist = c(node_dist, distance)
    }
    return( node_dist )
}

# loop through results files and build dataframe 
for ( i in 1:length(sim_types) ) {

    load( paste(file_path, sim_types[i], ".RData", sep="") )

    # simulation parameters
    sim_type = list()
    birth = list()
    death = list()
    beta = list()
    ppeak_birth = list()
    sigma_birth = list()
    jump = list()
    slide = list()
    density_dep = list()
    sim_time = list()

    d = data.frame()

    # loop through each simulation scenario
    for ( j in 1:length(results) ) {

        d[j, "sim_type"] = sim_types[i]
        d[j, "sim_hours"] = as.numeric( results[j][[1]]$processing_time, units="hours" )
        d[j, "n_replicates"] = length( results[j][[1]]$data )

        # tree shape stats
        n_too_many = 0
        n_extinct = 0
        n_error = 0
        n_single = 0 # record trees with a single surviving tip
        num_survivors = vector()
        node_ages = vector()
        colless_indices = vector()
        sackins_indices = vector()
        pd = vector()
        gamma = vector()

        # character simulation stats
        tip_sds = vector()
        tip_ranges = vector()
        tip_skews = vector()
        tip_kurtosis = vector()
        corr_tip_values_to_nodal_distance = vector()

        # ancestral state estimates
        sim_root = vector()
        bm_root_est = vector()
        bm_root_diff = vector()
        bm_mean_node_diff = vector()
        bm_root_ci_min = vector()
        bm_root_ci_max = vector()
        bm_root_diff_div_tip_range = vector()
        bm_root_diff_div_tip_sd = vector()
        lp_root_est = vector()
        lp_root_diff = vector()
        lp_mean_node_diff = vector()
        lp_root_diff_div_tip_range = vector()
        lp_root_diff_div_tip_sd = vector()


        # loop through each replicate
        for ( k in 1:length( results[j][[1]]$data ) ) {
           
            rep = results[j][[1]]$data[k][[1]]
            
            if (rep == "error when simulating tree") {
                
                n_error = n_error + 1

            } else if (rep == "Lineage threshold of 250 exceeded") {
              
              n_too_many = n_too_many + 1
              
            } else if (rep$error == "extinct") { 

                n_extinct = n_extinct + 1

            } else if (rep$error == "single") {
            
                n_single = n_single + 1

            } else {

                # tree shape stats
                num_survivors = c(num_survivors, length( rep$tip_states ) )
                node_ages = c(node_ages, rep$branch_times )        
                treeshape = as.treeshape( rep$trimmed_tree )
                if (length(rep$tip_states) > 2) {
                    colless_indices = c(colless_indices, colless( treeshape, norm="yule") )  
                    sackins_indices = c(sackins_indices, sackin( treeshape, norm="yule" ) )  
                }
                gamma = c(gamma, gammaStat( rep$trimmed_tree ) )
                pd = c(pd, sum(rep$trimmed_tree$edge.length))

                # character simulation stats
                tip_sds = c(tip_sds, sd( rep$tip_states ))
                tip_ranges = c(tip_ranges, max( rep$tip_states ) - min( rep$tip_states ) )
                tip_skews = c(tip_skews, skewness( rep$tip_states ) )
                tip_kurtosis = c(tip_kurtosis, kurtosis( rep$tip_states ) )
                node_dist = get_node_dist( rep$trimmed_tree )
                corr_tip_values_to_nodal_distance = c(corr_tip_values_to_nodal_distance, cor(node_dist, rep$tip_states))
 
                # ancestral state estimates
                sim_root = c(sim_root, rep$sim_anc_states[1])
                bm_root_est = c(bm_root_est, rep$est_anc_states_bm[1])
                bm_root_diff = c(bm_root_diff, rep$root_difference_bm)
                bm_mean_node_diff = c(bm_mean_node_diff, mean(rep$node_differences_bm, na.rm=TRUE))
  		    bm_root_ci_min = c(bm_root_ci_min, rep$est_data_bm$CI95[1,1])
		    bm_root_ci_max = c(bm_root_ci_max, rep$est_data_bm$CI95[1,2])
                bm_root_diff_div_tip_range = c(bm_root_diff_div_tip_range, rep$root_difference_bm / (max( rep$tip_states ) - min( rep$tip_states ) ) )
                bm_root_diff_div_tip_sd = c(bm_root_diff_div_tip_sd, rep$root_difference_bm / sd( rep$tip_states ) )
                lp_root_est = c(lp_root_est, rep$est_anc_states_lp[1])
                lp_root_diff = c(lp_root_diff, rep$root_difference_lp)
                lp_mean_node_diff = c(lp_mean_node_diff, mean(rep$node_differences_lp, na.rm=TRUE))
                lp_root_diff_div_tip_range = c(lp_root_diff_div_tip_range, rep$root_difference_lp / (max( rep$tip_states ) - min( rep$tip_states) ) )
                lp_root_diff_div_tip_sd = c(lp_root_diff_div_tip_sd, rep$root_difference_lp / sd( rep$tip_states ) )

           }
        } # end looping thru replicates

        # add data to dataframe

        # tree shape stats (means)
        d[j, "num_rep_over_250_lineages"] = n_too_many
        d[j, "num_rep_extinct"] = n_extinct
        d[j, "num_rep_error"] = n_error
        d[j, "num_rep_single"] = n_single  # record trees with a single surviving tip
        d[j, "mean_num_survivor_lineages"] = mean(num_survivors, na.rm=TRUE)
        d[j, "median_node_age"] = median(node_ages)
        d[j, "mean_tree_imbalance_colless"] = mean( colless_indices, na.rm=TRUE )
        d[j, "mean_tree_imbalance_sackin"] = mean( sackins_indices, na.rm=TRUE )
        d[j, "mean_gamma_stat"] = mean( gamma, na.rm=TRUE )
        d[j, "mean_phylo_diversity"] = mean( pd, na.rm=TRUE )

        # character simulation stats (means)
        d[j, "mean_tip_range"] = mean( tip_ranges, na.rm=TRUE )
        d[j, "mean_tip_skewness"] = mean( tip_skews, na.rm=TRUE )
        d[j, "mean_tip_kurtosis"] = mean( tip_kurtosis, na.rm=TRUE )
        d[j, "mean_corr_tip_values_to_nodal_distance"] = mean( corr_tip_values_to_nodal_distance, na.rm=TRUE )

        # ancestral state estimates (means)
        d[j, "mean_sim_root"] = mean(sim_root, na.rm=TRUE )
        d[j, "mean_bm_root_est"] = mean( bm_root_est, na.rm=TRUE )
        d[j, "mean_bm_root_diff"] = mean( bm_root_diff, na.rm=TRUE )
        d[j, "mean_bm_node_diff"] = mean( bm_mean_node_diff, na.rm=TRUE )
        d[j, "mean_bm_root_ci_min"] = mean( bm_root_ci_min, na.rm=TRUE )
        d[j, "mean_bm_root_ci_max"] = mean( bm_root_ci_max, na.rm=TRUE )
        d[j, "median_bm_root_diff_div_tip_range"] = median( bm_root_diff_div_tip_range, na.rm=TRUE )
        d[j, "median_bm_root_diff_div_tip_sd"] = median( bm_root_diff_div_tip_sd, na.rm=TRUE )
        d[j, "mean_lp_root_est"] = mean( lp_root_est, na.rm=TRUE )
        d[j, "mean_lp_root_diff"] = mean( lp_root_diff, na.rm=TRUE )
        d[j, "mean_lp_mean_node_diff"] = mean( lp_mean_node_diff, na.rm=TRUE )
        d[j, "median_lp_root_diff_div_tip_range"] = median( lp_root_diff_div_tip_range, na.rm=TRUE )
        d[j, "median_lp_root_diff_div_tip_sd"] = median( lp_root_diff_div_tip_sd, na.rm=TRUE )

    } # end looping through all simulation scenarios
    write.table(d, sep=",", row.names=FALSE, file=paste(summary_path, sim_types[i], ".csv", sep=""))
    print(paste("Wrote file: ", summary_path, sim_types[i], ".csv", sep=""))
}
print("Done.")