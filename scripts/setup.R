# This file loads the libraries and custom functions necessary to
# run simulations of character evolution and phylogeny under a moving
# selective regime.

library(diversitree)
library(phytools)
library(picante)
library(geiger)


# Main function that runs simulations and saves results.
simulate_trees <- function(replicates, lambda, lambda_d, mu, mu_d, char, regimes, max_time, root_state, include_extinct=TRUE) {

    # progress bar
    pb <- txtProgressBar(style = 3)
    setTxtProgressBar(pb, 0.0)

    # we will time the simulation
    start_time <- Sys.time()

    # a list to store each simulation's data in:
    simulations <- list()

    # run the simulations
    for (rep in 1:replicates) {
    
        # simulate tree and character
      sim_tree <- try(tree.quasse.regimes(list(lambda, lambda_d, mu, mu_d, char), regimes=regimes,max.t=max_time, x0=root_state, single.lineage=F, include.extinct=T))

        if (is.null(sim_tree)){simulations$data[[rep]] <- list("error when simulating tree")}    

        else if (class(sim_tree) == "try-error"){
          if (sim_tree=="Error in make.tree.quasse.regimes(pars, regimes, max.taxa, max.t, x0,  : \n  Lineage threshold of 250 exceeded\n") { # either error or null tree
            simulations$data[[rep]] <- list("Lineage threshold of 250 exceeded")
        # check if all lineages went extinct
        }  
      else {
        simulations$data[[rep]] <- list("error when simulating tree")
      }
        }
      else {

            # remember the finishing time of this simulation
            finishing_time <- max(attr(sim_tree$orig, "ages"))
            
            # NAMED vector with traits for both tips and internal nodes on the full tree
            all_traits = rep(NA,(length(sim_tree$tip.label)+length(sim_tree$node.label)));names(all_traits)=c(sim_tree$tip.label,sim_tree$node.label)
            # fill trait vector for tips
            for (i in 1:length(sim_tree$tip.label)){
              all_traits[i] = sim_tree$orig$state[which(sim_tree$orig$name==names(all_traits)[i])]
            }
            # now fill for internal nodes: node labels match idx in orig table
            all_traits["nd1"] = root_state # the true root state in the simulations
            if (length(sim_tree$node.label)>1){
              for (i in (length(sim_tree$tip.label)+2):length(all_traits)){
                all_traits[i] = sim_tree$orig$state[which(sim_tree$orig$idx2==i)]
              }
            }
            else {}
            
            # now trim the tree and get the trait values on the tips and internal node of the trimmed tree
            trimmed_tree = try(drop.extinct(sim_tree))
            tip.state=sapply(sim_tree$tip.label,FUN=function(x){strsplit(x,split="")[[1]][1]})
            if (sum(tip.state=='s')==1){simulations$data[[rep]] <- list(error="single",sim_tree=sim_tree)} # condition when only 1 lineage survived: we return the original tree in case it helps
            else {
              if (sum(tip.state=='e')==length(tip.state)){simulations$data[[rep]] <- list(error="extinct",sim_tree=sim_tree)} # condition when the tree went extinct: we return the original tree in case it helps
              else {
            branch_times=branching.times(trimmed_tree)    
            tip_states=all_traits[trimmed_tree$tip.label]
		sim_anc_states_full=all_traits[sim_tree$node.label] ## testing
		sim_anc_states=all_traits[trimmed_tree$node.label]

            # infer ancestral states using
            # Brownian motion (BM) model fitted by residual maximum likelihood
            est_data_bm <- ace(tip_states, trimmed_tree, type="continuous")
            est_root_state_bm <- as.vector(est_data_bm$ace)[1]
            est_anc_states_bm <- est_data_bm$ace # here I modified the script since the root value was binded to the vector again
            names(est_anc_states_bm)=trimmed_tree$node.label
            
            # infer ancestral states using
            # linear parsimony (LP)
            est_data_lp <- ace_lp(tip_states, trimmed_tree)
            est_root_state_lp <- est_data_lp[1]
            est_anc_states_lp <- est_data_lp
            names(est_anc_states_lp)=trimmed_tree$node.label

            # calculate contrasts between simulated and BM estimated ancestral states
            node_differences_bm <- sim_anc_states[trimmed_tree$node.label] - est_anc_states_bm[trimmed_tree$node.label]
            root_difference_bm <- sim_anc_states[1] - est_anc_states_bm[1]
            
            # calculate contrasts between simulated and LP estimated ancestral states
            node_differences_lp <- sim_anc_states[trimmed_tree$node.label] - est_anc_states_lp[trimmed_tree$node.label]
            root_difference_lp <- sim_anc_states[1] - est_anc_states_lp[1]

            # save all the data for this simulation
            simulations$data[[rep]] <- list(sim_tree=sim_tree, trimmed_tree=trimmed_tree,tip_states=tip_states, branch_times=branch_times,all_traits=all_traits,
                                          sim_anc_states_full=sim_anc_states_full, sim_anc_states=sim_anc_states, est_data_bm=est_data_bm, est_anc_states_bm=est_anc_states_bm,
                                          node_differences_bm=node_differences_bm, root_difference_bm=root_difference_bm, 
                                          est_data_lp=est_data_lp, est_anc_states_lp=est_anc_states_lp,
                                          node_differences_lp=node_differences_lp, root_difference_lp=root_difference_lp, 
                                          finishing_time=finishing_time,error="none")

              }
        }
}
        # update progress bar
        setTxtProgressBar(pb, rep/replicates)

    }

    # record the final time
    simulations$processing_time <- Sys.time() - start_time
    
    # remeber the simulations parameters
    simulations$replicates <- replicates 
    simulations$lambda <- lambda
    simulations$lambda_d <- lambda_d
    simulations$mu <- mu 
    simulations$mu_d <- mu_d 
    simulations$char <- char 
    simulations$regimes <- regimes
    simulations$max_time <- max_time
    simulations$root_state <- root_state

    writeLines("\n")

    simulations
}


# Plots various summaries of the results.
plot_simulations <- function(replicates, simulations, est_type="bm") {

    if (est_type == "bm") {
        est_anc_states = simulations[[1]]$est_anc_states_bm
        all_est_anc_states <- unlist( sapply(simulations[1:replicates], function(x){x$est_anc_states_bm}) )
    } else {
        est_anc_states = simulations[[1]]$est_anc_states_lp
        all_est_anc_states <- unlist( sapply(simulations[1:replicates], function(x){x$est_anc_states_lp}) )
    }

    par(mfrow=c(2,2))

    # plot two traitgrams on top of one another
    traitgram_given(simulations[[1]]$tip_states, simulations[[1]]$sim_anc_states, simulations[[1]]$sim_tree, simulations[[1]]$finishing_time, lab="trait", method="sim")
    traitgram_given(simulations[[1]]$tip_states, simulations[[1]]$est_anc_states_bm, simulations[[1]]$sim_tree, simulations[[1]]$finishing_time, color="red", method="est")
    traitgram_given(simulations[[1]]$tip_states, simulations[[1]]$est_anc_states_lp, simulations[[1]]$sim_tree, simulations[[1]]$finishing_time, color="orange", method="est")

    # plot rescaled branch times against difference between simulated and estimated states
    branch_times <- unlist( sapply(simulations[1:replicates], function(x){x$finishing_time - x$branch_times}) )
    sim_anc_states <- unlist( sapply(simulations[1:replicates], function(x){x$sim_anc_states}) )
    state_differences <- sim_anc_states[1:length(all_est_anc_states)] - all_est_anc_states
    plot(branch_times, state_differences, xlab="time", ylab="ancestral state differences")

    # plot lineage through time curve
    lineages <- attr(simulations[[1]]$sim_tree$orig, "lineages_thru_time")
    ages <- attr(simulations[[1]]$sim_tree$orig, "ages")
    plot(ages, lineages, type="l", xlab="time")

    # view tree unbalance
    plot(ladderize(simulations[[1]]$sim_tree), root.edge=TRUE)
    axisPhylo(backward=FALSE, root.time=round(simulations[[1]]$sim_tree$root.edge, 2))

    # plot the root differences for all simulations
    #root_differences <- sapply(simulations, function(x){x$root_difference})
    #hist(root_differences, breaks=20)

}


# Calculates ancestral states using linear parsimony as
# described in Swofford and Maddison 1987
ace_lp <- function(tip_states, tree) {

    # 1: pass down the tree and get the state sets (Farris intervals) for each internal node
    tree <- reorder(tree, "postorder")
    state_sets = list()
    for (i in 1:length(tree$edge[,1])) {
        node_id <- tree$edge[i,2]
        if (node_id <= length(tip_states))
            state_sets[[ node_id ]] <- c(tip_states[ node_id ])
        else {
            desc_id <- tree$edge[,2][ tree$edge[,1] == node_id ]
            mins <- c(min(state_sets[[ desc_id[1] ]]), min(state_sets[[ desc_id[2] ]]))
            maxs <- c(max(state_sets[[ desc_id[1] ]]), max(state_sets[[ desc_id[2] ]]))
            state_sets[[ node_id ]] <- c( min( maxs ), max( mins ) )
        }
    }
    root_id <- tree$edge[,1][length(tree$edge[,1])]
    desc_id <- tree$edge[,2][ tree$edge[,1] == root_id ]
    mins <- c(min(state_sets[[ desc_id[1] ]]), min(state_sets[[ desc_id[2] ]]))
    maxs <- c(max(state_sets[[ desc_id[1] ]]), max(state_sets[[ desc_id[2] ]]))
    state_sets[[ root_id ]] <- c( min( maxs ), max( mins ) )
    # 2: preorder traversal and for each node calculate median of ancestral node state and Farris intervals.
    anc_states <- c( median( state_sets[[ root_id ]] ) )
    for (i in length(tree$edge[,1]):1) {
        node_id <- tree$edge[i,2]
        if (node_id > root_id) {
            anc_id <- tree$edge[i,1]
            anc_state <- median( c(state_sets[[ node_id ]], anc_states[ anc_id - root_id + 1 ]) )
            anc_states[node_id - root_id + 1] <- anc_state
        }
    }
    anc_states
}


# Modified tree.quasse function to accept regimes. 
# lambda, mu, and char are each input as a list of functions where
# each function corresponds to a regime. The regime times are
# input as the vector 'regimes'.
tree.quasse.regimes <- function(pars, regimes=NA, max.taxa=Inf, max.t=Inf,
                        include.extinct=TRUE, x0=NA,
                        single.lineage=TRUE, verbose=FALSE) {
    if ( is.na(x0) )
        stop("x0 must be specified")
    else if ( length(x0) != 1 )
        stop("x0 must be of length 1")
    if ( is.na(regimes) || !is.vector(regimes) )
        stop("regimes must be a vector of times.")
    stopifnot(is.list(pars), all(sapply(pars, is.list)))
  
    info <- make.tree.quasse.regimes(pars, regimes, max.taxa, max.t, x0, single.lineage,
                             verbose)
    if ( single.lineage )
        info <- info[-1,]
    phy <- me.to.ape.quasse(info)
    if ( include.extinct || is.null(phy) )
        phy
    else {
        phy2 <- prune(phy)
        # after pruning extinct lineages, check if any survived
        if (class(phy2) != "phylo")
            NULL
        else {
            # add root stem
            finishing_time <- max(attr(phy2$orig, "ages"))
            phy2$root.edge <- finishing_time - max(branching.times(phy2))
            phy2
        }
    }
}

# Helper function to simulate quasse trees with regimes
# modified from package diversitree
make.tree.quasse.regimes <- function(pars, regimes, max.taxa=Inf, max.t=Inf, x0,
                             single.lineage=TRUE,
                             verbose=FALSE, k=500, ...) {
  lambda   <- pars[[1]]
  lambda_d <- pars[[2]]
  mu       <- pars[[3]]
  mu_d     <- pars[[4]]
  char     <- pars[[5]]
  
  if ( single.lineage ) {
    info <- data.frame(idx=1, len=1e-8, parent=0, state=x0,
                       extinct=FALSE, split=FALSE)
  } else {
    info <- data.frame(idx=1:2, len=1e-8, parent=0, state=x0,
                       extinct=FALSE, split=FALSE)
  }

  lineages <- which(!info$extinct & !info$split)
  n.taxa <- c(length(lineages))
  ages <- c(0)
  t <- 0
  t.left_regime <- regimes[1]
  t.left_total <- max.t
  i <- 1
  j <- 2
  while ( n.taxa[1] <= max.taxa && n.taxa[1] > 0 && t.left_total > 0 ) {
      verbose = FALSE
      while ( t.left_regime > 0 && t.left_total > 0 ) {
        x <- run.until.change(lineages, info, k, lambda[[i]], lambda_d[[i]], mu[[i]], mu_d[[i]], char[[i]], t.left_regime)
        lineages <- x[[1]]
#        if (length(lineages) > 250)
#            stop("Lineage threshold of 250 exceeded")
        info <- x[[2]]
        n.taxa <- c(n.taxa, length(lineages))
        t <- t + x[[4]]
        ages <- c(ages, t)
        t.left_total <- t.left_total - x[[4]]
        t.left_regime <- t.left_regime - x[[4]]
        if ( verbose )
          cat(sprintf("%s: %d [%2.3f]\n",
                      c("-", " ", "+")[sign(x[[3]])+2], n.taxa[j], t))
        j <- j + 1
      }
      i <- i + 1
      t.left_regime <- regimes[i]
  }

  if ( n.taxa[i] > max.taxa ) {
    ## Drop final speciation event.
    drop <- info[nrow(info)-1,]
    info$split[drop$parent] <- FALSE
    info$state[drop$parent] <- drop$state
    info$len[drop$parent] <- info$len[drop$parent] + drop$len
    info <- info[seq_len(nrow(info)-2),]
  }

  attr(info, "t") <- t
  attr(info, "lineages_thru_time") <- n.taxa
  attr(info, "ages") <- ages
  info
}


# Function to evolve the tree under a given regime
# modified to accept diversity dependent lambda and mu
# from package diversitree.
run.until.change <- function(lineages, info, k, lambda, lambda_d, mu, mu_d, char,
                             max.t) {
  i <- 1
  time <- 0
  n.extant <- length(lineages)
  p.change <- 1/k
  niter <- 1
  repeat {
    state <- info$state[lineages]
    lx <- lambda(state) + lambda_d(length(lineages))
    mx <- mu(state) + mu_d(length(lineages))
    r <- sum(lx + mx)
    if (r < 1.0e-3) {
        r <- 1.0e-3
    }
    dt <- 1/(r*k)
    if ( runif(1) < p.change ) {
      if ( runif(1) < sum(lx)/r ) { # speciation
        i <- sample(n.extant, 1, prob=lx)
        info <- speciate(info, lineages[i])
        lineages <- c(lineages[-i], c(-1,0) + nrow(info))
      } else {
        i <- sample(n.extant, 1, prob=mx)
        info$extinct[lineages[i]] <- TRUE
        lineages <- lineages[-i]
      }
      info$len[lineages]   <- info$len[lineages] + dt
      info$state[lineages] <- char(info$state[lineages], dt)
      time <- time + dt
      break
    }

    info$len[lineages]   <- info$len[lineages] + dt
    info$state[lineages] <- char(info$state[lineages], dt)
    niter <- niter + 1
    time <- time + dt

    if ( time > max.t )
      break
  }

  list(lineages, info, length(lineages) - n.extant, time)
}


# This function rejigs the 'lineages' structure when speciation
# happens, creating new species.
# from package diversitree
speciate <- function(info, i) {
  j <- 1:2 + nrow(info)
  info[j,"idx"] <- j
  info[j,-1] <- list(len=0, parent=i,
                     state=info$state[i],
                     extinct=FALSE, split=FALSE)
  info$split[i] <- TRUE
  info
}


# Transform the lineages structure produced by sim.tree into an ape
# phylogeny.
# from package diversitree
me.to.ape.quasse <- function(info) {
  if ( nrow(info) == 0 )
    return(NULL)
  Nnode <- sum(!info$split) - 1
  n.tips <- sum(!info$split)

  info$idx2 <- NA
  info$idx2[!info$split] <- 1:n.tips
  info$idx2[ info$split] <- order(info$idx[info$split]) + n.tips + 1

  i <- match(info$parent, info$idx)
  info$parent2 <- info$idx2[i]
  info$parent2[is.na(info$parent2)] <- n.tips + 1

  tip.label <- ifelse(subset(info, !split)$extinct,
                      sprintf("ex%d", 1:n.tips),
                      sprintf("sp%d", 1:n.tips))
  node.label <- sprintf("nd%d", 1:Nnode)

  info$name <- NA
  info$name[!info$split] <- tip.label

  tip.state <- info$state[match(1:n.tips, info$idx2)]
  names(tip.state) <- tip.label
  phy <- reorder(structure(list(edge=cbind(info$parent2, info$idx2),
                               Nnode=Nnode,
                               tip.label=tip.label,
                               tip.state=tip.state,
                               node.label=node.label,
                               edge.length=info$len,
                               orig=info),
                          class="phylo"))

  phy$edge.state <- info$state[match(phy$edge[,2], info$idx2)]
  phy
}

# Modified traitgram
traitgram_given_Flo=function (x,col_sim="black",col_bm="red",col_lp="blue",ymax=1,yoffset=0.1){ 
  # x is a single element of a sim$data object, e.g. sim$data[[1]]
  # ymax and yoffset can eventually be adjusted so that tip labels can be seen
  # colors for the simulated trait values and the two estimations can be customized
  tree=x$trimmed_tree
  node_labels_full=c(tree$tip.label,tree$node.label)
  branching_times_full=c(rep(0,length(tree$tip.label)),-x$branch_times) ; names(branching_times_full)[1:length(tree$tip.label)]=tree$tip.label
  tip_states=x$tip_states
  sim=c(tip_states,x$sim_anc_states)
  bm=c(tip_states,x$est_anc_states_bm)
  lp=c(tip_states,x$est_anc_states_lp)
  plot(x=sim[tree$edge[1,1]],y=-(x$finishing_time),xlim=c(min(c(sim,bm,lp)),max(c(sim,bm,lp))),ylim=c(-max(x$branch_times),ymax),pch=19,col=col_sim,cex=0.01,xlab='Trait',ylab='Time')
  for (b in 1:dim(tree$edge)[1]){
    segments(x0=sim[tree$edge[b,1]],x1=sim[tree$edge[b,2]],y0=branching_times_full[tree$edge[b,1]],y1=branching_times_full[tree$edge[b,2]],col=col_sim)
    segments(x0=bm[tree$edge[b,1]],x1=bm[tree$edge[b,2]],y0=branching_times_full[tree$edge[b,1]],y1=branching_times_full[tree$edge[b,2]],col=col_bm,lty='dashed')
    segments(x0=lp[tree$edge[b,1]],x1=lp[tree$edge[b,2]],y0=branching_times_full[tree$edge[b,1]],y1=branching_times_full[tree$edge[b,2]],col=col_lp,lty='dashed')
  }
  text(x=tip_states,y=yoffset,labels=tree$tip.label,srt=90,pos=3)
}
