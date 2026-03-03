#################################################################################
# TITLE: Run_MSEbased.R
#
# PURPOSE: Run all relaxed MEM simulations and store results. 
#
# OUTPUT: MSEbased_results.csv
#
# SECTIONS: Section 0 - load source file
#           Section 1 - Continuous MEM simulations
#           Section 2 - Binary MEM simulations 
#           Section 3 - combine and export results
#
# AUTHOR: Shannon Thomas
# DATE CREATED: DEC 16, 2025
# NOTES: Takes ~8 days to run on my DELL laptop 
#        (Intel(R) Core(TM) i7-10750H CPU @ 2.60GHz, 
#         6 cores, 12 logical processors)
#################################################################################

##########################################################
############### SECTION 0: LOAD SOURCE FILE ##############
##########################################################

source("Functions_RelaxedMEM.R")
library(plyr)

library(future)
library(future.apply)

##########################################################
################ SECTION 1: CONTINUOUS MEM ###############
##########################################################


###########REMOVES set.seed() in the function while parallelizing with future package!!
#options(future.apply.debug=TRUE)
#plan(multisession, workers = 10)
plan(multisession)
inputs <- list(
  list(x = c(1,1), S = c(1,1), N = c(30,30)),
  list(x = c(1,1.1), S = c(1,1), N = c(30,30)),
  list(x = c(1,2), S = c(1,1), N = c(30,30)),
  list(x = c(1,1,1), S = c(1,1,1), N = c(30,30,30)),
  list(x = c(1,1,1.1), S = c(1,1,1), N = c(30,30,30)),
  list(x = c(1,1,2), S = c(1,1,1), N = c(30,30,30)),
  list(x = c(1,1.1,1.1), S = c(1,1,1), N = c(30,30,30)),
  list(x = c(1,2,2), S = c(1,1,1), N = c(30,30,30)),
  list(x = c(1,1.1,0.9), S = c(1,1,1), N = c(30,30,30)),
  list(x = c(1,2,0), S = c(1,1,1), N = c(30,30,30)),
  list(x = c(1,2,0.9), S = c(1,1,1), N = c(30,30,30))
)

system.time({
results_cont_cp <- future_lapply(
    inputs,
    FUN = function(vals) {findcriticalpoints(x = vals$x, S = vals$S,
                                            N = vals$N, niter = 5000,
                                            zeroweightcutoff = 0.01,
                                            stepsize = 0.1,
                                            MSEtolerance = 0.01)},
    future.seed = TRUE
)
})

plan(sequential)


inputs_cp <- inputs
for(i in 1:length(inputs)){
  inputs_cp[[i]]$cp <- results_cont_cp[[i]]
}

system.time({
results_cont <- lapply(
  inputs_cp,
  FUN = function(vals) {relaxed.MEM_sim_calc_critpts(x = vals$x, S = vals$S, 
                                                     N = vals$N, niter = 1000,
                                                     zeroweightcutoff = 0.01, 
                                                     criticalpoints = vals$cp)}
)
})

results_cont_df <- as.data.frame(do.call(rbind,results_cont))
results_cont_df$outcome_type <- "continuous"



##########################################################
################## SECTION 2: BINARY MEM #################
##########################################################

inputs_bin <- list(
  list(xvec = c(15,15,15), nvec = c(30,30,30),
       avec = c(1,1,1), bvec = c(1,1,1)),
  list(xvec = c(15,15,16), nvec = c(30,30,30),
       avec = c(1,1,1), bvec = c(1,1,1)),
  list(xvec = c(15,15,25), nvec = c(30,30,30),
       avec = c(1,1,1), bvec = c(1,1,1)),
  list(xvec = c(15,16,16), nvec = c(30,30,30),
       avec = c(1,1,1), bvec = c(1,1,1)),
  list(xvec = c(15,25,25), nvec = c(30,30,30),
       avec = c(1,1,1), bvec = c(1,1,1)),
  list(xvec = c(15,16,14), nvec = c(30,30,30),
       avec = c(1,1,1), bvec = c(1,1,1)),
  list(xvec = c(15,25,5), nvec = c(30,30,30),
       avec = c(1,1,1), bvec = c(1,1,1)),
  list(xvec = c(15,25,14), nvec = c(30,30,30),
       avec = c(1,1,1), bvec = c(1,1,1)),
  list(xvec = c(15,15), nvec = c(30,30),
       avec = c(1,1), bvec = c(1,1)),
  list(xvec = c(15,16), nvec = c(30,30),
       avec = c(1,1), bvec = c(1,1)),
  list(xvec = c(15,25), nvec = c(30,30),
       avec = c(1,1), bvec = c(1,1))
)

system.time({
results_bin_cp <- lapply(
  inputs_bin,
  FUN = function(vals) {findcriticalpoints_binary(xvec = vals$xvec, nvec = vals$nvec,
                                                  avec = vals$avec, bvec = vals$bvec,
                                                  niter = 5000, zeroweightcutoff = 0.01,
                                                  MSEtolerance = 0.005)
                       }
)
})

inputs_bin_cp <- inputs_bin
for(i in 1:length(inputs_bin)){
  inputs_bin_cp[[i]]$cp <- results_bin_cp[[i]]
}

system.time({
results_bin <- lapply(
  inputs_bin_cp,
  FUN = function(vals) {relaxed_sim_calc.betabin_critpts(x = vals$xvec, nvec = vals$nvec,
                                                         avec = vals$avec, bvec = vals$bvec,
                                                         niter = 1000, zeroweightcutoff = 0.01,
                                                         prior = "equal", criticalpoints = vals$cp)}
)
})

results_bin_df <- as.data.frame(do.call(rbind,results_bin))
results_bin_df$outcome_type <- "binary"

# alldat <- results_bin_df
# rownames(alldat) <- NULL
# 
# alldat <- alldat %>% arrange(outcome_type, numsources, diff1, diff2, model)
# alldat$diff1_c <- paste("Mean Difference SS1 =", round(alldat$diff1,2))
# alldat$diff2_c <- paste("Mean Difference SS2 =", round(alldat$diff2,2))
# 
# write.csv(alldat, file = "MSEbased_binaryresults.csv", row.names = FALSE)

##########################################################
############ SECTION 3: COMBINE AND EXPORT DATA ##########
##########################################################

alldat <- rbind.fill(results_cont_df, results_bin_df)
rownames(alldat) <- NULL

alldat <- alldat %>% arrange(outcome_type, numsources, diff1, diff2, model)
alldat$diff1_c <- paste("Mean Difference SS1 =", round(alldat$diff1,2))
alldat$diff2_c <- paste("Mean Difference SS2 =", round(alldat$diff2,2))

write.csv(alldat, file = "MSEbased_results.csv", row.names = FALSE)



