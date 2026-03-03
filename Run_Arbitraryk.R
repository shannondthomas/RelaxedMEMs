#################################################################################
# TITLE: Run_Arbitraryk.R
#
# PURPOSE: Run all relaxed MEM simulations and store results. 
#
# OUTPUT: arbitraryk_results.csv
#
# SECTIONS: Section 0 - load source file
#           Section 1 - Continuous MEM simulations
#           Section 2 - Binary MEM simulations 
#           Section 3 - combine and export results
#
# AUTHOR: Shannon Thomas
# DATE CREATED: DEC 1, 2025
#################################################################################

##########################################################
############### SECTION 0: LOAD SOURCE FILE ##############
##########################################################

source("Functions_RelaxedMEM.R")
library(plyr)


##########################################################
################ SECTION 1: CONTINUOUS MEM ###############
##########################################################

## THREE SOURCE
system.time({
  #ALL EQUAL
  allequal_cont <- relaxed.MEM_sim_calc(x = c(1,1,1), S = c(1,1,1), 
                                   N = c(30,30,30), niter = 1000, k = 0.9,
                                   zeroweightcutoff = 0.01,func_seed=TRUE)
})
system.time({
  #ONE SMALL DIFFERENCE
  onesmalldiff_cont <- relaxed.MEM_sim_calc(x = c(1,1,1.1), S = c(1,1,1), 
                                       N = c(30,30,30), niter = 1000, k = 0.9,
                                       zeroweightcutoff = 0.01,func_seed=TRUE)
})
system.time({
  #ONE LARGE DIFFERENCE
  onelargediff_cont <- relaxed.MEM_sim_calc(x = c(1,1,2), S = c(1,1,1), 
                                       N = c(30,30,30), niter = 1000, k = 0.9,
                                       zeroweightcutoff = 0.01,func_seed=TRUE)
})
system.time({
  #PRIMARY UNEQUAL - SMALL DIFFERENCE
  primaryunequalsmall_cont <- relaxed.MEM_sim_calc(x = c(1,1.1,1.1), S = c(1,1,1), 
                                              N = c(30,30,30), niter = 1000, k = 0.9,
                                              zeroweightcutoff = 0.01,func_seed=TRUE)
})
system.time({
  #PRIMARY UNEQUAL - LARGE DIFFERENCE
  primaryunequallarge_cont <- relaxed.MEM_sim_calc(x = c(1,2,2), S = c(1,1,1), 
                                              N = c(30,30,30), niter = 1000, k = 0.9,
                                              zeroweightcutoff = 0.01,func_seed=TRUE)
})
system.time({
  #PRIMARY AS MIDPOINT - SMALL DIFFERENCE
  primarymidpointsmall_cont <- relaxed.MEM_sim_calc(x = c(1,1.1,0.9), S = c(1,1,1), 
                                               N = c(30,30,30), niter = 1000, k = 0.9,
                                               zeroweightcutoff = 0.01,func_seed=TRUE)
})
system.time({
  #PRIMARY AS MIDPOINT - LARGE DIFFERENCE
  primarymidpointlarge_cont <- relaxed.MEM_sim_calc(x = c(1,2,0), S = c(1,1,1), 
                                               N = c(30,30,30), niter = 1000, k = 0.9,
                                               zeroweightcutoff = 0.01,func_seed=TRUE)
})
system.time({
  #PRIMARY BETWEEN 
  primarybetween_cont <- relaxed.MEM_sim_calc(x = c(1,2,0.9), S = c(1,1,1), 
                                         N = c(30,30,30), niter = 1000, k = 0.9,
                                         zeroweightcutoff = 0.01,func_seed=TRUE)
})

## TWO SOURCE
system.time({
  #ALL EQUAL
  allequal_2source_cont <- relaxed.MEM_sim_calc(x = c(1,1), S = c(1,1),
                                           N = c(30,30), niter = 1000, k = 0.9,
                                           zeroweightcutoff = 0.01,func_seed=TRUE)
})
system.time({
  #ONE SMALL DIFFERENCE
  onesmalldiff_2source_cont <- relaxed.MEM_sim_calc(x = c(1,1.1), S = c(1,1),
                                               N = c(30,30), niter = 1000, k = 0.9,
                                               zeroweightcutoff = 0.01,func_seed=TRUE)
})
system.time({
  #ONE LARGE DIFFERENCE
  onelargediff_2source_cont <- relaxed.MEM_sim_calc(x = c(1,2), S = c(1,1),
                                               N = c(30,30), niter = 1000, k = 0.9,
                                               zeroweightcutoff = 0.01,func_seed=TRUE)
})

## COMBINE CONTINUOUS DATA
alldat_cont <- rbind(allequal_cont, onesmalldiff_cont, onelargediff_cont,
                     primarybetween_cont, primarymidpointsmall_cont, primarymidpointlarge_cont,
                     primaryunequalsmall_cont, primaryunequallarge_cont,
                     allequal_2source_cont, onesmalldiff_2source_cont, onelargediff_2source_cont) %>% relocate(model)

alldat_cont$outcome_type <- "continuous"


##########################################################
################## SECTION 2: BINARY MEM #################
##########################################################

## THREE SOURCE
system.time({
  #ALL EQUAL
  allequal_bin <- relaxed_sim_calc.betabin(xvec = c(15,15,15), nvec = c(30,30,30), 
                                           avec = c(1,1,1), bvec = c(1,1,1),
                                           prior = "equal", constraint = 1,
                                           k = 0.9, zeroweightcutoff = 0.01,
                                           niter = 1000)
})
system.time({
  #ONE SMALL DIFFERENCE
  onesmalldiff_bin <- relaxed_sim_calc.betabin(xvec = c(15,15,16), nvec = c(30,30,30), 
                                               avec = c(1,1,1), bvec = c(1,1,1),
                                               prior = "equal", constraint = 1,
                                               k = 0.9, zeroweightcutoff = 0.01,
                                               niter = 1000)
})
system.time({
  #ONE LARGE DIFFERENCE
  onelargediff_bin <- relaxed_sim_calc.betabin(xvec = c(15,15,25), nvec = c(30,30,30), 
                                               avec = c(1,1,1), bvec = c(1,1,1),
                                               prior = "equal", constraint = 1,
                                               k = 0.9, zeroweightcutoff = 0.01,
                                               niter = 1000)
})
system.time({
  #PRIMARY UNEQUAL - SMALL DIFFERENCE
  primaryunequalsmall_bin <- relaxed_sim_calc.betabin(xvec = c(15,16,16), nvec = c(30,30,30), 
                                                      avec = c(1,1,1), bvec = c(1,1,1),
                                                      prior = "equal", constraint = 1,
                                                      k = 0.9, zeroweightcutoff = 0.01,
                                                      niter = 1000)
})
system.time({
  #PRIMARY UNEQUAL - LARGE DIFFERENCE
  primaryunequallarge_bin <- relaxed_sim_calc.betabin(xvec = c(15,25,25), nvec = c(30,30,30), 
                                                      avec = c(1,1,1), bvec = c(1,1,1),
                                                      prior = "equal", constraint = 1,
                                                      k = 0.9, zeroweightcutoff = 0.01,
                                                      niter = 1000)
})
system.time({
  #PRIMARY AS MIDPOINT - SMALL DIFFERENCE
  primarymidpointsmall_bin <- relaxed_sim_calc.betabin(xvec = c(15,16,14), nvec = c(30,30,30), 
                                                       avec = c(1,1,1), bvec = c(1,1,1),
                                                       prior = "equal", constraint = 1,
                                                       k = 0.9, zeroweightcutoff = 0.01,
                                                       niter = 1000)
})
system.time({
  #PRIMARY AS MIDPOINT - LARGE DIFFERENCE
  primarymidpointlarge_bin <- relaxed_sim_calc.betabin(xvec = c(15,25,5), nvec = c(30,30,30), 
                                                       avec = c(1,1,1), bvec = c(1,1,1),
                                                       prior = "equal", constraint = 1,
                                                       k = 0.9, zeroweightcutoff = 0.01,
                                                       niter = 1000)
})
system.time({
  #PRIMARY BETWEEN 
  primarybetween_bin <- relaxed_sim_calc.betabin(xvec = c(15,25,14), nvec = c(30,30,30), 
                                                 avec = c(1,1,1), bvec = c(1,1,1),
                                                 prior = "equal", constraint = 1,
                                                 k = 0.9, zeroweightcutoff = 0.01,
                                                 niter = 1000)
})


## TWO SOURCE
system.time({
  #ALL EQUAL
  allequal_2source_bin <- relaxed_sim_calc.betabin(xvec = c(15,15), nvec = c(30,30), 
                                                   avec = c(1,1), bvec = c(1,1),
                                                   prior = "equal", constraint = 1,
                                                   k = 0.9, zeroweightcutoff = 0.01,
                                                   niter = 1000)
})
system.time({
  #ONE SMALL DIFFERENCE
  onesmalldiff_2source_bin <- relaxed_sim_calc.betabin(xvec = c(15,16), nvec = c(30,30), 
                                                       avec = c(1,1), bvec = c(1,1),
                                                       prior = "equal", constraint = 1,
                                                       k = 0.9, zeroweightcutoff = 0.01,
                                                       niter = 1000)
})
system.time({
  #ONE LARGE DIFFERENCE
  onelargediff_2source_bin <- relaxed_sim_calc.betabin(xvec = c(15,25), nvec = c(30,30), 
                                                       avec = c(1,1), bvec = c(1,1),
                                                       prior = "equal", constraint = 1,
                                                       k = 0.9, zeroweightcutoff = 0.01,
                                                       niter = 1000)
})

alldat_bin <- rbind(allequal_bin, onesmalldiff_bin, onelargediff_bin,
                primarybetween_bin, primarymidpointsmall_bin, primarymidpointlarge_bin,
                primaryunequalsmall_bin, primaryunequallarge_bin,
                allequal_2source_bin, onesmalldiff_2source_bin, onelargediff_2source_bin) %>% relocate(model)

alldat_bin$outcome_type <- "binary"

##########################################################
############ SECTION 3: COMBINE AND EXPORT DATA ##########
##########################################################

alldat <- rbind.fill(alldat_cont, alldat_bin)
rownames(alldat) <- NULL

alldat <- alldat %>% arrange(outcome_type, numsources, diff1, diff2, k, model)
alldat$diff1_c <- paste("Mean Difference SS1 =", round(alldat$diff1,2))
alldat$diff2_c <- paste("Mean Difference SS2 =", round(alldat$diff2,2))

write.csv(alldat, file = "arbitraryk_results.csv", row.names = FALSE)



