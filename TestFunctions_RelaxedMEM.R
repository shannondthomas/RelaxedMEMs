#################################################################################
# TITLE: Functions2_RelaxedMEM.R
#
# PURPOSE: Store functions to calculate weights and summary performance measures 
#          for both continuous and binary relaxed MEMs based on MSE instead of 
#          arbitrary k. 
#
# OUTPUT: NA
#
# SECTIONS: Section 0 - load dependencies
#           Section 1 - Continuous MEMs
#              boa.hpd() - Calculate HPD Intervals
#              general.calc.weights_MEM() - calculate continuous MEM weights for 
#                                           any number of sources
#              relaxed.MEM_sim_cal() - run all continuous MEMs (relaxed and 
#                                      original) and summarize performance
#           Section 2 - Binary MEMs 
#              BinaryESS2 - calculate ESS for binary data
#              calc.MEM.betabin - calculate binary MEM weights for any number of
#                                 sources
#              relaxed_summary_calc.betabin - run all binary MEMs (relaxed and
#                                             original) and summarize 
#                                             performance
#
# AUTHOR: Shannon Thomas
# DATE CREATED: DEC 15, 2025
#################################################################################

##########################################################
############## SECTION 0: LOAD DEPENDENCIES ##############
##########################################################

source("Functions_RelaxedMEM.R")
library(moments)

##########################################################
##### SECTION 1: CONTINUOUS MEM FUNCTIONS FOR TESTING ####
##########################################################

###Run all continuous MEMs (relaxed and original) and summarize performance
TEST_relaxed.MEM_sim_calc = function(x,S,N,prior = "pi_e",niter,k,zeroweightcutoff,RAND=TRUE,func_seed=FALSE){
  ###x: vector mean values for primary and then supplemental cohorts
  #S: vector of standard deviations for primary and then supplemental cohorts
  #N: sample size for primary and then supplemental cohorts
  #prior: prior to use for calculation (must be "pi_e" or a vector of length 2^H)
  #niter: number of iterations
  
  mu = x[1]
  est.mu <- est.var <- est.mse <- est.bias <- est.esss <-  est.xbar <- cov.ind <- hpd.width <- data.frame(MEM = rep(NA,niter),
                                                                                                          relaxedMEM_max = rep(NA, niter),
                                                                                                          relaxedMEM_sqrtmax = rep(NA, niter),
                                                                                                          relaxedMEM_allornothing = rep(NA, niter),
                                                                                                          relaxedMEM_twosource = rep(NA, niter),
                                                                                                          relaxedMEM_stretchsqrt = rep(NA, niter),
                                                                                                          relaxedMEM_stretch1 = rep(NA, niter)) #initialize df to store results
  
  simx <- rep(NA, niter)
  
  var.out <- rep(NA,7)
  
  numsources <- length(x)
  H <- numsources - 1                              #number of supplemental sources
  s_exchindicator <- expand.grid(rep(list(0:1),H)) #data frame with indicators for supp source exchangeability for each model
  nummodels <- 2^H                                 #number of models = nrow(s_exchindicator)
  
  s_exchindicator_wprimary <- cbind(rep(1,nummodels), s_exchindicator)
  
  #############RELAXED MEM ADDITION
  ###Calculate theoretical max for each weight (or sums of weight)
  ##Gives a data frame, max_w, containing the maximum weight ($objective)
  ##and corresponding x value ($maximum) for full exchangeability (first row)
  ##and each individual source (all remaining rows)
  max_w_pexch <- as.data.frame(optimize(function(z) general.calc.weights_MEM(xvec = c(z,x[-1]), svec = S, nvec = N)$weight[nummodels], 
                                        lower = 0.5*min(x), upper = 2*max(x), maximum = TRUE, tol = 0.00001))
  
  max_w_source <- as.data.frame(t(mapply(c(1:H),
                                         FUN = function(j) optimize(function(z) sum(s_exchindicator[,j]*general.calc.weights_MEM(xvec = c(z,x[-1]), svec = S, nvec = N)$weight), 
                                                                    lower = 0.5*min(x), upper = 2*max(x), maximum = TRUE, tol = 0.00001))))
  max_w <- rbind(max_w_pexch, max_w_source)
  
  max_w$objective <- ifelse(max_w$objective < zeroweightcutoff, 0, max_w$objective)
  
  for(i in 1:niter){
    if(func_seed == TRUE){
      set.seed(515+i)
    }
    if((RAND == TRUE) & (niter != 1)){
      xbar = rnorm(1,mean=mu,sd=S[1]/sqrt(N[1])) 
    }else{
      xbar <- x[1]
    }
    
    simx[i] <- xbar
    
    V <- S^2/N
    twosource_for_threesource_ind <- 0
    
    ##Calculate weights and derivative terms needed for variance calculation
    w_and_der <- general.calc.weights_MEM(xvec=c(xbar,x[-1]),svec=S,nvec=N,prior=prior)
    w <- w_and_der$weight
    b <- w_and_der$expderiv
    
    # Remove all-zero row (already covered by M1 = xbar) and not needed for weight replacement
    s_exchindicator <- s_exchindicator[rowSums(s_exchindicator) > 0, , drop = FALSE]
    
    #############RELAXED MEM ADDITION
    ### Check weights and redefine if within threshold
    threshold_ind <- (w[nummodels] > k*as.numeric(max_w$objective[1]))
    threshold_ind <- c(threshold_ind,
                       as.data.frame(mapply(c(1:H),
                                            FUN = function(j) sum(s_exchindicator[,j]*w[2:nummodels]))) > k*as.numeric(max_w$objective[2:H]))
    threshold_ind <- ifelse(max_w$objective == 0, FALSE, threshold_ind)
    
    
    ## Step 1: Check P(EXCH) (full exchangeability model)
    if((threshold_ind[1] == TRUE) ||
       ((numsources == 3) & (threshold_ind[2] == TRUE) & (threshold_ind[3] == TRUE))
    ){
      #get weights at maximum x value
      #w_and_der2 <-  general.calc.weights_MEM(xvec = c(as.numeric(max_w$maximum[nummodels-1]),x[-1]), svec = S, nvec = N, prior=prior)
      w_and_der2 <-  general.calc.weights_MEM(xvec = c(as.numeric(max_w$maximum[1]),x[-1]), svec = S, nvec = N, prior=prior)
      w2 <- w_and_der2$weight
      b2 <- w_and_der2$expderiv
      
      #calculate weights for P(EXCH) = sqrt(P(EXCH)_{max})
      w3 <- w2
      w3[nummodels] <- sqrt(w2[nummodels])
      w3[1:(nummodels-1)] <- (1-w3[nummodels])*(w2[1:(nummodels-1)]/sum(w2[1:(nummodels-1)]))
      
      #define all or nothing borrowing weights
      w4 <- c(rep(0,nummodels - 1),1)
      
      #two source method uses OG three source results in this setting
      w5 <- w
      
      #define stretch weight 1 (sqrt(w))
      w6 <- w
      w6[nummodels] <- sqrt(w[nummodels])
      w6[1:(nummodels - 1)] <- (1-w6[nummodels])*(w[1:(nummodels-1)]/sum(w[1:(nummodels-1)]))
      
      #define stretch weight 2 (maximum is 1)
      w7 <- w
      w7[nummodels] <- (1/w2[nummodels])*w[nummodels]
      w7[1:(nummodels - 1)] <- (1-w7[nummodels])*(w[1:(nummodels-1)]/sum(w[1:(nummodels-1)]))
      
      allweights <- rbind(w,w2,w3,w4,w5,w6,w7)
      
    } else if((numsources == 3) & 
              ((threshold_ind[2] == TRUE) | (threshold_ind[3] == TRUE))
    ){ ## Step 2: Check P(EXCH)_s for external sources (if numsources > 2)
      exch_s_ind <- which(threshold_ind == TRUE)
      
      xval <- as.numeric(max_w$maximum[exch_s_ind])
      
      #get weights at maximum x value
      w_and_der2 <-  general.calc.weights_MEM(xvec = c(xval,x[-1]), svec = S, nvec = N, prior=prior)
      w2 <- w_and_der2$weight
      b2 <- w_and_der2$expderiv
      
      #calculate weights for P(EXCH) = sqrt(P(EXCH)_{max})
      w3 <- w2
      w3[exch_s_ind] <- sqrt(w2[exch_s_ind])
      w3[c(1:4)[-exch_s_ind]] <- (1-w3[exch_s_ind])*(w2[-exch_s_ind]/sum(w2[-exch_s_ind]))
      
      #define all or nothing borrowing weights
      w4 <- rep(0,nummodels)
      w4[exch_s_ind] <- 1
      
      #get two-source model weights 
      twosource_for_threesource_ind <- 1
      w_and_der5 <-  general.calc.weights_MEM(xvec = c(xbar,x[exch_s_ind]), svec = S[c(1,exch_s_ind)], nvec = N[c(1,exch_s_ind)], prior=prior)
      w5 <- rep(NA,4)
      w5[1] <- w_and_der5$weight[1]
      w5[exch_s_ind] <- w_and_der5$weight[2]
      b5 <- w_and_der5$expderiv
      
      #define stretch weight 1 (sqrt(w))
      w6 <- w
      w6[exch_s_ind] <- sqrt(w[exch_s_ind])
      w6[c(1:4)[-exch_s_ind]] <- (1-w6[exch_s_ind])*(w[-exch_s_ind]/sum(w[-exch_s_ind]))
      
      #define stretch weight 2 (maximum is 1)
      w7 <- w
      totalweightforsuppsource <- w2[exch_s_ind] + w2[nummodels]
      propweightforONLYsuppsource <- w2[exch_s_ind]/(w2[nummodels]+w2[exch_s_ind])
      w7[exch_s_ind] <- (1/totalweightforsuppsource)*propweightforONLYsuppsource*(w[exch_s_ind] + w[nummodels])
      w7[c(1:4)[-exch_s_ind]] <- (1-w7[exch_s_ind])*(w[-exch_s_ind]/sum(w[-exch_s_ind]))
      
      allweights <- rbind(w,w2,w3,w4,w5,w6,w7)
      
    } else{
      allweights <- rbind(w,w,w,w,w,w,w) 
    }
    
    #getting an error with very small negative weights: ex: -7e-12, so reset those to 0
      # if(sum(allweights<0)>0){
      #   print(paste("i=",i,"\n x = ",x,"\n S=",S,"\n N=",N,"\n allweights=",allweights))
      # }
    allweights[abs(allweights)<1e-10] <- 0
    
    for(l in 1:nrow(allweights)){
      ############LOOP THROUGH ALL WEIGHTS TO GET PERFORMANCE FOR EACH WEIGHT ADJUSTMENT
      
      V1 <- V[1]
      v_rest <- V[-1]
      x_rest <- x[-1]
      
      s_exchindicator2 <- expand.grid(rep(list(0:1),H))
      s_exchindicator2 <- s_exchindicator2[rowSums(s_exchindicator2) > 0, , drop = FALSE]
      currentweight <- allweights[l,]
      currentb <- b
      currentnummodels <- nummodels
      
      #FOR THE TWO-SOURCE FOR THREE SOURCE METHOD, ONLY USE DATA FROM THE TWO SOURCES
      if((l == 5) & (twosource_for_threesource_ind == 1)){
        #edit inputs to only include two sources
        v_rest <- V[exch_s_ind]
        x_rest <- x[exch_s_ind]
        s_exchindicator2 <- data.frame(Var = s_exchindicator2[1,1])
        currentweight <- allweights[l,c(1,exch_s_ind)]
        currentb <- b5
        currentnummodels <- 2
      }
      
      #############BACK TO STANDARD MEM CODE
      
      ###Posterior components
      ##Means for each model
      
      # Compute all products for each combination
      v_prod <- apply(s_exchindicator2, 1, function(ind) prod(v_rest ^ ind))
      
      numerators <- mapply(function(ind_row, vprod) {
        ind <- as.numeric(ind_row)
        idxs <- which(ind == 1)
        num <- vprod * xbar
        lo_prod <- sapply(idxs, function(j) {
          other_idxs <- setdiff(idxs, j)
          if (length(other_idxs) == 0) 1 else prod(v_rest[other_idxs])
        })
        num + sum(V1 * lo_prod * x_rest[idxs])
      }, split(s_exchindicator2, seq_len(nrow(s_exchindicator2))), v_prod)
      
      denominators <- mapply(function(ind_row, vprod) {
        ind <- as.numeric(ind_row)
        idxs <- which(ind == 1)
        denom <- vprod
        lo_prod <- sapply(idxs, function(j) {
          other_idxs <- setdiff(idxs, j)
          if (length(other_idxs) == 0) 1 else prod(v_rest[other_idxs])
        })
        denom + sum(V1 * lo_prod)
      }, split(s_exchindicator2, seq_len(nrow(s_exchindicator2))), v_prod)
      
      #Means for each model
      M <- c(xbar, numerators / denominators)
      #Derivative of the posterior mean for delta method variance calculation wrt xbar
      dM <- c(1, v_prod/denominators)
      
      
      ##Variances for each model
      # Compute sums of 1/v0j for each row
      inv_sums <- apply(s_exchindicator2, 1, function(ind) sum(ind / v_rest))
      
      # Add 1/v and take reciprocal
      Vmod_rest <- 1 / (1 / V1 + inv_sums)
      Vmod <- c(V1, Vmod_rest)  
      
      ##Effective historical sample size for each model
      E <- S[1]^2/Vmod - N[1]
      
      ###Calculate MSE
      ##Variance component
      #Weight denominator needed for delta method g function calculation (weight for each model):
      ##Weight is w1 = 1/a1, calculations as a1 for ease of inclusion with quotient rule for derivatives
      #note that 1/a1 equals the weight for model 1, and 1/a2 for model 2, etc.
      a <- sum(currentweight)/currentweight
      
      #Derivative of weight portion for delta method variance calculation (i.e., c1=the derivative of a1 wrt xbar)
      # Create a matrix where each row is b - b[i]
      B <- matrix(rep(currentb, each = currentnummodels), nrow = currentnummodels)
      diff_mat <- B - t(B)  # diff_mat[i, j] = b[j] - b[i]
      # Multiply each column by corresponding weight w[j]
      weighted_diff <- diff_mat * matrix(rep(currentweight, each = currentnummodels), nrow = currentnummodels)
      # Divide row sums by currentweight[i]
      C <- rowSums(weighted_diff) / currentweight
      
      
      #Derivative of g for delta method variance calculation
      g = (a*dM - M*C)/a^2
      gmat = g%*%t(g)
      
      mu.est <- sum(currentweight*M)
      
      var.comp = rep(1,currentnummodels)%*%gmat%*%rep(1,currentnummodels) * V[1] #accounting for covariance, variance of models
      
      bias.orig = (mu.est-mu)
      mse = var.comp + bias.orig^2
      
      esss.orig = sum(E*currentweight) 
      
      post = t(rmultinom(1,30000,currentweight))
      #post = t(rmultinom(1,5000,currentweight))
      
      post.mod <- mapply(function(n, mu, sigma2) {
        rnorm(n, mean = mu, sd = sqrt(sigma2))
      }, n = post[1, ], mu = M, sigma2 = Vmod)
      
      post.mu = unlist(post.mod)
      hpd.bma_rep = boa.hpd(x=post.mu,alpha=.05)
      
      include.mu = (hpd.bma_rep[1] < mu) & (mu < hpd.bma_rep[2])
      diff.hpd = hpd.bma_rep[2] - hpd.bma_rep[1]
      
      var.out[l] <- var.comp #store theoretical variance for RAND = FALSE and niter = 1
      
      est.mu[i,l] <- mu.est      #store simulation results for mean
      est.bias[i,l] <- bias.orig #store simulation results for bias
      est.esss[i,l] <- esss.orig #store simulation results for esss
      cov.ind[i,l] <- include.mu #store simulation results for coverage
      hpd.width[i,l] <- diff.hpd #store simulation results for HPD interval
      
      
    }
  }
  
  
  return(list(simx = simx,
              est.mu = est.mu, 
              est.bias = est.bias,
              est.esss = est.esss,
              cov.ind = cov.ind,
              hpd.width = hpd.width
              ))
  
  # est.var = sapply(est.mu, var) #estimate variance of estimated mu-hat from simulation
  # if((RAND==FALSE)&(niter==1)){
  #   est.var <- var.out
  # }
  # est.mean = sapply(est.mu, mean) #estimate mu-hat
  # est.bias = est.mean - mu #bias
  # mse.part = est.bias^2 + est.var #MSE
  # mse.sum.part = sapply(est.mu, FUN = function(z) (1/niter)*sum((z - mu)^2)) #MSE sum equation
  # mae.part =  sapply(est.mu, FUN = function(z) (1/niter)*sum(abs(z - mu))) #MAE
  # est.esss = sapply(est.esss, median) #ESSS
  # abserr.part = abs(est.bias) #absolute error of bias for summary stat
  # cov.part = sapply(cov.ind, mean) #coverage
  # hpd.part = sapply(hpd.width, mean) #HPD width
  # 
  # mat = cbind(mu,est.mean,est.var,est.bias,mse.part,mse.sum.part,mae.part,est.esss,abserr.part,cov.part,hpd.part)
  # 
  # df <- as.data.frame(mat)
  # df$num_iterations <- niter
  # df$model <- rownames(df)
  # df$numsources <- numsources
  # if(numsources == 3){
  #   df$diff1 <- x[2] - x[1]
  #   df$diff2 <- x[3] - x[1]
  # }else if(numsources == 2){
  #   df$diff1 <- x[2] - x[1]
  #   df$diff2 <- NA
  # }
  # df$k <- k
  # return(df)
}



##########################################################
####### SECTION 3: BINARY MEM FUNCTIONS FOR TESTING ######
##########################################################


###Run all binary MEMs (relaxed and original) and summarize performance
TEST_relaxed_sim_calc.betabin <- function(xvec,nvec,avec,bvec,prior,constraint=1, niter, zeroweightcutoff, k, RAND=TRUE){
  #p: probability of the outcome for the primary cohort
  #xvec: vector with counts of those with event (primary source first, then supplemental)
  #nvec: vector of total number of subjects (primary source first, supplemental afterwards)
  #avec: vector of alpha parameters for beta priors for each source
  #bvec: vector of beta parameters for beta priors for each source
  #prior: prior to use on MEM source inclusion probability: equal
  #constraint: value to limit maximum value of empirical Bayes prior to
  
  
  #############RELAXED MEM ADDITION
  ###Calculate theoretical max for each weight (or sums of weight)
  ##Gives a data frame, max_w, containing the maximum weight ($objective)
  ##and corresponding x value ($maximum) for full exchangeability (first row)
  ##and each individual source (all remaining rows)
  
  #initialize df to store results
  est.esss <- est.mu <- est.bias <- var.est <- data.frame(MEM = rep(NA,niter),
                                                          relaxedMEM_max = rep(NA, niter),
                                                          relaxedMEM_sqrtmax = rep(NA, niter),
                                                          relaxedMEM_allornothing = rep(NA, niter),
                                                          relaxedMEM_twosource = rep(NA, niter),
                                                          relaxedMEM_stretchsqrt = rep(NA, niter),
                                                          relaxedMEM_stretch1 = rep(NA, niter)) #initialize df to store results
  
  simp <- rep(NA, niter)
  
  truep <- xvec[1]/nvec[1]
  numsources <- length(xvec)
  H <- numsources - 1                              #number of supplemental sources
  s_exchindicator <- expand.grid(rep(list(0:1),H)) #data frame with indicators for supp source exchangeability for each model
  nummodels <- 2^H  
  s_exchindicator_wprimary <- cbind(rep(1,nummodels), s_exchindicator)
  
  max_w_pexch <- as.data.frame(optimize(function(z) calc.MEM.betabin(xvec = c(z,xvec[-1]), nvec = nvec, 
                                                                     avec = avec, bvec = bvec, 
                                                                     prior = prior)$q[nummodels], 
                                        lower = 0, upper = 2*max(xvec), maximum = TRUE, tol = 0.00001))
  
  max_w_source <- as.data.frame(t(mapply(c(1:H),
                                         FUN = function(j) optimize(function(z) sum(s_exchindicator[,j]*calc.MEM.betabin(xvec = c(z,xvec[-1]), 
                                                                                                                         nvec = nvec,
                                                                                                                         avec = avec, bvec = bvec, 
                                                                                                                         prior = prior)$q), 
                                                                    lower = 0, upper = 2*max(xvec), maximum = TRUE, tol = 0.00001))))
  max_w <- rbind(max_w_pexch, max_w_source)
  
  max_w$objective <- ifelse(max_w$objective < zeroweightcutoff, 0, max_w$objective)
  
  for(i in 1:niter){
    
    #simulate number of events in primary source
    set.seed(10+i)
    if((RAND==TRUE)&(niter!=1)){
      simx <- rbinom(n = 1, size = nvec[1], prob = truep)
    }else{
      simx <- xvec[1]
    }
    
    simp[i] <- simx
    
    #Calculate standard MEM weights
    res <- calc.MEM.betabin(xvec=c(simx, xvec[-1]), nvec=nvec, avec=avec, bvec=bvec, prior=prior, constraint=constraint)
    w <- c(res$q)
    
    #############RELAXED MEM ADDITION
    ### Check weights and redefine if within threshold
    s_exchindicator <- s_exchindicator[rowSums(s_exchindicator) > 0, , drop = FALSE]
    
    threshold_ind <- (w[nummodels] > k*as.numeric(max_w$objective[1]))
    threshold_ind <- c(threshold_ind,
                       as.data.frame(mapply(c(1:H),
                                            FUN = function(j) sum(s_exchindicator[,j]*w[2:nummodels]))) > k*as.numeric(max_w$objective[2:H]))
    threshold_ind <- ifelse(max_w$objective == 0, FALSE, threshold_ind)
    
    twosource_for_threesource_ind <- 0
    
    ## Step 1: Check P(EXCH) (full exchangeability model)
    if((threshold_ind[1] == TRUE) ||
       ((numsources == 3) & (threshold_ind[2] == TRUE) & (threshold_ind[3] == TRUE))
    ){
      #get weights at maximum x value
      #w_and_der2 <-  general.calc.weights_MEM(xvec = c(as.numeric(max_w$maximum[nummodels-1]),x[-1]), svec = S, nvec = N, prior=prior)
      w2_all <- calc.MEM.betabin(xvec = c(as.numeric(max_w$maximum[1]),xvec[-1]), nvec = nvec,
                                 avec = avec, bvec = bvec, prior=prior)
      w2 <- c(w2_all$q)
      
      #calculate weights for P(EXCH) = sqrt(P(EXCH)_{max})
      w3 <- w2
      w3[nummodels] <- sqrt(w2[nummodels])
      w3[1:(nummodels-1)] <- (1-w3[nummodels])*(w2[1:(nummodels-1)]/sum(w2[1:(nummodels-1)]))
      
      #define all or nothing borrowing weights
      w4 <- c(rep(0,nummodels - 1),1)
      
      #two source method uses OG three source results in this setting
      w5 <- w
      
      #define stretch weight 1 (sqrt(w))
      w6 <- w
      w6[nummodels] <- sqrt(w[nummodels])
      w6[1:(nummodels - 1)] <- (1-w6[nummodels])*(w[1:(nummodels-1)]/sum(w[1:(nummodels-1)]))
      
      #define stretch weight 2 (maximum is 1)
      w7 <- w
      w7[nummodels] <- (1/w2[nummodels])*w[nummodels]
      w7[1:(nummodels - 1)] <- (1-w7[nummodels])*(w[1:(nummodels-1)]/sum(w[1:(nummodels-1)]))
      
      allweights <- rbind(w,w2,w3,w4,w5,w6,w7)
      
    } else if((numsources == 3) & 
              ((threshold_ind[2] == TRUE) | (threshold_ind[3] == TRUE))
    ){ ## Step 2: Check P(EXCH)_s for external sources (if numsources > 2)
      exch_s_ind <- which(threshold_ind == TRUE)
      
      xval <- as.numeric(max_w$maximum[exch_s_ind])
      
      #get weights at maximum x value
      w2_all <-  calc.MEM.betabin(xvec = c(xval,xvec[-1]), nvec = nvec,
                                  avec = avec, bvec = bvec, prior=prior)
      w2 <- c(w2_all$q)
      
      #calculate weights for P(EXCH) = sqrt(P(EXCH)_{max})
      w3 <- w2
      w3[exch_s_ind] <- sqrt(w2[exch_s_ind])
      w3[c(1:4)[-exch_s_ind]] <- (1-w3[exch_s_ind])*(w2[-exch_s_ind]/sum(w2[-exch_s_ind]))
      
      #define all or nothing borrowing weights
      w4 <- rep(0,nummodels)
      w4[exch_s_ind] <- 1
      
      #get two-source model weights 
      twosource_for_threesource_ind <- 1
      w5_all <-  calc.MEM.betabin(xvec = c(simx,xvec[exch_s_ind]), nvec = nvec[c(1,exch_s_ind)],
                                  avec = avec[c(1,exch_s_ind)], bvec = bvec[c(1,exch_s_ind)], prior=prior)
      w5 <- rep(NA, 4)
      w5[1] <- w5_all$q[1]
      w5[exch_s_ind] <- w5_all$q[2]
      
      #define stretch weight 1 (sqrt(w))
      w6 <- w
      w6[exch_s_ind] <- sqrt(w[exch_s_ind])
      w6[c(1:4)[-exch_s_ind]] <- (1-w6[exch_s_ind])*(w[-exch_s_ind]/sum(w[-exch_s_ind]))
      
      #define stretch weight 2 (maximum is 1)
      w7 <- w
      totalweightforsuppsource <- w2[exch_s_ind] + w2[nummodels]
      propweightforONLYsuppsource <- w2[exch_s_ind]/(w2[nummodels]+w2[exch_s_ind])
      w7[exch_s_ind] <- (1/totalweightforsuppsource)*propweightforONLYsuppsource*(w[exch_s_ind] + w[nummodels])
      w7[c(1:4)[-exch_s_ind]] <- (1-w7[exch_s_ind])*(w[-exch_s_ind]/sum(w[-exch_s_ind]))
      
      allweights <- rbind(w,w2,w3,w4,w5,w6,w7)
      
    } else{
      allweights <- rbind(w,w,w,w,w,w,w) 
    }
    
    for(l in 1:nrow(allweights)){
      
      #FOR THE TWO-SOURCE FOR THREE SOURCE METHOD, ONLY USE DATA FROM THE TWO SOURCES
      if((l == 5) & (twosource_for_threesource_ind == 1)){
        #edit inputs to only include two sources
        avec.o <- avec[c(1,exch_s_ind)]
        bvec.o <- bvec[c(1,exch_s_ind)]
        xvec.o <- c(simx, xvec[exch_s_ind])
        nvec.o <- nvec[c(1,exch_s_ind)]
        mod.mat.o <- res$mod.mat[1:2, 1:2]
        allweights.o <- allweights[l,c(1,exch_s_ind)]
        
      } 
      else{
        avec.o <- avec
        bvec.o <- bvec
        xvec.o <- c(simx, xvec[-1])
        nvec.o <- nvec
        mod.mat.o <- res$mod.mat
        allweights.o <- allweights[l,]
        
      }
      
      ###Calculate MEM mean
      a <- avec.o[1] + mod.mat.o %*% xvec.o #Beta distribution alpha term for each model
      b <- bvec.o[1] + mod.mat.o %*% (nvec.o - xvec.o) #Beta distribution beta term for each model
      model.mean <- a/(a+b) #Posterior mean from Beta distribution
      
      est.var <- (a*b)/((a+b)^2 * (a+b+1)) #Posterior variance from Beta distribution
      
      est.mu[i,l] <- sum(allweights.o * model.mean)
      est.q <- allweights.o #estimated model weights, res$mod.mat for which sources are included in each model
      
      est.bias[i,l] <- est.mu[i,l] - truep
      
      ###Calculate ESSS
      est.esss[i,l] <- sum((BinaryESS2(M = model.mean, V = est.var) - nvec[1]) * allweights.o)  
    }
  }
  
  return(list(simp = simp,
              est.mu = est.mu, 
              est.bias = est.bias,
              est.esss = est.esss
  ))
  
  # est.var = sapply(est.mu, var) #estimate variance of estimated mu-hat from simulation
  # est.mean = sapply(est.mu, mean) #estimate mu-hat
  # est.bias = est.mean - truep #bias
  # mse.part = est.bias^2 + est.var #MSE
  # mse.sum.part = sapply(est.mu, FUN = function(z) (1/niter)*sum((z - truep)^2)) #MSE sum equation
  # mae.part =  sapply(est.mu, FUN = function(z) (1/niter)*sum(abs(z - truep))) #MAE
  # est.esss = sapply(est.esss, median) #ESSS
  # abserr.part = abs(est.bias) #absolute error of bias for summary stat
  # 
  # ###Create output df
  # mat = cbind(truep,est.mean,est.var,est.bias,mse.part,mse.sum.part,mae.part,est.esss,abserr.part)
  # 
  # df <- as.data.frame(mat)
  # df$mu <- xvec[1]/nvec[1]
  # df$num_iterations <- niter
  # df$model <- rownames(df)
  # df$numsources <- numsources
  # if(numsources == 3){
  #   df$diff1 <- xvec[2]/nvec[2] - xvec[1]/nvec[1]
  #   df$diff2 <- xvec[3]/nvec[3] - xvec[1]/nvec[1]
  # }else if(numsources == 2){
  #   df$diff1 <- xvec[2]/nvec[2] - xvec[1]/nvec[1]
  #   df$diff2 <- NA
  # }
  # df$k <- k
  # 
  # return(df)
  
}

##########################################################
####### SECTION 4: RUN TESTS #######
##########################################################


#we expect the absolute bias to ALWAYS BE LOWER for the relaxed methods compared to traditional MEM
#expect all abs(test$est.bias[,2:7]) <= abs(test$est.bias[,1]), 
# so which(abs(test$est.bias[,2:7]) > abs(test$est.bias[,1])) should be empty
#BUT the overall absolute mean bias for the relaxed methods is HIGHER
test_2scont <- TEST_relaxed.MEM_sim_calc(x = c(1,1), S = c(1,1), N = c(30,30), niter = 1000, k = 0.9, 
                                         zeroweightcutoff = 0.01, func_seed = TRUE)
which(abs(test_2scont$est.bias[,2:7]) > abs(test_2scont$est.bias[,1])) #=integer(0)
abs(mean(test_2scont$est.bias[,2])) > abs(mean(test_2scont$est.bias[,1])) #=TRUE
mean(abs(test_2scont$est.bias[,2])) > mean(abs(test_2scont$est.bias[,1])) #=FALSE

test_3scont <- TEST_relaxed.MEM_sim_calc(x = c(1,1,1), S = c(1,1,1), N = c(30,30,30), niter = 1000, k = 0.9, 
                                         zeroweightcutoff = 0.01, func_seed = TRUE)
which(abs(test_3scont$est.bias[,2:7]) > abs(test_3scont$est.bias[,1])) #=integer(0)
abs(mean(test_3scont$est.bias[,2])) > abs(mean(test_3scont$est.bias[,1])) #=FALSE
mean(abs(test_3scont$est.bias[,2])) > mean(abs(test_3scont$est.bias[,1])) #=FALSE

test_2sbin <- TEST_relaxed_sim_calc.betabin(xvec = c(15,15), nvec = c(30,30), avec = c(1,1), bvec = c(1,1),
                                             prior = "equal", niter = 1000, zeroweightcutoff = 0.01, k = 0.9)
which(abs(test_2sbin$est.bias[,2:7]) > abs(test_2sbin$est.bias[,1])) #=integer(0)
abs(mean(test_2sbin$est.bias[,2])) > abs(mean(test_2sbin$est.bias[,1])) #=TRUE
mean(abs(test_2sbin$est.bias[,2])) > mean(abs(test_2sbin$est.bias[,1])) #=FALSE

test_3sbin <- TEST_relaxed_sim_calc.betabin(xvec = c(15,15,15), nvec = c(30,30,30), avec = c(1,1,1), bvec = c(1,1,1),
                                             prior = "equal", niter = 1000, zeroweightcutoff = 0.01, k = 0.9)
which(abs(test_3sbin$est.bias[,2:7]) > abs(test_3sbin$est.bias[,1])) #=integer(0)
abs(mean(test_3sbin$est.bias[,2])) > abs(mean(test_3sbin$est.bias[,1])) #=TRUE
mean(abs(test_3sbin$est.bias[,2])) > mean(abs(test_3sbin$est.bias[,1])) #=FALSE


#compare the distributions
par(mfrow = c(4,4))
hist(test_2scont$est.bias[,1], labels = TRUE,ylim = c(0,550),xlim = c(-0.6,0.4),
     main = paste("2SCont, trad, skewness=",
                  round(skewness(test_2scont$est.bias[,1]),3),
                  ",\nmean=",
                  signif(mean(test_2scont$est.bias[,1]),3),
                  ", median=",
                  signif(median(test_2scont$est.bias[,1]),3)))
hist(test_3scont$est.bias[,1], labels = TRUE,ylim = c(0,550),xlim = c(-0.6,0.4), 
     main = paste("3SCont, trad, skewness=",
                  round(skewness(test_3scont$est.bias[,1]),3),
                  ",\nmean=",
                  signif(mean(test_3scont$est.bias[,1]),3),
                  ", median=",
                  signif(median(test_3scont$est.bias[,1]),3)))
hist(test_2sbin$est.bias[,1], labels = TRUE,ylim = c(0,550),xlim = c(-0.6,0.4), 
     main = paste("2SBin, trad, skewness=",
                  round(skewness(test_2sbin$est.bias[,1]),3),
                  ",\nmean=",
                  signif(mean(test_2sbin$est.bias[,1]),3),
                  ", median=",
                  signif(median(test_2sbin$est.bias[,1]),3)))
hist(test_3sbin$est.bias[,1], labels = TRUE,ylim = c(0,550),xlim = c(-0.6,0.4), 
     main = paste("3SBin, trad, skewness=",
                  round(skewness(test_3sbin$est.bias[,1]),3),
                  ",\nmean=",
                  signif(mean(test_3sbin$est.bias[,1]),3),
                  ", median=",
                  signif(median(test_3sbin$est.bias[,1]),3)))


hist(test_2scont$est.bias[,2], labels = TRUE,ylim = c(0,550),xlim = c(-0.6,0.4), 
     main = paste("2SCont, max, skewness=",
                  round(skewness(test_2scont$est.bias[,2]),3),
                  ",\nmean=",
                  signif(mean(test_2scont$est.bias[,2]),3),
                  ", median=",
                  signif(median(test_2scont$est.bias[,2]),3)))
hist(test_3scont$est.bias[,2], labels = TRUE,ylim = c(0,550),xlim = c(-0.6,0.4), 
     main = paste("3SCont, max, skewness=",
                  round(skewness(test_3scont$est.bias[,2]),3),
                  ",\nmean=",
                  signif(mean(test_3scont$est.bias[,2]),3),
                  ", median=",
                  signif(median(test_3scont$est.bias[,2]),3)))
hist(test_2sbin$est.bias[,2], labels = TRUE,ylim = c(0,550),xlim = c(-0.6,0.4), 
     main = paste("2SBin, max, skewness=",
                  round(skewness(test_2sbin$est.bias[,2]),3),
                  ",\nmean=",
                  signif(mean(test_2sbin$est.bias[,2]),3),
                  ", median=",
                  signif(median(test_2sbin$est.bias[,2]),3)))
hist(test_3sbin$est.bias[,2], labels = TRUE,ylim = c(0,550),xlim = c(-0.6,0.4),
     main = paste("3SBin, max, skewness=",
                  round(skewness(test_3sbin$est.bias[,2]),3),
                  ",\nmean=",
                  signif(mean(test_3sbin$est.bias[,2]),3),
                  ", median=",
                  signif(median(test_3sbin$est.bias[,2]),3)))


hist(test_2scont$est.bias[,7], labels = TRUE,ylim = c(0,550),xlim = c(-0.6,0.4),
     main = paste("2SCont, stretch, skewness=",
                  round(skewness(test_2scont$est.bias[,7]),3),
                  ",\nmean=",
                  signif(mean(test_2scont$est.bias[,7]),3),
                  ", median=",
                  signif(median(test_2scont$est.bias[,7]),3)))
hist(test_3scont$est.bias[,7], labels = TRUE,ylim = c(0,550),xlim = c(-0.6,0.4), 
     main = paste("3SCont, stretch, skewness=",
                  round(skewness(test_3scont$est.bias[,7]),3),
                  ",\nmean=",
                  signif(mean(test_3scont$est.bias[,7]),3),
                  ", median=",
                  signif(median(test_3scont$est.bias[,7]),3)))
hist(test_2sbin$est.bias[,7], labels = TRUE,ylim = c(0,550),xlim = c(-0.6,0.4),
     main = paste("2SBin, stretch, skewness=",
                  round(skewness(test_2sbin$est.bias[,7]),3),
                  ",\nmean=",
                  signif(mean(test_2sbin$est.bias[,7]),3),
                  ", median=",
                  signif(median(test_2sbin$est.bias[,7]),3)))
hist(test_3sbin$est.bias[,7], labels = TRUE,ylim = c(0,550),xlim = c(-0.6,0.4), 
     main = paste("3SBin, stretch, skewness=",
                  round(skewness(test_3sbin$est.bias[,7]),3),
                  ",\nmean=",
                  signif(mean(test_3sbin$est.bias[,7]),3),
                  ", median=",
                  signif(median(test_3sbin$est.bias[,7]),3)))


hist(test_2scont$est.bias[,4], labels = TRUE,ylim = c(0,550),xlim = c(-0.6,0.4), 
     main = paste("2SCont, all, skewness=",
                  round(skewness(test_2scont$est.bias[,4]),3),
                  ",\nmean=",
                  signif(mean(test_2scont$est.bias[,4]),3),
                  ", median=",
                  signif(median(test_2scont$est.bias[,4]),3)))
hist(test_3scont$est.bias[,4], labels = TRUE,ylim = c(0,550),xlim = c(-0.6,0.4),
     main = paste("3SCont, all, skewness=",
                  round(skewness(test_3scont$est.bias[,4]),3),
                  ",\nmean=",
                  signif(mean(test_3scont$est.bias[,4]),3),
                  ", median=",
                  signif(median(test_3scont$est.bias[,4]),3)))
hist(test_2sbin$est.bias[,4], labels = TRUE,ylim = c(0,550),xlim = c(-0.6,0.4), 
     main = paste("2SBin, all, skewness=",
                  round(skewness(test_2sbin$est.bias[,4]),3),
                  ",\nmean=",
                  signif(mean(test_2sbin$est.bias[,4]),3),
                  ", median=",
                  signif(median(test_2sbin$est.bias[,4]),3)))
hist(test_3sbin$est.bias[,4], labels = TRUE,ylim = c(0,550),xlim = c(-0.6,0.4),
     main = paste("3SBin, all, skewness=",
                  round(skewness(test_3sbin$est.bias[,4]),3),
                  ",\nmean=",
                  signif(mean(test_3sbin$est.bias[,4]),3),
                  ", median=",
                  signif(median(test_3sbin$est.bias[,4]),3)))


### ABOLUTE BIAS PLOT
par(mfrow = c(4,4))
hist(abs(test_2scont$est.bias[,1]), labels = TRUE,ylim = c(0,700),xlim = c(0,0.7),
     main = paste("2SCont, trad, skewness=",
                  round(skewness(abs(test_2scont$est.bias[,1])),3),
                  ",\nmean=",
                  signif(mean(abs(test_2scont$est.bias[,1])),3),
                  ", median=",
                  signif(median(abs(test_2scont$est.bias[,1])),3)))
hist(abs(test_3scont$est.bias[,1]), labels = TRUE,ylim = c(0,700),xlim = c(0,0.7), 
     main = paste("3SCont, trad, skewness=",
                  round(skewness(abs(test_3scont$est.bias[,1])),3),
                  ",\nmean=",
                  signif(mean(abs(test_3scont$est.bias[,1])),3),
                  ", median=",
                  signif(median(abs(test_3scont$est.bias[,1]),3))))
hist(abs(test_2sbin$est.bias[,1]), labels = TRUE,ylim = c(0,700),xlim = c(0,0.7), 
     main = paste("2SBin, trad, skewness=",
                  round(skewness(abs(test_2sbin$est.bias[,1])),3),
                  ",\nmean=",
                  signif(mean(abs(test_2sbin$est.bias[,1])),3),
                  ", median=",
                  signif(median(abs(test_2sbin$est.bias[,1])),3)))
hist(abs(test_3sbin$est.bias[,1]), labels = TRUE,ylim = c(0,700),xlim = c(0,0.7), 
     main = paste("3SBin, trad, skewness=",
                  round(skewness(abs(test_3sbin$est.bias[,1])),3),
                  ",\nmean=",
                  signif(mean(abs(test_3sbin$est.bias[,1])),3),
                  ", median=",
                  signif(median(abs(test_3sbin$est.bias[,1])),3)))


hist(abs(test_2scont$est.bias[,2]), labels = TRUE,ylim = c(0,700),xlim = c(0,0.7), 
     main = paste("2SCont, max, skewness=",
                  round(skewness(abs(test_2scont$est.bias[,2])),3),
                  ",\nmean=",
                  signif(mean(abs(test_2scont$est.bias[,2])),3),
                  ", median=",
                  signif(median(abs(test_2scont$est.bias[,2])),3)))
hist(abs(test_3scont$est.bias[,2]), labels = TRUE,ylim = c(0,700),xlim = c(0,0.7), 
     main = paste("3SCont, max, skewness=",
                  round(skewness(abs(test_3scont$est.bias[,2])),3),
                  ",\nmean=",
                  signif(mean(abs(test_3scont$est.bias[,2])),3),
                  ", median=",
                  signif(median(abs(test_3scont$est.bias[,2])),3)))
hist(abs(test_2sbin$est.bias[,2]), labels = TRUE,ylim = c(0,700),xlim = c(0,0.7), 
     main = paste("2SBin, max, skewness=",
                  round(skewness(abs(test_2sbin$est.bias[,2])),3),
                  ",\nmean=",
                  signif(mean(abs(test_2sbin$est.bias[,2])),3),
                  ", median=",
                  signif(median(abs(test_2sbin$est.bias[,2])),3)))
hist(abs(test_3sbin$est.bias[,2]), labels = TRUE,ylim = c(0,700),xlim = c(0,0.7), 
     main = paste("3SBin, max, skewness=",
                  round(skewness(abs(test_3sbin$est.bias[,2])),3),
                  ",\nmean=",
                  signif(mean(abs(test_3sbin$est.bias[,2])),3),
                  ", median=",
                  signif(median(abs(test_3sbin$est.bias[,2])),3)))


hist(abs(test_2scont$est.bias[,7]), labels = TRUE,ylim = c(0,700),xlim = c(0,0.7), 
     main = paste("2SCont, stretch, skewness=",
                  round(skewness(abs(test_2scont$est.bias[,7])),3),
                  ",\nmean=",
                  signif(mean(abs(test_2scont$est.bias[,7])),3),
                  ", median=",
                  signif(median(abs(test_2scont$est.bias[,7])),3)))
hist(abs(test_3scont$est.bias[,7]), labels = TRUE,ylim = c(0,700),xlim = c(0,0.7), 
     main = paste("3SCont, stretch, skewness=",
                  round(skewness(abs(test_3scont$est.bias[,7])),3),
                  ",\nmean=",
                  signif(mean(abs(test_3scont$est.bias[,7])),3),
                  ", median=",
                  signif(median(abs(test_3scont$est.bias[,7])),3)))
hist(abs(test_2sbin$est.bias[,7]), labels = TRUE,ylim = c(0,700),xlim = c(0,0.7),
     main = paste("2SBin, stretch, skewness=",
                  round(skewness(abs(test_2sbin$est.bias[,7])),3),
                  ",\nmean=",
                  signif(mean(abs(test_2sbin$est.bias[,7])),3),
                  ", median=",
                  signif(median(abs(test_2sbin$est.bias[,7])),3)))
hist(abs(test_3sbin$est.bias[,7]), labels = TRUE,ylim = c(0,700),xlim = c(0,0.7), 
     main = paste("3SBin, stretch, skewness=",
                  round(skewness(abs(test_3sbin$est.bias[,7])),3),
                  ",\nmean=",
                  signif(mean(abs(test_3sbin$est.bias[,7])),3),
                  ", median=",
                  signif(median(abs(test_3sbin$est.bias[,7])),3)))


hist(abs(test_2scont$est.bias[,4]), labels = TRUE,ylim = c(0,700),xlim = c(0,0.7), 
     main = paste("2SCont, all, skewness=",
                  round(skewness(abs(test_2scont$est.bias[,4])),3),
                  ",\nmean=",
                  signif(mean(abs(test_2scont$est.bias[,4])),3),
                  ", median=",
                  signif(median(abs(test_2scont$est.bias[,4])),3)))
hist(abs(test_3scont$est.bias[,4]), labels = TRUE,ylim = c(0,700),xlim = c(0,0.7),
     main = paste("3SCont, all, skewness=",
                  round(skewness(abs(test_3scont$est.bias[,4])),3),
                  ",\nmean=",
                  signif(mean(abs(test_3scont$est.bias[,4])),3),
                  ", median=",
                  signif(median(abs(test_3scont$est.bias[,4])),3)))
hist(abs(test_2sbin$est.bias[,4]), labels = TRUE,ylim = c(0,700),xlim = c(0,0.7), 
     main = paste("2SBin, all, skewness=",
                  round(skewness(abs(test_2sbin$est.bias[,4])),3),
                  ",\nmean=",
                  signif(mean(abs(test_2sbin$est.bias[,4])),3),
                  ", median=",
                  signif(median(abs(test_2sbin$est.bias[,4])),3)))
hist(abs(test_3sbin$est.bias[,4]), labels = TRUE,ylim = c(0,700),xlim = c(0,0.7),
     main = paste("3SBin, all, skewness=",
                  round(skewness(abs(test_3sbin$est.bias[,4])),3),
                  ",\nmean=",
                  signif(mean(abs(test_3sbin$est.bias[,4])),3),
                  ", median=",
                  signif(median(abs(test_3sbin$est.bias[,4])),3)))


