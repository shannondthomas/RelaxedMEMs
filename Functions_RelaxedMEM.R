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

library(gtools)
library(matrixStats)
library(xtable)
library(tidyverse)



##########################################################
########### SECTION 1: CONTINUOUS MEM FUNCTIONS ##########
##########################################################

boa.hpd <- function(x, alpha){
  ###Calculate HPD Intervals given vector of values and alpha
  n <- length(x)
  m <- max(1, ceiling(alpha * n))
  y <- sort(x)
  a <- y[1:m]
  b <- y[(n - m + 1):n]
  i <- order(b - a)[1]
  structure(c(a[i], b[i]), names = c("Lower Bound", "Upper Bound"))
}

#function to calculate MEM weights for any number of sources
general.calc.weights_MEM = function(xvec,svec,nvec,prior = "pi_e"){
  ###function to calculate model weights for MEM approach with "correct" calculations which don't assume conditional independence
  #xvec: means for sources
  #svec: standard deviation for sources
  #nvec: sample size for sources
  #prior: prior to use for calculations
  
  numsources <- length(xvec)
  H <- numsources - 1                              #number of supplemental sources
  s_exchindicator <- expand.grid(rep(list(0:1),H)) #data frame with indicators for supp source exchangeability for each model
  nummodels <- 2^H                                 #number of models = nrow(s_exchindicator)
  
  vvec <- (svec^2)/nvec             #define v vector
  prob_model <- rep(NA, nummodels)  #initialize weight vector
  expderiv <- rep(NA, nummodels)      #initialize derivative vector
  
  suppvvec <- vvec[2:length(vvec)]
  suppxvec <- xvec[2:length(xvec)]
  
  #loop through each exchangeability model to calculate weights individually
  for(i in 1:nummodels){
    
    evec <- s_exchindicator[i,] #define exchangeability vector
    
    
    term11 <- (sqrt(2*pi))^(H+1 - sum(evec))
    term12 <- (1/vvec[1]) + sum(evec/suppvvec)
    term13 <- prod((1/suppvvec)^(1-evec))
    
    term1 <- term11/sqrt(term12*term13)
    
    #EXPONENTIAL TERM
    #vector of length H giving sum_{m \neq l} evec[m]/suppvvec[m]
    term21_denom_sum1 <- as.numeric(lapply(as.data.frame(abs(diag(nrow = H)-1)), 
                                           FUN = function(x){sum(x*evec/suppvvec)}))
    #vector of length H giving the first denominator in the exponential func
    term21_denom1 <- vvec[1] + suppvvec + vvec[1]*suppvvec*term21_denom_sum1
    
    #vector of length H giving evec(xvec[1] - xvec[i])^2
    term21_num1 <- evec*as.numeric((xvec[1] - suppxvec%*%diag(nrow = H))^2)
    
    #derivative numerator - vector of length H giving evec(xvec[1] - xvec[i])
    term21_num1_deriv <- evec*as.numeric(xvec[1] - suppxvec%*%diag(nrow = H))
    
    #vector of length H for first term in exponential func
    term21 <- term21_num1/term21_denom1
    
    #define pairs of supp. sources that are exchangeable for term 2
    exch_pairs <- t(as.matrix(evec))%*%as.matrix(evec)
    #keep only the lower triangle
    exch_pairs_lt <- lower.tri(exch_pairs)*exch_pairs
    #get list of sources that will be in term22
    term22_sources <- which(exch_pairs_lt != 0, arr.ind = TRUE)
    
    #vector of length = nrow(term22_sources) giving (x_l - x_r)^2
    term22_num <- as.numeric(lapply(as.data.frame(t(term22_sources)),
                                    FUN = function(x){(suppxvec[x[1]] - suppxvec[x[2]])^2}))
    
    #vector of length = nrow(term22_sources) giving (x_l - x_r)^2
    term22_den <- as.numeric(lapply(as.data.frame(t(term22_sources)),
                                    FUN = function(x){suppvvec[x[1]] + suppvvec[x[2]] + 
                                        suppvvec[x[1]]*suppvvec[x[2]]*((1/vvec[1]) + sum(evec[-x]/suppvvec[-x]))}))
    #vector of length = nrow(term22_sources)
    term22 <- term22_num/term22_den
    
    
    prob_model[i] <- term1*exp(-0.5*(sum(term21) + sum(term22)))
    
    #Derivative of exponential part of each w_i weight component for delta method variance calculation
    ##Note that the exp() part is left off because it is contained within the w_i part included in the next step
    ##(i.e., (w1+w2+w3+...+w8)/w2 will contain all the necessary info besides the deriv of the exp part calculated for b here)
    expderiv[i] <- sum(-1*term21_num1_deriv/term21_denom1)
    
  }                                                                                                                      
  
  if(length(prior) == 1){
    if(prior == "pi_e"){
      #note: don't need to use prior here since they're all equal to 1/2 and cancel out
      weights <- prob_model/sum(prob_model)
    }else if(prior != "pi_e"){stop("Prior must be pi_e or a vector of length 2^H.")}
  }else if (length(prior) == nummodels){
    weights <- prior*prob_model/sum(prior*prob_model)
  }else{stop("Prior must be pi_e or a vector of length 2^H.")}
  
  return(data.frame(weight = weights, expderiv = expderiv))
  
}

###Run all continuous MEMs (relaxed and original) and summarize performance
relaxed.MEM_sim_calc = function(x,S,N,prior = "pi_e",niter,k,zeroweightcutoff,RAND=TRUE,func_seed=FALSE){
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
  
  
  est.var = sapply(est.mu, var) #estimate variance of estimated mu-hat from simulation
  if((RAND==FALSE)&(niter==1)){
    est.var <- var.out
  }
  est.mean = sapply(est.mu, mean) #estimate mu-hat
  est.bias = est.mean - mu #bias
  mse.part = est.bias^2 + est.var #MSE
  mse.sum.part = sapply(est.mu, FUN = function(z) (1/niter)*sum((z - mu)^2)) #MSE sum equation
  mae.part =  sapply(est.mu, FUN = function(z) (1/niter)*sum(abs(z - mu))) #MAE
  est.esss = sapply(est.esss, median) #ESSS
  abserr.part = abs(est.bias) #absolute error of bias for summary stat
  cov.part = sapply(cov.ind, mean) #coverage
  hpd.part = sapply(hpd.width, mean) #HPD width
  
  mat = cbind(mu,est.mean,est.var,est.bias,mse.part,mse.sum.part,mae.part,est.esss,abserr.part,cov.part,hpd.part)
  
  df <- as.data.frame(mat)
  df$num_iterations <- niter
  df$model <- rownames(df)
  df$numsources <- numsources
  if(numsources == 3){
    df$diff1 <- x[2] - x[1]
    df$diff2 <- x[3] - x[1]
  }else if(numsources == 2){
    df$diff1 <- x[2] - x[1]
    df$diff2 <- NA
  }
  df$k <- k
  return(df)
}


###Function to determine x value at which k=1 is better than k=0.
findcriticalpoints <- function(x,S,N,prior = "pi_e",niter,zeroweightcutoff,RAND=TRUE,stepsize = NA,MSEtolerance = 0.002){
  if(is.na(stepsize)){
    stepsize <- min(S)/10
  }
  
  diffx <- x[-1] - x[1]
  
  minx <- x[1] + min(diffx) - 1.5*max(S)
  maxx <- x[1] + max(diffx) + 1.5*max(S)
  
  numx <- round((maxx - minx)/stepsize,0)
  
  df <- data.frame(primarymean = rep(NA,7*numx),
                   model = rep(NA,7*numx),
                   MSE_k0 = rep(NA,7*numx),
                   MSE_k1 = rep(NA,7*numx),
                   k1gtk0 = rep(NA,7*numx))
  
  criticalvals <- data.frame(model = c("MEM",
                                       "relaxedMEM_max",
                                       "relaxedMEM_sqrtmax",
                                       "relaxedMEM_allornothing",
                                       "relaxedMEM_twosource",
                                       "relaxedMEM_stretchsqrt",
                                       "relaxedMEM_stretch1"),
                             min1 = rep(NA,7),
                             max1 = rep(NA,7),
                             min2 = rep(NA,7),
                             max2 = rep(NA,7),
                             min3 = rep(NA,7),
                             max3 = rep(NA,7),
                             min_simx = rep(minx,7),
                             max_simx = rep(maxx,7)
                             )
  
  newx <- minx
  pb <- txtProgressBar(min = 0,      
                       max = numx, 
                       style = 3,    
                       width = 50,  
                       char = "=")   
  system.time({
  for(i in 1:numx){
    k0 <- relaxed.MEM_sim_calc(x=c(newx, x[-1]),S=S,N=N,prior = "pi_e",niter=niter,k=0,zeroweightcutoff=zeroweightcutoff,RAND=RAND)
    k1 <- relaxed.MEM_sim_calc(x=c(newx, x[-1]),S=S,N=N,prior = "pi_e",niter=niter,k=1,zeroweightcutoff=zeroweightcutoff,RAND=RAND)
    
    df$MSE_k0[(7*(i-1)+1):(7*(i-1)+7)] <- k0$mse.part
    df$MSE_k1[(7*(i-1)+1):(7*(i-1)+7)] <- k1$mse.part
    df$model[(7*(i-1)+1):(7*(i-1)+7)] <- k1$model
    df$primarymean[(7*(i-1)+1):(7*(i-1)+7)] <- newx
    
    criticalvals$min1[(is.na(criticalvals$min1))&
                        (df$MSE_k1[(7*(i-1)+1):(7*(i-1)+7)] >= (df$MSE_k0[(7*(i-1)+1):(7*(i-1)+7)] - MSEtolerance))] <- newx
    criticalvals$max1[(!is.na(criticalvals$min1))&(is.na(criticalvals$max1))&
                        (df$MSE_k1[(7*(i-1)+1):(7*(i-1)+7)] < (df$MSE_k0[(7*(i-1)+1):(7*(i-1)+7)] - MSEtolerance))] <- newx-stepsize
    criticalvals$min2[(!is.na(criticalvals$max1))&(is.na(criticalvals$min2))&
                        (df$MSE_k1[(7*(i-1)+1):(7*(i-1)+7)] >= (df$MSE_k0[(7*(i-1)+1):(7*(i-1)+7)] - MSEtolerance))] <- newx
    criticalvals$max2[(!is.na(criticalvals$min2))&(is.na(criticalvals$max2))&
                        (df$MSE_k1[(7*(i-1)+1):(7*(i-1)+7)] < (df$MSE_k0[(7*(i-1)+1):(7*(i-1)+7)] - MSEtolerance))] <- newx-stepsize
    criticalvals$min3[(!is.na(criticalvals$max2))&(is.na(criticalvals$min3))&
                        (df$MSE_k1[(7*(i-1)+1):(7*(i-1)+7)] >= (df$MSE_k0[(7*(i-1)+1):(7*(i-1)+7)] - MSEtolerance))] <- newx
    criticalvals$max3[(!is.na(criticalvals$min3))&(is.na(criticalvals$max3))&
                        (df$MSE_k1[(7*(i-1)+1):(7*(i-1)+7)] < (df$MSE_k0[(7*(i-1)+1):(7*(i-1)+7)] - MSEtolerance))] <- newx-stepsize
    
    newx <- newx + stepsize
    
    setTxtProgressBar(pb, i)
  }
  })
  
  criticalvals$max1 <- ifelse(is.na(criticalvals$max1), maxx, criticalvals$max1)
  criticalvals$max2 <- ifelse(is.na(criticalvals$max2) & (!is.na(criticalvals$min2)), maxx, criticalvals$max2)
  criticalvals$max3 <- ifelse(is.na(criticalvals$max3) & (!is.na(criticalvals$min3)), maxx, criticalvals$max3)
  
  # df$k1gtek0 <- df$MSE_k1 >= df$MSE_k0
  close(pb)
  
  # #wide df to check that there is a continuous range of primarymean for which k1gtek0 == TRUE
  # df_wide <- pivot_wider(df, id_cols = c(primarymean),
  #                        names_from= model,
  #                        values_from = c(MSE_k0, MSE_k1,k1gtek0),
  #                        names_glue = "{model}_{.value}", names_vary = "slowest")
  # 
  # #plot difference in MSEs
  # ggplot(data = df, aes(x = primarymean, color = model, group = model)) +
  #   geom_point(aes(y=MSE_k1 - MSE_k0)) +
  #   geom_line(aes(y=MSE_k1 - MSE_k0)) +
  #   geom_vline(xintercept = x[2], color = "red") + 
  #   geom_vline(xintercept = x[3], color = "red") +
  #   ggtitle(paste("Difference in MSE for never relaxing (k=1) or always relaxing (k=0) \nPositive -> Relaxed MEMs are better \n SS1 Mean =",x[2],", SS2 Mean =",x[3])) +
  #   ylim(-0.1,0.1)
  
  
  return(criticalvals)
}

###Run all continuous MEMs (relaxed and original) and summarize performance GIVEN CRITICAL POINTS
relaxed.MEM_sim_calc_critpts = function(x,S,N,prior = "pi_e",niter,zeroweightcutoff,RAND=TRUE,criticalpoints=NULL){
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
    set.seed(515+i)
    if((RAND == TRUE) & (niter != 1)){
      xbar = rnorm(1,mean=mu,sd=S[1]/sqrt(N[1])) 
    }else{
      xbar <- x[1]
    }
    V <- S^2/N

    ##Calculate weights and derivative terms needed for variance calculation
    w_and_der <- general.calc.weights_MEM(xvec=c(xbar,x[-1]),svec=S,nvec=N,prior=prior)
    w <- w_and_der$weight
    b <- w_and_der$expderiv
    
    # Remove all-zero row (already covered by M1 = xbar) and not needed for weight replacement
    s_exchindicator <- s_exchindicator[rowSums(s_exchindicator) > 0, , drop = FALSE]
    
    #############RELAXED MEM ADDITION
    ### Create indicator for if xbar is in critical range for each method 
    criticalrange_ind <- ifelse(((xbar >= criticalpoints$min1)&(xbar <= criticalpoints$max1))|
                                ((xbar >= criticalpoints$min2)&(xbar <= criticalpoints$max2))|
                                ((xbar >= criticalpoints$min3)&(xbar <= criticalpoints$max3)), TRUE, FALSE)
    
    ## If any lie in the critical range, then calculate relaxed weights.
    if(sum(criticalrange_ind,na.rm = TRUE)){
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
      
      allweights <- rbind(w,
                          if(criticalrange_ind[2] %in% TRUE) w2 else(w),
                          if(criticalrange_ind[3] %in% TRUE) w3 else(w),
                          if(criticalrange_ind[4] %in% TRUE) w4 else(w),
                          if(criticalrange_ind[5] %in% TRUE) w5 else(w),
                          if(criticalrange_ind[6] %in% TRUE) w6 else(w),
                          if(criticalrange_ind[7] %in% TRUE) w7 else(w))
      
    } else{
      allweights <- rbind(w,w,w,w,w,w,w) 
    }
    
    
    
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
  
  
  est.var = sapply(est.mu, var) #estimate variance of estimated mu-hat from simulation
  if((RAND==FALSE)&(niter==1)){
    est.var <- var.out
  }
  est.mean = sapply(est.mu, mean) #estimate mu-hat
  est.bias = est.mean - mu #bias
  mse.part = est.bias^2 + est.var #MSE
  mse.sum.part = sapply(est.mu, FUN = function(z) (1/niter)*sum((z - mu)^2)) #MSE sum equation
  mae.part =  sapply(est.mu, FUN = function(z) (1/niter)*sum(abs(z - mu))) #MAE
  est.esss = sapply(est.esss, median) #ESSS
  abserr.part = abs(est.bias) #absolute error of bias for summary stat
  cov.part = sapply(cov.ind, mean) #coverage
  hpd.part = sapply(hpd.width, mean) #HPD width
  
  mat = cbind(mu,est.mean,est.var,est.bias,mse.part,mse.sum.part,mae.part,est.esss,abserr.part,cov.part,hpd.part,criticalpoints)
  
  df <- as.data.frame(mat)
  df$num_iterations <- niter
  df$model <- rownames(df)
  df$numsources <- numsources
  if(numsources == 3){
    df$diff1 <- x[2] - x[1]
    df$diff2 <- x[3] - x[1]
  }else if(numsources == 2){
    df$diff1 <- x[2] - x[1]
    df$diff2 <- NA
  }
  return(df)
}


##########################################################
############# SECTION 3: BINARY MEM FUNCTIONS ############
##########################################################

###Function to calculate ESS for binary data
BinaryESS2 <- function(V,M){
  #M: posterior mean
  #V: posterior variance
  
  u <- (M*(1-M))/V - 1
  a <- M*u
  b <- (1-M)*u
  return(a+b)
  
}

###Function to calculate the MEM model weights for binomial data with beta(alpha,beta) prior
calc.MEM.betabin <- function(xvec, nvec, avec, bvec, prior, constraint=1){
  #xvec: vector with counts of those with event (primary source first, supplemental afterwards)
  #nvec: vector of total number of subjects
  #avec: vector of alpha parameters for beta priors for each source
  #bvec: vector of beta parameters for beta priors for each source
  #prior: prior to use on MEM source inclusion probability: equal (pi_e), pool (naive pooling), opt.source (pi_EB), opt.source_constrain (pi_EBc)
  #constraint: value to limit maximum value of empirical Bayes prior to for pi_EBc
  
  mod.mat <- as.matrix(expand.grid( rep(list(c(0,1)), length(xvec)-1) )) #create matrix of all models with potential combinations of supplemental sources
  mod.mat <- mod.mat[order(rowSums(mod.mat)),] #group by number of supplemental sources
  mod.mat <- cbind(1, mod.mat) #add column before 1st column for primary source indicator
  colnames(mod.mat) <- c('p',paste0('s',seq(1,length(xvec)-1)))
  
  prod.vec <- beta(xvec + avec, nvec + bvec - xvec) / beta(avec, bvec) #calculate the product portion of integrated marginal likelihood
  p.vec <- apply( t(sapply(1:dim(mod.mat)[1], function(x) prod.vec^(1-mod.mat[x,]))), 1, prod) #calculate the product portion of the integrated marginal likelihood corresponding to each model by taking power w/indicator for each source in a model and then using apply over rows of resulting matrix
  
  ###Calculate the integrated marginal likelihood given data:
  marg.vec <- (beta(avec[1] + mod.mat%*%xvec, bvec[1] + mod.mat%*%(nvec-xvec)) / beta(avec[1],bvec[1]) ) * p.vec
  
  ###Calculate prior:
  if(prior=='equal'){
    prior1 <- rep(.5, length(xvec)-1)
    prior0 <- rep(.5, length(xvec)-1)
  }else if(prior=='pool'){
    prior1 <- rep(1, length(xvec)-1)
    prior0 <- rep(0, length(xvec)-1)
  }else if(prior=='opt.source'){
    #this prior identifies the optimal MEM and gives it a weight of 1
    
    s.mat_opt.source <- mod.mat[,-1]
    min.vec <- -log( marg.vec ) 
    if(length(xvec)>2){
      if( length(which(min.vec == min(min.vec))) > 1){
        prior1.sum <- colSums(s.mat_opt.source[ which(min.vec == min(min.vec)) ,])
        prior1.sum[which(prior1.sum != 0)] <- 1
        prior1 <- prior1.sum
      }else{
        prior1 <- s.mat_opt.source[ which(min.vec == min(min.vec)) ,]
      }
      prior0 <- 1 - prior1
    }else{
      prior1 <- s.mat_opt.source[ which(min.vec == min(min.vec)) ]
      prior0 <- 1 - prior1
    }
  }else if(prior=='opt.source_constrain'){
    #this prior identifies the optimal MEM and gives it a weight of 1
    
    s.mat_opt.source <- mod.mat[,-1]
    min.vec <- -log( marg.vec ) 
    if(length(xvec)>2){
      if( length(which(min.vec == min(min.vec))) > 1){
        prior1.sum <- colSums(s.mat_opt.source[ which(min.vec == min(min.vec)) ,])
        prior1.sum[which(prior1.sum != 0)] <- 1
        prior1 <- prior1.sum*constraint
      }else{
        prior1 <- s.mat_opt.source[ which(min.vec == min(min.vec)) ,]*constraint
      }
      prior0 <- 1 - prior1
    }else{
      prior1 <- s.mat_opt.source[ which(min.vec == min(min.vec)) ]*constraint
      prior0 <- 1 - prior1
    }
  }else{print('Prior not valid, please enter a valid prior option.')}
  
  ###Calculate model weights given priors
  
  if(length(xvec)==2){mps <- matrix( rbind(prior0,prior1)[ paste0('prior',(mod.mat[,2:length(xvec)])), 1], ncol=1)} #extract priors for sources for each MEM
  if(length(xvec)>2){mps <- sapply(1:(length(xvec)-1), function(x) rbind(prior0,prior1)[ paste0('prior',(mod.mat[,2:length(xvec)])[,x]), x])} #extract priors for sources for each MEM
  
  q.vec <- marg.vec * ( rowProds(mps)/sum(rowProds(mps)) ) / sum(marg.vec * ( rowProds(mps)/sum(rowProds(mps)) ) ) #model weights, note: need to include sum(rowProds(mps)) otherwise some cases result in NaN because values are so small
  
  ret <- list(q = q.vec, mod.mat = mod.mat, prior = prior )
  
  return(ret)
}


###Run all binary MEMs (relaxed and original) and summarize performance
relaxed_sim_calc.betabin <- function(xvec,nvec,avec,bvec,prior,constraint=1, niter, zeroweightcutoff, k, RAND=TRUE){
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
  
  est.var = sapply(est.mu, var) #estimate variance of estimated mu-hat from simulation
  est.mean = sapply(est.mu, mean) #estimate mu-hat
  est.bias = est.mean - truep #bias
  mse.part = est.bias^2 + est.var #MSE
  mse.sum.part = sapply(est.mu, FUN = function(z) (1/niter)*sum((z - truep)^2)) #MSE sum equation
  mae.part =  sapply(est.mu, FUN = function(z) (1/niter)*sum(abs(z - truep))) #MAE
  est.esss = sapply(est.esss, median) #ESSS
  abserr.part = abs(est.bias) #absolute error of bias for summary stat
  
  ###Create output df
  mat = cbind(truep,est.mean,est.var,est.bias,mse.part,mse.sum.part,mae.part,est.esss,abserr.part)
  
  df <- as.data.frame(mat)
  df$mu <- xvec[1]/nvec[1]
  df$num_iterations <- niter
  df$model <- rownames(df)
  df$numsources <- numsources
  if(numsources == 3){
    df$diff1 <- xvec[2]/nvec[2] - xvec[1]/nvec[1]
    df$diff2 <- xvec[3]/nvec[3] - xvec[1]/nvec[1]
  }else if(numsources == 2){
    df$diff1 <- xvec[2]/nvec[2] - xvec[1]/nvec[1]
    df$diff2 <- NA
  }
  df$k <- k
  
  return(df)
  
}



###Function to determine x value at which k=1 is better than k=0.
findcriticalpoints_binary <- function(xvec,nvec,avec,bvec,prior="equal",constraint=1,niter,zeroweightcutoff,RAND=TRUE,MSEtolerance=0.002){
  
  minx <- 0
  maxx <- nvec[1]
  numx <- nvec[1]+1
  
  df <- data.frame(primarymean = rep(NA,7*numx),
                   model = rep(NA,7*numx),
                   MSE_k0 = rep(NA,7*numx),
                   MSE_k1 = rep(NA,7*numx),
                   k1gtk0 = rep(NA,7*numx))
  
  criticalvals <- data.frame(model = c("MEM",
                                       "relaxedMEM_max",
                                       "relaxedMEM_sqrtmax",
                                       "relaxedMEM_allornothing",
                                       "relaxedMEM_twosource",
                                       "relaxedMEM_stretchsqrt",
                                       "relaxedMEM_stretch1"),
                             min1 = rep(NA,7),
                             max1 = rep(NA,7),
                             min2 = rep(NA,7),
                             max2 = rep(NA,7),
                             min3 = rep(NA,7),
                             max3 = rep(NA,7),
                             min_simx = rep(minx,7),
                             max_simx = rep(maxx,7)
  )
  
  newx <- minx
  pb <- txtProgressBar(min = 0,      
                       max = numx, 
                       style = 3,    
                       width = 50,  
                       char = "=")   
  system.time({
    for(i in 1:numx){
      k0 <- relaxed_sim_calc.betabin(xvec = c(newx,xvec[-1]), nvec = nvec, 
                                     avec = avec, bvec = bvec,
                                     prior = "equal", constraint = constraint,
                                     k = 0, zeroweightcutoff = zeroweightcutoff,
                                     niter = niter, RAND=RAND)
      k1 <- relaxed_sim_calc.betabin(xvec = c(newx,xvec[-1]), nvec = nvec, 
                                     avec = avec, bvec = bvec,
                                     prior = "equal", constraint = constraint,
                                     k = 1, zeroweightcutoff = zeroweightcutoff,
                                     niter = niter, RAND=RAND)
      
      df$MSE_k0[(7*(i-1)+1):(7*(i-1)+7)] <- k0$mse.part
      df$MSE_k1[(7*(i-1)+1):(7*(i-1)+7)] <- k1$mse.part
      df$model[(7*(i-1)+1):(7*(i-1)+7)] <- k1$model
      df$primarymean[(7*(i-1)+1):(7*(i-1)+7)] <- newx/nvec[1]
      
      criticalvals$min1[(is.na(criticalvals$min1))&
                          (df$MSE_k1[(7*(i-1)+1):(7*(i-1)+7)] >= (df$MSE_k0[(7*(i-1)+1):(7*(i-1)+7)] - MSEtolerance))] <- newx
      criticalvals$max1[(!is.na(criticalvals$min1))&(is.na(criticalvals$max1))&
                          (df$MSE_k1[(7*(i-1)+1):(7*(i-1)+7)] < (df$MSE_k0[(7*(i-1)+1):(7*(i-1)+7)] - MSEtolerance))] <- newx-1
      criticalvals$min2[(!is.na(criticalvals$max1))&(is.na(criticalvals$min2))&
                          (df$MSE_k1[(7*(i-1)+1):(7*(i-1)+7)] >= (df$MSE_k0[(7*(i-1)+1):(7*(i-1)+7)] - MSEtolerance))] <- newx
      criticalvals$max2[(!is.na(criticalvals$min2))&(is.na(criticalvals$max2))&
                          (df$MSE_k1[(7*(i-1)+1):(7*(i-1)+7)] < (df$MSE_k0[(7*(i-1)+1):(7*(i-1)+7)] - MSEtolerance))] <- newx-1
      criticalvals$min3[(!is.na(criticalvals$max2))&(is.na(criticalvals$min3))&
                          (df$MSE_k1[(7*(i-1)+1):(7*(i-1)+7)] >= (df$MSE_k0[(7*(i-1)+1):(7*(i-1)+7)] - MSEtolerance))] <- newx
      criticalvals$max3[(!is.na(criticalvals$min3))&(is.na(criticalvals$max3))&
                          (df$MSE_k1[(7*(i-1)+1):(7*(i-1)+7)] < (df$MSE_k0[(7*(i-1)+1):(7*(i-1)+7)] - MSEtolerance))] <- newx-1
      
      newx <- newx + 1
      
      setTxtProgressBar(pb, i)
    }
  })
  
  criticalvals$max1 <- ifelse(is.na(criticalvals$max1), maxx, criticalvals$max1)
  criticalvals$max2 <- ifelse(is.na(criticalvals$max2) & (!is.na(criticalvals$min2)), maxx, criticalvals$max2)
  criticalvals$max3 <- ifelse(is.na(criticalvals$max3) & (!is.na(criticalvals$min3)), maxx, criticalvals$max3)
  
  df$k1gtek0 <- df$MSE_k1 >= df$MSE_k0
  close(pb)
  
  # #wide df to check that there is a continuous range of primarymean for which k1gtek0 == TRUE
  # df_wide <- pivot_wider(df, id_cols = c(primarymean),
  #                        names_from= model,
  #                        values_from = c(MSE_k0, MSE_k1,k1gtek0),
  #                        names_glue = "{model}_{.value}", names_vary = "slowest")
  # 
  # #plot difference in MSEs
  # ggplot(data = df, aes(x = primarymean, color = model, group = model)) +
  #   geom_point(aes(y=MSE_k1 - MSE_k0)) +
  #   geom_line(aes(y=MSE_k1 - MSE_k0)) +
  #   geom_vline(xintercept = xvec[2]/nvec[2], color = "red") +
  #   geom_vline(xintercept = xvec[3]/nvec[3], color = "red") +
  #   ggtitle(paste("Difference in MSE for never relaxing (k=1) or always relaxing (k=0) \nPositive -> Relaxed MEMs are better \n SS1 Mean =",round(xvec[2]/nvec[2],2),", SS2 Mean =",round(xvec[3]/nvec[3],2))) #+
  #   ylim(-0.01,0.01)
  
  
  return(criticalvals)
}



###Run all binary MEMs (relaxed and original) and summarize performance
relaxed_sim_calc.betabin_critpts <- function(xvec,nvec,avec,bvec,prior,constraint=1, niter, zeroweightcutoff, RAND=TRUE, criticalpoints){
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
    
    #Calculate standard MEM weights
    res <- calc.MEM.betabin(xvec=c(simx, xvec[-1]), nvec=nvec, avec=avec, bvec=bvec, prior=prior, constraint=constraint)
    w <- c(res$q)
    
    #############RELAXED MEM ADDITION
    ### Create indicator for if simx is in critical range for each method 
    criticalrange_ind <- ifelse(((simx >= criticalpoints$min1)&(simx <= criticalpoints$max1))|
                                  ((simx >= criticalpoints$min2)&(simx <= criticalpoints$max2))|
                                  ((simx >= criticalpoints$min3)&(simx <= criticalpoints$max3)), TRUE, FALSE)
    
    ## If any lie in the critical range, then calculate relaxed weights.
    if(sum(criticalrange_ind,na.rm = TRUE)){
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
      
      allweights <- rbind(w,
                          if(criticalrange_ind[2] %in% TRUE) w2 else(w),
                          if(criticalrange_ind[3] %in% TRUE) w3 else(w),
                          if(criticalrange_ind[4] %in% TRUE) w4 else(w),
                          if(criticalrange_ind[5] %in% TRUE) w5 else(w),
                          if(criticalrange_ind[6] %in% TRUE) w6 else(w),
                          if(criticalrange_ind[7] %in% TRUE) w7 else(w))
      
    } else{
      allweights <- rbind(w,w,w,w,w,w,w) 
    }
    
    for(l in 1:nrow(allweights)){
      
        avec.o <- avec
        bvec.o <- bvec
        xvec.o <- c(simx, xvec[-1])
        nvec.o <- nvec
        mod.mat.o <- res$mod.mat
        allweights.o <- allweights[l,]
        
      
      
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
  
  est.var = sapply(est.mu, var) #estimate variance of estimated mu-hat from simulation
  est.mean = sapply(est.mu, mean) #estimate mu-hat
  est.bias = est.mean - truep #bias
  mse.part = est.bias^2 + est.var #MSE
  mse.sum.part = sapply(est.mu, FUN = function(z) (1/niter)*sum((z - truep)^2)) #MSE sum equation
  mae.part =  sapply(est.mu, FUN = function(z) (1/niter)*sum(abs(z - truep))) #MAE
  est.esss = sapply(est.esss, median) #ESSS
  abserr.part = abs(est.bias) #absolute error of bias for summary stat
  
  ###Create output df
  mat = cbind(truep,est.mean,est.var,est.bias,mse.part,mse.sum.part,mae.part,est.esss,abserr.part,criticalpoints)
  
  df <- as.data.frame(mat)
  df$mu <- xvec[1]/nvec[1]
  df$num_iterations <- niter
  df$model <- rownames(df)
  df$numsources <- numsources
  if(numsources == 3){
    df$diff1 <- xvec[2]/nvec[2] - xvec[1]/nvec[1]
    df$diff2 <- xvec[3]/nvec[3] - xvec[1]/nvec[1]
  }else if(numsources == 2){
    df$diff1 <- xvec[2]/nvec[2] - xvec[1]/nvec[1]
    df$diff2 <- NA
  }

  return(df)
  
}





