#################################################################################
# TITLE: CaseStudies.R
#
# PURPOSE: Run CF case studies (arbitrary k and MSE based)
#
# OUTPUT: case study results as .csv files and
#         figures of MSE vs k for each case study as .jpeg files
#
# SECTIONS: Section 0 - load source file
#           Section 1 - binary example
#           Section 2 - continuous example
#
# AUTHOR: Shannon Thomas
# DATE CREATED: FEB 1, 2026
#################################################################################

##########################################################
############### SECTION 0: LOAD SOURCE FILE ##############
##########################################################


library(xtable)
library(future)
library(future.apply)
library(tidyverse)

source('Functions_RelaxedMEM.R') #set source code file


##########################################################
################ SECTION 1: BINARY EXAMPLE ###############
##########################################################

#########CYSTIC FIBROSIS EXAMPLE 
#https://clinicaltrials.gov/study/NCT05719311?cond=Cystic%20Fibrosis&aggFilters=docs:sap,phase:2%203,results:with,studyType:int&lead=Entero%20Therapeutics&rank=1&tab=results
#https://clinicaltrials.gov/study/NCT04375878?cond=Cystic%20Fibrosis&aggFilters=docs:sap,phase:2%203,results:with,studyType:int&lead=Entero%20Therapeutics&rank=2&tab=results
#https://clinicaltrials.gov/study/NCT03746483?cond=Cystic%20Fibrosis&aggFilters=docs:sap,phase:2%203,results:with,studyType:int&lead=Entero%20Therapeutics&rank=3&tab=results 
#Counts for MS1819 studies (# pts with 1 or more AE)
p.count <- 2
s1.count <- 1
s2.count <- 13

#N's (sample sizes) for MS1819 studies (# pts with 1 or more AE)
p.n.bin <- 13
s1.n.bin <- 14
s2.n.bin <- 40

#store results
res.mat2 <- matrix(nrow=2,ncol=8)
rownames(res.mat2) <- c('Mean','ESSS')
colnames(res.mat2) <- c("No_Borrowing", 
                       "MEM",
                       "relaxedMEM_max",
                       "relaxedMEM_sqrtmax",
                       "relaxedMEM_allornothing",
                       "relaxedMEM_twosource",
                       "relaxedMEM_stretchsqrt",
                       "relaxedMEM_stretch1")

#no borrowing:
res.mat2[,1] <- c(round(p.count/p.n.bin,3),
                 0)

#MEM
pr.use <- "pi_e"

cf_bincp <- findcriticalpoints_binary(xvec = c(p.count,s1.count,s2.count), 
                                      nvec = c(p.n.bin,s1.n.bin,s2.n.bin),
                                      avec=c(1,1,1),bvec = c(1,1,1),
                                      niter = 5000, zeroweightcutoff = 0.01,
                                      MSEtolerance = 0.005)

cf_binex <- relaxed_sim_calc.betabin_critpts(x = c(p.count,s1.count,s2.count), 
                                             nvec = c(p.n.bin,s1.n.bin,s2.n.bin),
                                             avec=c(1,1,1),bvec = c(1,1,1),
                                             niter = 1, RAND = FALSE, zeroweightcutoff = 0.01,
                                             prior = "equal", criticalpoints = cf_bincp)

write.csv(cf_binex, file = "casestudy_bin_5000iter_msetol005.csv", row.names = FALSE)

res.mat2[1,2:8] <- paste0(round(cf_binex$est.mean,3))   #est.mean (est.var) for trt
res.mat2[2,2:8] <- round(cf_binex$est.esss,1)  #ESSS for trt
res.mat2
write.csv(res.mat2, file = "casestudy_bin_5000iter_msetol005_results.csv", row.names = FALSE)

#arbitrary k method:
cf_binex_arbk <- relaxed_sim_calc.betabin(xvec=c(p.count,s1.count,s2.count),
                                          nvec=c(p.n.bin,s1.n.bin,s2.n.bin),
                                          avec=c(1,1,1),bvec = c(1,1,1),
                                          prior="equal",
                                          zeroweightcutoff = 0.01,RAND=FALSE,niter=1,
                                          k=0.82) #no relaxing for k >= 0.85


#run full range of k
kvals <- seq(0.01,0.99,0.01)

plan(multisession)
system.time({
  casestudy_bin_rangek <- Reduce(rbind, 
                                 future_lapply(1:length(kvals), 
                                               FUN = function(j) {relaxed_sim_calc.betabin(x = c(p.count,s1.count,s2.count),
                                                                                           nvec=c(p.n.bin,s1.n.bin,s2.n.bin), 
                                                                                           avec = c(1,1,1), bvec = c(1,1,1),
                                                                                           prior = "equal", constraint = 1,
                                                                                           niter = 20000, k = kvals[j],
                                                                                           zeroweightcutoff = 0.01)
                                               },
                                               future.seed = TRUE
  )
  )
})
plan(sequential)

casestudy_bin_rangek$model <- factor(casestudy_bin_rangek$model, levels = c("MEM","relaxedMEM_twosource",
                                                                            "relaxedMEM_max",
                                                                            "relaxedMEM_stretchsqrt","relaxedMEM_sqrtmax",
                                                                            "relaxedMEM_stretch1","relaxedMEM_allornothing"))

jpeg("Bin_CaseStudy_MSEvsk.jpeg", width = 1600, height = 800, quality = 100)
ggplot(data = casestudy_bin_rangek, 
       aes(x = k, y = mse.part, col = model, id = model)) +
  facet_wrap(~diff1+diff2, nrow = 1)+
  geom_point() + geom_line() +
  ggtitle("Case Study Three Source, Binary Outcome: MSE vs k") +
  theme(text = element_text(size = 15))
dev.off()



##########################################################
############## SECTION 2: CONTINUOUS EXAMPLE #############
##########################################################

#Coefficient of Fat Absorption (CFA) - Continuous example
#only NCT04375878 and NCT03746483 because NCT05719311 only reported prop. > 80%
p.trt.mean <- 66
s1.trt.mean <- 55.6
p.ctrl.mean <- 83.3
s1.ctrl.mean <- 86.2

p.trt.sd <- 19.1
s1.trt.sd <- 21.44
p.ctrl.sd <- 13.1
s1.ctrl.sd <- 7.39

p.trt.n <- 12
p.ctrl.n <- 12
s1.trt.n <- 33
s1.ctrl.n <- 35

#store results - BORROWING IN BOTH ARMS
res.mat3 <- matrix(nrow=5,ncol=8)
rownames(res.mat3) <- c('trt','esss.trt','ctrl','esss.ctrl','trt.effect')
colnames(res.mat3) <- c("No_Borrowing", 
                       "MEM",
                       "relaxedMEM_max",
                       "relaxedMEM_sqrtmax",
                       "relaxedMEM_allornothing",
                       "relaxedMEM_twosource",
                       "relaxedMEM_stretchsqrt",
                       "relaxedMEM_stretch1")

#no borrowing:
res.mat3[,1] <- c(paste0(round(p.trt.mean,3),' (',round(p.trt.sd/sqrt(p.trt.n),3),')'),
                 0,
                 paste0(round(p.ctrl.mean,3),' (',round(p.ctrl.sd/sqrt(p.ctrl.n),3),')'),
                 0, 
                 paste0(round(p.trt.mean-p.ctrl.mean,3),' (',round(sqrt(p.trt.sd^2/p.trt.n + p.ctrl.sd^2/p.ctrl.n),3),')'))

#MSE based method:
set.seed(1)
system.time({
cf_trt_contcp <- findcriticalpoints(x = c(p.trt.mean,s1.trt.mean), 
                                    S = c(p.trt.sd,s1.trt.sd),
                                    N = c(p.trt.n,s1.trt.n), 
                                    niter = 5000,
                                    zeroweightcutoff = 0.01,
                                    stepsize = 1,
                                    MSEtolerance = 5)
})

cf_trt_contex <- relaxed.MEM_sim_calc_critpts(x = c(p.trt.mean,s1.trt.mean), 
                                              S = c(p.trt.sd,s1.trt.sd),
                                              N = c(p.trt.n,s1.trt.n),
                                              prior="pi_e",
                                              niter =  1, RAND = FALSE,
                                              zeroweightcutoff = 0.01,
                                              criticalpoints = cf_trt_contcp)

write.csv(cf_trt_contex, file = "casestudy_cont_trt_5000iter_msetol5.csv", row.names = FALSE)

set.seed(2)
system.time({
cf_ctrl_contcp <- findcriticalpoints(x=c(p.ctrl.mean,s1.ctrl.mean),
                                     S=c(p.ctrl.sd,s1.ctrl.sd),
                                     N=c(p.ctrl.n,s1.ctrl.n), 
                                     niter = 5000,
                                     zeroweightcutoff = 0.01,
                                     stepsize = 1,
                                     MSEtolerance = 5)
})

cf_ctrl_contex <- relaxed.MEM_sim_calc_critpts(x=c(p.ctrl.mean,s1.ctrl.mean),
                                               S=c(p.ctrl.sd,s1.ctrl.sd),
                                               N=c(p.ctrl.n,s1.ctrl.n),
                                               prior="pi_e",
                                               niter = 1, RAND = FALSE,
                                               zeroweightcutoff = 0.01,
                                               criticalpoints = cf_ctrl_contcp)

write.csv(cf_ctrl_contex, file = "casestudy_cont_ctrl_5000iter_msetol5.csv", row.names = FALSE)

res.mat3[1,2:8] <- paste0(round(cf_trt_contex$est.mean,3),' (',round(sqrt(cf_trt_contex$est.var),3),')')   #est.mean (est.var) for trt
res.mat3[2,2:8] <- round(cf_trt_contex$est.esss,1)  #ESSS for trt
res.mat3[3,2:8] <- paste0(round(cf_ctrl_contex$est.mean,3),' (',round(sqrt(cf_ctrl_contex$est.var),3),')') #est.mean (est.var) for ctrl
res.mat3[4,2:8] <- round(cf_ctrl_contex$est.esss,1) #ESSS for ctrl
res.mat3[5,2:8] <- paste0(round(cf_trt_contex$est.mean-cf_ctrl_contex$est.mean,3),' (',round(sqrt(cf_trt_contex$est.var+cf_ctrl_contex$est.var),3),')') #difference in est.mean (trt - ctrl), variance of difference in est.mean (sum)
res.mat3

write.csv(res.mat3, file = "casestudy_cont_5000iter_msetol5_results.csv", row.names = FALSE)

#arbitrary k method:
trt.cf.cont_arbk <- relaxed.MEM_sim_calc(x=c(p.trt.mean,s1.trt.mean),S=c(p.trt.sd,s1.trt.sd),N=c(p.trt.n,s1.trt.n),prior="pi_e",zeroweightcutoff = 0.01,RAND=FALSE,niter=1,
                                        k=0.3) #no relaxing for k >= 0.35
ctrl.cf.cont_arbk <- relaxed.MEM_sim_calc(x=c(p.ctrl.mean,s1.ctrl.mean),S=c(p.ctrl.sd,s1.ctrl.sd),N=c(p.ctrl.n,s1.ctrl.n),prior="pi_e",zeroweightcutoff = 0.01,RAND=FALSE,niter=1,
                                         k=0.75) #no relaxing for k >= 0.8

#run full range of k
plan(multisession)
system.time({
  caststudy_cont_trt_rangek <- Reduce(rbind, 
                                      future_lapply(1:length(kvals), 
                                                FUN = function(j) {relaxed.MEM_sim_calc(x = c(p.trt.mean,s1.trt.mean), 
                                                                                        S = c(p.trt.sd,s1.trt.sd),
                                                                                        N = c(p.trt.n,s1.trt.n),
                                                                                        niter = 20000, k = kvals[j],
                                                                                        zeroweightcutoff = 0.01)
                                                                   }, future.seed = TRUE
                                             )
  )
})
plan(sequential)

caststudy_cont_trt_rangek$model <- factor(caststudy_cont_trt_rangek$model, levels = c("MEM","relaxedMEM_twosource",
                                                                              "relaxedMEM_max",
                                                                              "relaxedMEM_stretchsqrt","relaxedMEM_sqrtmax",
                                                                              "relaxedMEM_stretch1","relaxedMEM_allornothing"))

jpeg("Cont_CaseStudy_TrtArm_MSEvsk.jpeg", width = 1600, height = 800, quality = 100)
ggplot(data = caststudy_cont_trt_rangek[caststudy_cont_trt_rangek$model != "relaxedMEM_twosource",], 
       aes(x = k, y = mse.part, col = model, id = model)) +
  facet_wrap(~diff1)+
  geom_point() + geom_line() +
  ggtitle("Case Study Two Source, Continuous Outcome, Treatment Arm: MSE vs k") +
  theme(text = element_text(size = 15))
dev.off()

plan(multisession)
system.time({
  caststudy_cont_ctrl_rangek <- Reduce(rbind, 
                                       future_lapply(1:length(kvals), 
                                                     FUN = function(j) {relaxed.MEM_sim_calc(x=c(p.ctrl.mean,s1.ctrl.mean),
                                                                                             S=c(p.ctrl.sd,s1.ctrl.sd),
                                                                                             N=c(p.ctrl.n,s1.ctrl.n),
                                                                                             niter = 20000, k = kvals[j],
                                                                                             zeroweightcutoff = 0.01)
                                                                        }, future.seed = TRUE
  )
  )
})
plan(sequential)

caststudy_cont_ctrl_rangek$model <- factor(caststudy_cont_ctrl_rangek$model, levels = c("MEM","relaxedMEM_twosource",
                                                                                      "relaxedMEM_max",
                                                                                      "relaxedMEM_stretchsqrt","relaxedMEM_sqrtmax",
                                                                                      "relaxedMEM_stretch1","relaxedMEM_allornothing"))

jpeg("Cont_CaseStudy_CtrlArm_MSEvsk.jpeg", width = 1600, height = 800, quality = 100)
ggplot(data = caststudy_cont_ctrl_rangek[caststudy_cont_ctrl_rangek$model != "relaxedMEM_twosource",], 
       aes(x = k, y = mse.part, col = model, id = model)) +
  facet_wrap(~diff1)+
  geom_point() + geom_line() +
  ggtitle("Case Study Two Source, Continuous Outcome, Control Arm: MSE vs k") +
  theme(text = element_text(size = 15))
dev.off()