#################################################################################
# TITLE: Figures_MSEbased.R
#
# PURPOSE: Create figures for all relaxed MEM simulations. 
#
# OUTPUT: all MSE based figures as .jpeg files
#
# SECTIONS: Section 0 - load simulation results data set
#           Section 1 - Continuous MEM Figures
#           Section 2 - Binary MEM Figures 
#
# AUTHOR: Shannon Thomas
# DATE CREATED: DEC 16, 2025
#################################################################################

##########################################################
################### SECTION 0: LOAD DATA #################
##########################################################

library(tidyverse)
library(readr)
allresults <- read_csv("MSEbased_results.csv")
allresults$model <- factor(allresults$model, levels = c("MEM","relaxedMEM_twosource",
                                                "relaxedMEM_max",
                                                "relaxedMEM_stretchsqrt","relaxedMEM_sqrtmax",
                                                "relaxedMEM_stretch1","relaxedMEM_allornothing"))

#get relative values
allresults <- allresults %>% group_by(numsources, outcome_type, diff1_c, diff2_c) %>%
  mutate(absbias_diff = abs(est.bias) - abs(est.bias[1]),
         absbias_diff = replace(absbias_diff, row_number() == 1, NA),
         var_diff = est.var - est.var[1],
         var_diff = replace(var_diff, row_number() == 1, NA),
         esss_diff = est.esss - est.esss[1],
         esss_diff = replace(esss_diff, row_number() == 1, NA),
         mse_diff = mse.part - mse.part[1],
         mse_diff = replace(mse_diff, row_number() == 1, NA)) 
allresults$scenario <- interaction(allresults$diff1_c, allresults$diff2_c)

#define subset of models for paper figures
modelsub <- c("relaxedMEM_max",
              "relaxedMEM_stretch1","relaxedMEM_allornothing")

##########################################################
################ SECTION 1: CONTINUOUS MEM ###############
##########################################################

## THREE SOURCE
jpeg("MSEbased_Figures/MSEBased_Diff_Cont_3S_biasvsvar.jpeg", width = 1600, height = 800, quality = 100)
ggplot(data = allresults[(allresults$numsources == 3) & (allresults$outcome_type == "continuous") &
                           (allresults$model %in% modelsub),], 
       aes(y = absbias_diff, x = var_diff, col = model)) +
  geom_point(size = 4) +
  facet_wrap(~diff1_c + diff2_c, nrow = 1) +
  ggtitle("Three Source, Continuous Outcome: Difference in Absolute Bias vs Difference in Variance Relative to Traditional MEM") +
  theme(text = element_text(size = 15)) +
  xlab("Difference in Variance") + ylab("Difference in Absolute Bias")
dev.off()

jpeg("MSEbased_Figures/MSEBased_Cont_3S_biasvsvar.jpeg", width = 1600, height = 800, quality = 100)
ggplot(data = allresults[(allresults$numsources == 3) & (allresults$outcome_type == "continuous"),], 
       aes(y = est.bias, x = est.var, col = model)) +
  geom_point(size = 4) +
  facet_wrap(~diff1_c + diff2_c, nrow = 1) +
  ggtitle("Three Source, Continuous Outcome: Bias vs Variance") +
  theme(text = element_text(size = 15))
dev.off()

jpeg("MSEbased_Figures/MSEBased_Diff_Cont_3S_biasvsESSS.jpeg", width = 1600, height = 800, quality = 100)
ggplot(data = allresults[(allresults$numsources == 3) & (allresults$outcome_type == "continuous") &
                           (allresults$model %in% modelsub),], 
       aes(y = absbias_diff, x = esss_diff, col = model)) +
  geom_point(size = 4) +
  facet_wrap(~diff1_c + diff2_c, nrow = 1) +
  ggtitle("Three Source, Continuous Outcome: Difference in Absolute Bias vs Difference in ESSS Relative to Traditional MEM") +
  theme(text = element_text(size = 15)) +
  xlab("Difference in ESSS") + ylab("Difference in Absolute Bias")
dev.off()

jpeg("MSEbased_Figures/MSEBased_Cont_3S_biasvsESSS.jpeg", width = 1600, height = 800, quality = 100)
ggplot(data = allresults[(allresults$numsources == 3) & (allresults$outcome_type == "continuous"),], 
       aes(y = est.bias, x = est.esss, col = model)) +
  geom_point(size = 4) +
  facet_wrap(~diff1_c + diff2_c, nrow = 1) +
  ggtitle("Three Source, Continuous Outcome: Bias vs ESSS") +
  theme(text = element_text(size = 15))
dev.off()

jpeg("MSEbased_Figures/MSEBased_Diff_Cont_3S_MSEvsESSS.jpeg", width = 1600, height = 800, quality = 100)
ggplot(data = allresults[(allresults$numsources == 3) & (allresults$outcome_type == "continuous") &
                           (allresults$model %in% modelsub),], 
       aes(y = mse_diff, x = esss_diff, col = model)) +
  geom_point(size = 4) +
  facet_wrap(~diff1_c + diff2_c, nrow = 1) +
  ggtitle("Three Source, Continuous Outcome: Difference in MSE vs Difference in ESSS Relative to Traditional MEM") +
  theme(text = element_text(size = 15)) +
  xlab("Difference in ESSS") + ylab("Difference in MSE")
dev.off()

jpeg("MSEbased_Figures/MSEBased_Cont_3S_MSEvsESSS.jpeg", width = 1600, height = 800, quality = 100)
ggplot(data = allresults[(allresults$numsources == 3) & (allresults$outcome_type == "continuous"),], 
       aes(y = mse.part, x = est.esss, col = model)) +
  geom_point(size = 4) +
  facet_wrap(~diff1_c + diff2_c, nrow = 1) +
  ggtitle("Three Source, Continuous Outcome: MSE vs ESSS") +
  theme(text = element_text(size = 15))
dev.off()

## TWO SOURCE

jpeg("MSEbased_Figures/MSEBased_Diff_Cont_2S_biasvsvar.jpeg", width = 1600, height = 800, quality = 100)
ggplot(data = allresults[(allresults$numsources == 2) & (allresults$outcome_type == "continuous") &
                           (allresults$model %in% modelsub) & (allresults$model != "relaxedMEM_twosource"),], 
       aes(y = absbias_diff, x = var_diff, col = model)) +
  geom_point(size = 4) +
  facet_wrap(~diff1_c, nrow = 1) +
  ggtitle("Two Source, Continuous Outcome: Difference in Absolute Bias vs Difference in Variance Relative to Traditional MEM") +
  theme(text = element_text(size = 15)) +
  xlab("Difference in Variance") + ylab("Difference in Absolute Bias")
dev.off()

jpeg("MSEbased_Figures/MSEBased_Cont_2S_biasvsvar.jpeg", width = 1600, height = 800, quality = 100)
ggplot(data = allresults[(allresults$numsources == 2) & 
                           (allresults$outcome_type == "continuous") & 
                           (allresults$model != "relaxedMEM_twosource"),], 
       aes(y = est.bias, x = est.var, col = model)) +
  geom_point(size = 4) +
  facet_wrap(~diff1_c, nrow = 1) +
  ggtitle("Two Source, Continuous Outcome: Bias vs Variance") +
  theme(text = element_text(size = 15))
dev.off()

jpeg("MSEbased_Figures/MSEBased_Diff_Cont_2S_biasvsESSS.jpeg", width = 1600, height = 800, quality = 100)
ggplot(data = allresults[(allresults$numsources == 2) & (allresults$outcome_type == "continuous") &
                           (allresults$model %in% modelsub) & (allresults$model != "relaxedMEM_twosource"),], 
       aes(y = absbias_diff, x = esss_diff, col = model)) +
  geom_point(size = 4) +
  facet_wrap(~diff1_c, nrow = 1) +
  ggtitle("Two Source, Continuous Outcome: Difference in Absolute Bias vs Difference in ESSS Relative to Traditional MEM") +
  theme(text = element_text(size = 15)) +
  xlab("Difference in ESSS") + ylab("Difference in Absolute Bias")
dev.off()

jpeg("MSEbased_Figures/MSEBased_Cont_2S_biasvsESSS.jpeg", width = 1600, height = 800, quality = 100)
ggplot(data = allresults[(allresults$numsources == 2) & (allresults$outcome_type == "continuous") &
                           (allresults$model != "relaxedMEM_twosource"),], 
       aes(y = est.bias, x = est.esss, col = model)) +
  geom_point(size = 4) +
  facet_wrap(~diff1_c, nrow = 1)  +
  ggtitle("Two Source, Continuous Outcome: Bias vs ESSS") +
  theme(text = element_text(size = 15))
dev.off()

jpeg("MSEbased_Figures/MSEBased_Diff_Cont_2S_MSEvsESSS.jpeg", width = 1600, height = 800, quality = 100)
ggplot(data = allresults[(allresults$numsources == 2) & (allresults$outcome_type == "continuous") &
                           (allresults$model %in% modelsub) & (allresults$model != "relaxedMEM_twosource"),], 
       aes(y = mse_diff, x = esss_diff, col = model)) +
  geom_point(size = 4) +
  facet_wrap(~diff1_c, nrow = 1) +
  ggtitle("Two Source, Continuous Outcome: Difference in MSE vs Difference in ESSS Relative to Traditional MEM") +
  theme(text = element_text(size = 15)) +
  xlab("Difference in ESSS") + ylab("Difference in MSE")
dev.off()

jpeg("MSEbased_Figures/MSEBased_Cont_2S_MSEvsESSS.jpeg", width = 1600, height = 800, quality = 100)
ggplot(data = allresults[(allresults$numsources == 2) & (allresults$outcome_type == "continuous") &
                           (allresults$model != "relaxedMEM_twosource"),], 
       aes(y = mse.part, x = est.esss, col = model)) +
  geom_point(size = 4) +
  facet_wrap(~diff1_c, nrow = 1) +
  ggtitle("Two Source, Continuous Outcome: MSE vs ESSS") +
  theme(text = element_text(size = 15))
dev.off()


##########################################################
################## SECTION 2: BINARY MEM #################
##########################################################

## THREE SOURCE

jpeg("MSEbased_Figures/MSEBased_Diff_Bin_3S_biasvsvar.jpeg", width = 1600, height = 800, quality = 100)
ggplot(data = allresults[(allresults$numsources == 3) & (allresults$outcome_type == "binary") &
                           (allresults$model %in% modelsub),], 
       aes(y = absbias_diff, x = var_diff, col = model)) +
  geom_point(size = 4) +
  facet_wrap(~diff1_c + diff2_c, nrow = 1) +
  ggtitle("Three Source, Binary Outcome: Difference in Absolute Bias vs Difference in Variance Relative to Traditional MEM") +
  theme(text = element_text(size = 15)) +
  xlab("Difference in Variance") + ylab("Difference in Absolute Bias")
dev.off()

jpeg("MSEbased_Figures/MSEBased_Bin_3S_biasvsvar.jpeg", width = 1600, height = 800, quality = 100)
ggplot(data = allresults[(allresults$numsources == 3) & (allresults$outcome_type == "binary"),], 
       aes(y = est.bias, x = est.var, col = model)) +
  geom_point(size = 4) +
  facet_wrap(~diff1_c + diff2_c, nrow = 1) +
  ggtitle("Three Source, Binary Outcome: Bias vs ESSS") +
  theme(text = element_text(size = 15))
dev.off()

jpeg("MSEbased_Figures/MSEBased_Diff_Bin_3S_biasvsESSS.jpeg", width = 1600, height = 800, quality = 100)
ggplot(data = allresults[(allresults$numsources == 3) & (allresults$outcome_type == "binary") &
                           (allresults$model %in% modelsub),], 
       aes(y = absbias_diff, x = esss_diff, col = model)) +
  geom_point(size = 4) +
  facet_wrap(~diff1_c + diff2_c, nrow = 1) +
  ggtitle("Three Source, Binary Outcome: Difference in Absolute Bias vs Difference in ESSS Relative to Traditional MEM") +
  theme(text = element_text(size = 15)) +
  xlab("Difference in ESSS") + ylab("Difference in Absolute Bias")
dev.off()

jpeg("MSEbased_Figures/MSEBased_Bin_3S_biasvsESSS.jpeg", width = 1600, height = 800, quality = 100)
ggplot(data = allresults[(allresults$numsources == 3) & (allresults$outcome_type == "binary"),], 
       aes(y = est.bias, x = est.esss, col = model)) +
  geom_point(size = 4) +
  facet_wrap(~diff1_c + diff2_c, nrow = 1) +
  ggtitle("Three Source, Binary Outcome: Bias vs ESSS") +
  theme(text = element_text(size = 15))
dev.off()

jpeg("MSEbased_Figures/MSEBased_Diff_Bin_3S_MSEvsESSS.jpeg", width = 1600, height = 800, quality = 100)
ggplot(data = allresults[(allresults$numsources == 3) & (allresults$outcome_type == "binary") &
                           (allresults$model %in% modelsub),], 
       aes(y = mse_diff, x = esss_diff, col = model)) +
  geom_point(size = 4) +
  facet_wrap(~diff1_c + diff2_c, nrow = 1) +
  ggtitle("Three Source, Binary Outcome: Difference in MSE vs Difference in ESSS Relative to Traditional MEM") +
  theme(text = element_text(size = 15)) +
  xlab("Difference in ESSS") + ylab("Difference in MSE")
dev.off()

jpeg("MSEbased_Figures/MSEBased_Bin_3S_MSEvsESSS.jpeg", width = 1600, height = 800, quality = 100)
ggplot(data = allresults[(allresults$numsources == 3) & (allresults$outcome_type == "binary"),], 
       aes(y = mse.part, x = est.esss, col = model)) +
  geom_point(size = 4) +
  facet_wrap(~diff1_c + diff2_c, nrow = 1) +
  ggtitle("Three Source, Binary Outcome: MSE vs ESSS") +
  theme(text = element_text(size = 15))
dev.off()

## TWO SOURCE

jpeg("MSEbased_Figures/MSEBased_Diff_Bin_2S_biasvsvar.jpeg", width = 1600, height = 800, quality = 100)
ggplot(data = allresults[(allresults$numsources == 2) & (allresults$outcome_type == "binary") &
                           (allresults$model %in% modelsub) & (allresults$model != "relaxedMEM_twosource"),], 
       aes(y = absbias_diff, x = var_diff, col = model)) +
  geom_point(size = 4) +
  facet_wrap(~diff1_c, nrow = 1) +
  ggtitle("Two Source, Binary Outcome: Difference in Absolute Bias vs Difference in Variance Relative to Traditional MEM") +
  theme(text = element_text(size = 15)) +
  xlab("Difference in Variance") + ylab("Difference in Absolute Bias")
dev.off()

jpeg("MSEbased_Figures/MSEBased_Bin_2S_biasvsvar.jpeg", width = 1600, height = 800, quality = 100)
ggplot(data = allresults[(allresults$numsources == 2) & (allresults$outcome_type == "binary") &
                           (allresults$model != "relaxedMEM_twosource"),], 
       aes(y = est.bias, x = est.var, col = model)) +
  geom_point(size = 4) +
  facet_wrap(~diff1_c, nrow = 1) +
  ggtitle("Two Source, Binary Outcome: Bias vs ESSS") +
  theme(text = element_text(size = 15))
dev.off()

jpeg("MSEbased_Figures/MSEBased_Diff_Bin_2S_biasvsESSS.jpeg", width = 1600, height = 800, quality = 100)
ggplot(data = allresults[(allresults$numsources == 2) & (allresults$outcome_type == "binary") &
                           (allresults$model %in% modelsub) & (allresults$model != "relaxedMEM_twosource"),], 
       aes(y = absbias_diff, x = esss_diff, col = model)) +
  geom_point(size = 4) +
  facet_wrap(~diff1_c, nrow = 1) +
  ggtitle("Two Source, Binary Outcome: Difference in Absolute Bias vs Difference in ESSS Relative to Traditional MEM") +
  theme(text = element_text(size = 15)) +
  xlab("Difference in ESSS") + ylab("Difference in Absolute Bias")
dev.off()

jpeg("MSEbased_Figures/MSEBased_Bin_2S_biasvsESSS.jpeg", width = 1600, height = 800, quality = 100)
ggplot(data = allresults[(allresults$numsources == 2) & (allresults$outcome_type == "binary") &
                           (allresults$model != "relaxedMEM_twosource"),], 
       aes(y = est.bias, x = est.esss, col = model)) +
  geom_point(size = 4) +
  facet_wrap(~diff1_c, nrow = 1) +
  ggtitle("Two Source, Binary Outcome: Bias vs ESSS") +
  theme(text = element_text(size = 15))
dev.off()

jpeg("MSEbased_Figures/MSEBased_Diff_Bin_2S_MSEvsESSS.jpeg", width = 1600, height = 800, quality = 100)
ggplot(data = allresults[(allresults$numsources == 2) & (allresults$outcome_type == "binary") &
                           (allresults$model %in% modelsub) & (allresults$model != "relaxedMEM_twosource"),], 
       aes(y = mse_diff, x = esss_diff, col = model)) +
  geom_point(size = 4) +
  facet_wrap(~diff1_c, nrow = 1) +
  ggtitle("Two Source, Binary Outcome: Difference in MSE vs Difference in ESSS Relative to Traditional MEM") +
  theme(text = element_text(size = 15)) +
  xlab("Difference in ESSS") + ylab("Difference in MSE")
dev.off()

jpeg("MSEbased_Figures/MSEBased_Bin_2S_MSEvsESSS.jpeg", width = 1600, height = 800, quality = 100)
ggplot(data = allresults[(allresults$numsources == 2) & (allresults$outcome_type == "binary") &
                           (allresults$model != "relaxedMEM_twosource"),], 
       aes(y = mse.part, x = est.esss, col = model)) +
  geom_point(size = 4) +
  facet_wrap(~diff1_c, nrow = 1) +
  ggtitle("Two Source, Binary Outcome: MSE vs ESSS") +
  theme(text = element_text(size = 15))
dev.off()



