### Type 1 Error Simulation Study
# This can be used for any map type - just change the maptype from MDS, Sammon, IsoMDS, etc. 

# 3/2025


#library the packages
library(ggplot2)
library(readr)
library(tibble)
library(dplyr)
library(tidyr)
library(janitor)
library(magrittr)
#library(Rtsne)
library(robustHD)
library(GGally)
library(stringr)
#library(plotly)
library(R.utils)
library(vegan)
library(mnormt)
library(MASS)
library(ggpubr)

#source 
source("MultiPermanova Functions.R")



## Simulation Code (Runs 100 Simulations at a time)

plot_list <- list()
data_list <- list()
n_sims <- 100
n_runs <- 99
pval_vec <- rep(0, n_sims)

tic = Sys.time()

for(i in 1:n_sims){
  # Simulate the MVN observation:
  g1 <- mnormt::rmnorm(n = 30, mean = rep(40, 30), diag(50, 30)) %>%
    t() %>%
    data.frame()
  yes_structure <- g1 %>%
    data.frame() %>%
    mutate(Group = c(rep("1", 30)))
  
  yes_structure_scale <- lapply(yes_structure[,1:10], scale) %>% as_tibble()
  yes_structure_scale$Group <- yes_structure$Group
  #save the data to a list and saves data as yes_structure_final
  data_list[[i]] <- yes_structure_final <- yes_structure_scale
  
  ### Run multipermanova and get the largest F-statistic value
  xMDS <- MultiPermanova(yes_structure_final[,1:10], maptype = "isoMDS")
  ### Pulls the F-statistics from each of the holdout repititions
  max_F_inf <- lapply(xMDS$modellist, pullFstat) %>% unlist() %>% max()
  
  #### Get the reference distribution (with 99 observations):
  
  data_list <- out_mds_list <- list()
  max_F_vec <- rep(0, n_runs)
  
  # simulate with 60 "reference" observations
  for(j in 1:n_runs){
    
    ### Generate a Reference Group
    nvars = 10
    nsamps = 60 #size of reference group
    ### Uniform Version
    #samps <- replicate(nsamps, runif(10,-3,3)) %>% t() %>% data.frame() %>%
    #mutate(Group = rep("1",nsamps)) # draw samples on scaled data (all same group to align with other code)
    #Scaling step - might not be necessary but for consistency
    #yes_structure_scale <- lapply(samps[,1:10], scale) %>% as_tibble()
    
    ### Generate a Multivariate Normal reference group (already scaled)
    # n here is number of dimensions
    # nsamps is our number of rows in simulated data set
    refdist <- mnormt::rmnorm(n = nvars, mean = rep(0, nsamps), diag(1, nsamps)) %>%
      t() %>%
      data.frame()
    yes_structure_scale <- refdist %>% data.frame()
    
    
    
    data_list[[j]] <- yes_structure_final <- yes_structure_scale
    ### Run the diagnostic and save outputs to the list
    out_mds <- MultiPermanova(yes_structure_final[,1:10], maptype = "isoMDS")
    out_mds_list[[j]] <- out_mds
    
    max_F_vec[j] <- lapply(out_mds$modellist, pullFstat) %>% unlist() %>% max()
  }
  
  #store the pseudo-F values in a tibble and create a plot
  max_F_tib <- tibble(Rep = 1:n_runs, MaxF = max_F_vec)
  xupper <- max(max_F_tib$MaxF, max_F_inf) + 5  #for plot xlimit
  
  #generate an empirical p-value (i.e. observations greater than our observed p-value)
  pval = length(which(max_F_tib$MaxF > max_F_inf))/(n_runs + 1)
  pval_vec[i] <- pval
  
  
}#for loop (i)

toc = Sys.time()
toc - tic

pval_df <- tibble(Iteration = 1:n_sims, Pval = pval_vec)
write.csv(pval_vec, "IsoMDSpvalues.csv")

#length(which(pval_df$Pval < 0.05))



