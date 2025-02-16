#### Reference Distribution Function
## Note: This function is intended to be called after a call to MultiPermanova() and uses the outputs. 
## 
## Inputs: 
#   + dat_scale: Scaled Data input to CINDR
#   + max_F_inf: Max pseudo F of the potentially influential observation (from CINDR output)
#   + n_runs: Number of permutation runs for reference distribution
#   + dist_type: Specifies the reference distribution type. Currently either a Multivariate Normal or Uniform reference distribution are available. 
#   + 
#
####


CINDR_reference <- function(dat_scale, max_F_inf, n_runs= 99, dist_type = "MVN", map_type = "MDS", sig_level = 0.05){
  
  #
  ## Generate Reference Distribution
  #library(mnormt)
  #library(MASS)
  
 #Default n_runs is 99 so that we have 99 + 1 = 100 for easy p-value calculations
  
  data_list <- out_ref_list <- list()
  max_F_vec <- rep(0, n_runs)
  nvars = ncol(dat_scale)
  nsamps = nrow(dat_scale)
  
  # simulate with n "normal" observations
  for(j in 1:n_runs){
    
    
    
    if(dist_type %in% c("Uniform","U","u","uniform")){
      ### Generate a Uniform Reference Group
      samps <- replicate(nsamps, runif(n,-3,3)) %>% t() %>% data.frame() %>%
      mutate(Group = rep("1",nsamps)) # draw samples on scaled data (all same group to align with other code)
      }
    else{
    ### Generate a Multivariate Normal reference group (already scaled)
    # n here is number of dimensions (based on dat_scale)
    # nsamps is our number of rows in simulated data set (based on dat_scale)
    refdist <- mnormt::rmnorm(n = nvars, mean = rep(0, nsamps), diag(1, nsamps)) %>%
      t() %>%
      data.frame()
    reference_scale <- refdist %>% data.frame()
    }
    
    #Store scaled reference data in list
    data_list[[j]] <- reference_final <- reference_scale
    
    ### Run the diagnostic and save outputs to the list
    out_mds_ref <- MultiPermanova(reference_scale, maptype = map_type)
    out_ref_list[[j]] <- out_mds_ref
    
    #pulls index of the largest pseudo-F statistic
    max_F_vec[j] <- lapply(out_mds_ref$modellist, pullFstat) %>% unlist() %>% max()
    
  }
  
  #store the pseudo-F values in a tibble
  max_F_tib <- tibble(Rep = 1:n_runs, MaxF = max_F_vec)
  xupper <- max(max_F_tib$MaxF, max_F_inf) + 5
  
  #create a plot of reference distribution
  refplot <- ggplot(max_F_tib, aes(MaxF)) + geom_density(fill = "dodgerblue3", alpha = 0.8) + 
    ggtitle(paste("Max Pseudo F Statistics")) + geom_vline(xintercept = max_F_inf,size = 2, col = "orange3") + xlim(0, xupper)
  
  #generate an empirical p-value (i.e. observations greater than our observed p-value)
  pval <-  length(which(max_F_tib$MaxF > max_F_inf))/(n_runs + 1)
  
  #guidance 
  guidance <- ifelse(pval <= sig_level, c("There is evidence that this observation is influential."), c("There is a lack of evidence that this observation is influential."))
  
  
  #Function returns
  return(list(refplot = refplot, guidance = guidance, p_value = pval, Fstats = max_F_tib, Fstat_inf = max_F_inf))
  
  
}



## Testing Code: 


