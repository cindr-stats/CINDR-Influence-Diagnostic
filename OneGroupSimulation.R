### Simulations for 


### Basics ####

# Library stuff
#library the packages
library(ggplot2)
library(readr)
library(tibble)
library(dplyr)
library(tidyr)
library(janitor)
library(magrittr)
library(Rtsne)
library(robustHD)
library(GGally)
library(stringr)
#library(plotly)
library(R.utils)
library(vegan)
library(mnormt)
library(MASS)



# Multipermanova function

## Draft version of outlier comparison function agnostic to mapping technique


# preliminary functions:

## Vectorize Dist: Function that takes a dataframe and returns a vectorized version of a distance matrix applied to that data matrix.
## This function is used for generating vectorized data when we step through the Mantel test.

vectorize_dist <- function(df) {
  tmp <- dist(df)
  vec <- c(tmp)
  return(vec)}

## Pull P-Val and Pull F stat return the p-values and f-statstics from the list of output MultiPermanova puts out.

pullPval <- function(x){temp <- x$`Pr(>F)`[1]; return(temp)}

pullFstat <- function(x){temp <- x$`F`[1]; return(temp)}





#generate maps on holdout data (jth iteration is held out)
# this function should be able to operate on multiple mapping techniques or a single one

Generate_Holdout_Maps <- function(data, j, maptype = "MDS", perp = perp_val, group = NULL){
  
  #hold out a single row of the data (in the jth row)
  #note - have to remove, as t-SNE doesn't like NAs
  data_holdout <- data[-j,]
  
  # RE-Standardize and then
  #calculate distances FIRST
  data_dist <- dist(scale(data_holdout))
  
  #### Creates Maps (based on several methods: tSNE and MDS currently)
  
  if(maptype %in% c("MDS", "mds")){
    #generate the classical MDS in 2 dimensions
    map <- cmdscale(data_dist, k = 2)
  }
  else if(maptype == "lda"){
    #generate a linear discriminant analysis map
    
    #check for group response variable - if not present, just print a warning
    if(is.null(group)){return(print("Error, error! No grouping variable selected."))}
    else{
      
      data_holdout$Group <- group[-j] %>% factor() #should add the grouping variable to dataset - with holdout removed
      ld_out <- lda(Group ~ ., data = data_holdout) #equal prior class weights, and Group is the response
      ld_map <- predict(ld_out, data_holdout, dimen = 2) #generate predictions on original data to create map
      map <- ld_map$x[,1:2]
    }
  }
  else if(maptype == "tSNE"){
    #generate the t-SNE map and store x,y in map list
    map <- Rtsne(data_dist, perplexity = perp, pca = pca_options, is_distance = TRUE)$Y
  }
  else if (maptype %in% c("Sammon","sammon")){
    map <- sammon(data_dist, y = cmdscale(data_dist,2), k =2, trace = FALSE)$points #initialized with cmdscale for Sammon
  }
  else if (maptype %in% c("isomds","isoMDS")){
    map <- isoMDS(data_dist, y = cmdscale(data_dist,2), k =2, trace = FALSE)$points #initialized with cmdscale for isoMDS
  }
  else if (maptype %in% c("isomap")){
    data_dist2 <- vegdist(data_holdout, method = "euclidean") # for input to isomap, chooses relatively large neighborhood size as median of distances
    map <- isomap(data_dist2, ndim = 2, epsilon = median(data_dist))$points
  }
  else{print("No maptype selected.")}
  
  
  return(list(map = map))
}




#mapping method supported

MultiPermanova <- function(data, nperms = 1000, perp_val= 10, pca_options = FALSE, maptype = "MDS", lda_group = NULL){
  
  # initialize the size of the data and output lists
  n <- nrow(data)
  maplist <- list()
  
  
  #generate maps on holdout data
  for(j in 1:n){
    
    temp <- Generate_Holdout_Maps(data, j, maptype = maptype, group = lda_group, perp = perp_val)
    maplist[[j]] <- temp$map
    
  }
  
  
  # Insert NAs into the holdout indices
  maplist_new <- list()
  
  #
  for(j in 1:length(maplist)){
    #create temporary vectors with NAs for each of the XY coordingates
    tempX <- R.utils::insert(maplist[[j]][,1], j, NA)
    tempY <- R.utils::insert(maplist[[j]][,2], j, NA)
    
    
    maplist_new[[j]] <- tibble(X = tempX, Y = tempY)
    
  }
  
  #vectorizes the distances and outputs a list of vectorized distance matrices
  dist_list <- lapply(maplist_new, vectorize_dist)
  
  
  ## takes the vectorized distance matrices and smashes them into a single data frame - to build correlation matrix
  
  dist_mat <- dist_list %>% as_tibble(.name_repair = "minimal")
  
  #had to remove the NAs here to get this to run (as it removes when I try to use complete.obs later on)
  #dist_df_na <- sapply(dist_mat, na.omit) %>% as_tibble(.name_repair = "minimal")
  #since I'm using pairwise complete obs, I shouldn't need to do this
  
  cormat <- cor(dist_mat, use = "pairwise.complete.obs", method = "pearson")
  cormat %>% dim()
  
  
  ##generate distance on correlation matrix
  cordist <- sqrt(2*(1-cormat)) %>% as.dist()
  
  model_list <- list()
  
  
  #Step through the permanova method:
  
  for(j in 1: nrow(as.matrix(cordist))){
    #generate the indicator variable for each feature
    indicator_df <- rep(0, nrow(maplist_new[[1]])) %>% as.data.frame()
    
    
    names(indicator_df) <- "Indicator"
    indicator_df$Indicator[j] = 1
    
    ### create the permutation matrix
    #initialize
    permat1 <- matrix(0, nrow = length(indicator_df$Indicator), ncol = length(indicator_df$Indicator))
    permat1[1,] <- 1:length(indicator_df$Indicator)
    
    for(i in 2:nrow(permat1)){
      permat1[i,] <- c(i:length(indicator_df$Indicator), 1:(i-1))
    }
    
    
    #each item in the list is an Adonis test using the indicator for the holdout variable (removing the parallel calculations, with 1 permutation)
    model_list[[j]] <- adonis2(cordist ~ Indicator, data = indicator_df, permutations = 1, parallel = 1)
    
    
    #returns the list of ADONIS outputs
    
    
  }
  
  return(list(modellist = model_list, plot_map = maplist))
}



#### Loop ###############


start_time <- Sys.time()
n_simulations <- 100
pval_vec <- c() #vector to store pvalues
max_F_index_vec <- c()
#pval_mat <- matrix(0, ncol = kmax, nrow = n_simulations)
#generates a different value of covariates at sd of roughly 4

#Set K value (this can be changed to control number of "aberrant dimensions")
#k = 9

for(k in 1:10){

#generates 20 simulations at a given number of predictors
for(i in 1:n_simulations){
  # Simulate one group with 50 observations with a potential influential point:
  g1 <- mnormt::rmnorm(n = 10, mean = rep(20, 50), diag(35, 50)) %>%
    t() %>%
    data.frame()
  #adds a group variable and creates a data frame
  yes_structure <- rbind(g1) %>%
    data.frame() %>%
    mutate(Group = c(rep("1", 50)))
  
  yes_structure_scale <- lapply(yes_structure[,1:10], scale) %>% as_tibble()
  
  ### adds outliers - but not as a single group. Here we have an outlier that includes a really big one, a really small observation, and a very abnormal one
  add_outliers <- rbind(rnorm(10, 0, 1)) %>% as_tibble() #all on normal scale
  add_outliers[1:k] <- 4 #extreme point is outside the 3 standard deviations
  
  
  names(add_outliers) <- names(yes_structure_scale)
  
  yes_structure_final <- bind_rows(yes_structure_scale, add_outliers)
  yes_structure_final$Group <- c(yes_structure$Group, rep("Outlier",1))
  
  
  
  ### Run multipermanova and get the largest F-statistic value
  xMDS <- MultiPermanova(yes_structure_final[,1:10], maptype = "Sammon")
  ### Pulls the F-statistics from each of the holdout repititions
  max_F_inf <- lapply(xMDS$modellist, pullFstat) %>% unlist() %>% max() #return max value
  max_F_index <- lapply(xMDS$modellist, pullFstat) %>% unlist() %>% which.max() #return index of max value
  
  # store max f index for each run (one for each of the "i" runs)
  max_F_index_vec[i] <- max_F_index
  
  # Get the reference distribution (with 50 observations):
  
  n_runs = 99 #try 99 so that we have 99 + 1 = 100 for easy p-value calculations
  
  data_list <- out_ref_list <- list()
  max_F_vec <- rep(0, n_runs)
  
  # simulate with n "normal" observations
  for(j in 1:n_runs){
    
    ### Generate a Uniform Reference Group
    nvars = 10
    nsamps = 30
    #    samps <- replicate(nsamps, runif(10,-3,3)) %>% t() %>% data.frame() %>%
    #    mutate(Group = rep("1",30)) # draw samples on scaled data (all same group to align with other code)
    
    ### Generate a Multivariate Normal reference group (already scaled)
    # n here is number of dimensions
    # nsamps is our number of rows in simulated data set
    refdist <- mnormt::rmnorm(n = nvars, mean = rep(0, nsamps), diag(1, nsamps)) %>%
      t() %>%
      data.frame()
    reference_scale<- refdist %>% data.frame()
    #Scaling step - might not be necessary but for consistency
    #reference_scale <- lapply(samps[,1:10], scale) %>% as_tibble()
    #not necessary if using the multivariate normal (since already 0/1)
    
    
    data_list[[j]] <- reference_final <- reference_scale
    
    
    ### Run the diagnostic and save outputs to the list
    out_mds_ref <- MultiPermanova(reference_scale[,1:10], maptype = "Sammon")
    out_ref_list[[j]] <- out_mds_ref
    
    max_F_vec[j] <- lapply(out_mds_ref$modellist, pullFstat) %>% unlist() %>% max()
    
    
  }
  
  #store the pseudo-F values in a tibble
  max_F_tib <- tibble(Rep = 1:n_runs, MaxF = max_F_vec)
  xupper <- max(max_F_tib$MaxF, max_F_inf) + 5
  
  #create a plot
  # plot_list[[i]] <-
  # ggplot(max_F_tib, aes(MaxF)) + geom_density(fill = "dodgerblue3", alpha = 0.8) + theme_bw() +
  # ggtitle(paste("Iteration ","i", ": Max Pseudo F (Uniform Reference Distribution)")) + geom_vline(xintercept = max_F_inf,size = 2, col = "orange3") +
  # xlim(0, xupper)
  
  #generate an empirical p-value (i.e. observations greater than our observed p-value)
  pval = length(which(max_F_tib$MaxF > max_F_inf))/(n_runs + 1)
  pval_vec[i] <- pval
  
}#for loop (i)




## Store Values
#saveRDS(plot_list, "plot_list.RDS")
write.csv(pval_vec, paste0("OneGroupPvalues_K",k,"_Sammon_1204.csv"))
write.csv(max_F_index_vec, paste0("OneGroupMaxFIndex_K",k,"_Sammon_1204.csv"))

}


## Time check
end_time <- Sys.time()
end_time - start_time







