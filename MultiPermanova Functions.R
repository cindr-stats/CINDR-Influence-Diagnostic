## Multipermanova Function
## 
## 


## Helper Functions 

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





## Main Function 
MultiPermanova <- function(data, nperms = 1000, perp_val= 10, pca_options = FALSE, maptype = "MDS", lda_group = NULL){
  
  # initialize the size of the data and output lists
  n <- nrow(data)
  maplist <- list()
  
  
  #generate maps on holdout data
  for(j in 1:n){
    
    temp <- Generate_Holdout_Maps(data, j, maptype = maptype, group = lda_group)
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
  #cormat %>% dim()
  
  ## Use the lower triangle 
  cormat_lower <- lower.tri(cormat)
  
  
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
    
    
    #each item in the list is an Adonis test using the indicator for the holdout variable
    
    model_list[[j]] <- adonis2(cordist ~ Indicator, data = indicator_df, permutations = permat1, parallel = 1)
    
    
    #returns the list of ADONIS outputs
    
    
  }
  
  return(list(modellist = model_list, plot_map = maplist))
}
