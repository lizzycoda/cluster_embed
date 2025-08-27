library(torch)
source("./scripts/utils.R")


# applies loe to each cluster 
clusterwise_loe <- function(clusters, distances){
  clusters <- reindex_clusters(clusters)
  distances <- as.matrix(distances)
  n = nrow(distances)
  
  # get an loe embedding for each cluster
  embeddings <- matrix(NA, n,2)
  for (i in 1:max(clusters)){
    idxs = which(clusters == i)
    dist_i = distances[idxs,idxs]
    m = length(idxs)
    X_init = matrix(rnorm(m*2), nrow = m)
    emb_i = get_loe(X_init, dist_i, 10) 
    embeddings[idxs,] = emb_i$X_embedded
  }
  
  #scale each cluster embedding
  for(i in 1:max(clusters)){
    idxs = which(clusters == i)
    emb_distances = dist(embeddings[idxs,])
    dist_i <- as.dist(distances[idxs,idxs])
    scale_factor = sum(emb_distances*dist_i)/sum(emb_distances^2)
    embeddings[idxs,] = embeddings[idxs,]*scale_factor
    
  }
  
  return(embeddings)
  
}


get_loe <-function(X_init, distances, k, max_iter = 250, delta = 2e-6, thresh = 10e-6){
  distances <- as.matrix(distances)
  n = dim(X_init)[1]
  
  #compute triplets 
  all_sources  = c()
  all_neighbors = c()
  all_nonneighbors = c()
  
  for(i in 1:n){
    
    #remove sort later, use it now to debug
    t = Sys.time()
    distance_ordering <- order(distances[i,])
    neighbors = distance_ordering[2:(k+1)]
    non_neighbors = distance_ordering[(k+2):n]
    e = Sys.time()
    
    
    #is this part the bottleneck?
    all_sources = c(all_sources, rep(i, k*(n-k-1)))
    all_neighbors = c(all_neighbors, rep(neighbors, each = n-k-1))
    all_nonneighbors = c(all_nonneighbors, rep(non_neighbors, k))
    
  }
  
  #initialize embedding as a tensor
  X_embedded <- torch_tensor(X_init, requires_grad = TRUE)
  
  optim <- optim_adam(X_embedded, lr = 1)
  
  #no batches seems to be fine 
  losses = c()
  for(i in 1:max_iter){
    X_i <- X_embedded[all_sources,] 
    X_j <- X_embedded[all_neighbors,] 
    X_k <- X_embedded[all_nonneighbors,]
    
    loss <- loe_loss(X_i, X_j, X_k, delta)
    losses = c(losses,loss$item())
    
    if(loss$item() < thresh){
      break 
    }
    
    optim$zero_grad()
    loss$backward()
    optim$step()
    
  } 
  
  return(list(X_embedded = as.array(X_embedded), loss = losses))
  
}


loe_loss <-function(x_i, x_j, x_k, delta){
  loss = delta + torch_norm(x_i - x_j, p =2, dim = 2) - torch_norm(x_i - x_k, p=2,  dim = 2)
  loss = loss[loss > 0]
  return(torch_sum(loss))
}