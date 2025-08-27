library(ggplot2)
library(readr)
library(dplyr)
library(Rtsne)
library(stats)
library(vegan)
library(FreeSortR)
library(igraph)
library(dbscan)

#--------------------------------------------------------------------------
# Embedding functions  ----------------------------------------------------
#--------------------------------------------------------------------------

save_embeddings <- function(embeddings, gt_labels, path, cluster_labels = NULL, center = T) {
  if (center) {
    embeddings <- scale(embeddings, scale = F) # center the embeddings
  }
  df <- as_tibble(embeddings)
  colnames(df) <- c("V1", "V2")
  df$labels <- as.factor(gt_labels)
  if (!is.null(cluster_labels)) {
    df$clusters <- as.factor(cluster_labels)
  }
  write_csv(df, path)
}

get_tsne <- function(X, perplexity) {
  out <- Rtsne(X, dims = 2, perplexity = 20, pca = TRUE, check_duplicates = FALSE)$Y
  return(out)
}

get_pca <- function(euclidean_distances) {
  return(stats::cmdscale(distances))
}

get_isomap <- function(distances, dist_type = c("geodesic", "euclidean"), k = 10) {
  dist_type <- match.arg(dist_type)

  if (dist_type == "geodesic") {
    return(stats::cmdscale(distances))
  } else {
    return(vegan::isomap(distances, k = k, ndim = 2)$points)
  }
}


#--------------------------------------------------------------------------
# Clustering functions ----------------------------------------------------
#--------------------------------------------------------------------------

load_clusters <- function(path) {
  clusters <- read_csv(path)
  clusters <- pull(clusters, 1)
  return(clusters)
}

save_clusters <- function(clusters, path) {
  df <- as_tibble(clusters)
  write_csv(df, path)
}

# wrapper function for rand_index
rand_index <- function(labels1, labels2) {
  rand <- FreeSortR::RandIndex(labels1, labels2)
  return(rand$Rand)
}


# reindex cluster labels so labels are 1:num_clusters
reindex_clusters <- function(clusters) {
  vals <- unique(clusters)
  n_clusters <- length(vals)

  new_labels <- clusters
  for (i in 1:n_clusters) {
    idxs <- which(clusters == vals[i])
    new_labels[idxs] <- i
  }

  return(new_labels)
}

# takes in labels (taken as ground truth) and cluster labels from a clustering function to get confusion matrix
plot_confusion_matrix <- function(labels, clusters) {
  classes <- unique(labels)
  classes_labels <- reindex_clusters(labels) # make take on 1:num_labels

  clusters <- reindex_clusters(clusters) # make take on 1:num_clusters
  cluster_labels <- sort(unique(clusters))

  n <- length(classes)
  m <- length(cluster_labels)

  confusion <- matrix(0, nrow = n, ncol = m)

  for (j in 1:m) {
    counts <- table(classes_labels[which(clusters == j)])
    confusion[as.numeric(names(counts)), j] <- as.vector(counts) # already sorted
  }

  actual <- rep(classes, each = m) # want to take original names
  predicted <- rep(cluster_labels, n)
  counts <- as.vector(t(confusion))

  confusion_df <- as_tibble(cbind(actual, predicted, counts))
  confusion_df$predicted <- as.factor(predicted)
  confusion_df$actual <- as.factor(actual)
  confusion_df <- confusion_df |> arrange(actual)


  p <- ggplot(confusion_df, aes(actual, predicted, fill = counts)) +
    geom_tile() +
    geom_text(aes(label = counts)) +
    scale_fill_gradient(low = "white", high = "cyan3") +
    labs(x = "Reference", y = "Prediction")

  return(p)
}


# cluster the data with the Leiden algorithm
get_leiden_clustering <- function(distances, k, resolution) {
  n <- nrow(distances)
  nn <- kNN(distances, k)$id
  edges <- c()
  for (i in 1:n) {
    targets <- nn[i, ][which(nn[i, ] > i)] # convention is edge from source < target
    if (length(targets) > 0) {
      sources <- rep(i, length(targets))
      edges <- c(edges, c(rbind(sources, targets)))
    }
  }
  g <- make_graph(n = n, edges = edges, directed = F)


  ldc <- cluster_leiden(g, n_iterations = 100, objective_function = "modularity", resolution = resolution)
  clusters <- membership(ldc)

  return(clusters)
}


#--------------------------------------------------------------------------
# Geodesic distance calculation -------------------------------------------
#--------------------------------------------------------------------------

calculate_geodesic_distances <- function(euclidean_distances, k, save = TRUE, path = NULL) {
  geodesic_distances <- vegan::isomapdist(euclidean_distances, k = k)
  if (save) {
    df <- as_tibble(as.vector(geodesic_distances))
    write_csv(df, path)
  }
  return(geodesic_distances)
}


# load as a dist object
load_geodesic_distances <- function(path, n) {
  geodesic_distances <- read_csv(path)
  geodesic_distances <- pull(geodesic_distances, 1)
  
  # cast as a dist object
  geodesic_distances <- convert_to_dist(geodesic_distances, n)
  
  return(geodesic_distances)
}


# convert a vector of length n*(n-1)/2 to a dist object
convert_to_dist <- function(distance_vec, n, method = "lower") {
  m <- matrix(0, n, n)
  
  if (method == "lower") {
    lower_tri_indices <- which(lower.tri(m, diag = FALSE))
    m[lower_tri_indices] <- distance_vec
  } else if (method == "upper") {
    upper_tri_indices <- which(upper.tri(m, diag = FALSE))
    m[upper_tri_indices] <- distance_vec
  }
  m <- m + t(m) # make symmetric
  return(as.dist(m))
}
