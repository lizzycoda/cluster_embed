library(ggplot2)
library(tidyr)
library(dbscan)

# -------------------------------------------------------------------------
# Evaluation functions ----------------------------------------------------
# -------------------------------------------------------------------------

# calculates the spearman correlation of the original distances and embedded distances
# inputs should be dist objects
get_spearman <- function(orig_distances, emb_distances) {
  orig_distances <- as.dist(orig_distances)
  emb_distances <- as.dist(emb_distances)

  return(cor(orig_distances, emb_distances, method = "spearman"))
}

# OLD: calculates spearman correlation for a random sample of n distances
# get_spearman <- function(orig_distances, emb_distances, n = 5000, seed = NULL){
#   if(!is.null(seed)){
#     set.seed(seed)
#   }
#
#   #sample n idxs, each corresponding to distance between a pair of points
#   idxs = sample(length(orig_distances), size = n)
#
#   #calculate the spearman correlation between the distances
#   return(cor(orig_distances[idxs], emb_distances[idxs], method = 'spearman'))
#
# }


# calculates the overall stress and normalized stress
get_stress <- function(orig_distances, emb_distances) {
  orig_distances <- as.dist(orig_distances)
  emb_distances <- as.dist(emb_distances)

  stress <- sum(orig_distances - emb_distances)^2
  normalized_stress <- stress / sum(orig_distances^2)
  return(list(
    stress = stress,
    normalized_stress = sqrt(normalized_stress)
  ))
}



# for each eps in eps_vals:
# returns the precision and recall of an eps nbhd (averaged over all points)
# precision is the proportion of points in embedded eps nbhd also in original eps nbhd
# recall is the proportion of points in the original eps nbhd also in embedded eps nbhd
#
# This function enables comparing the precision and recall of different embedding methods
# It takes in the original distances, a list of the embedding distances, and method labels
# It returns a precision data frame and a recall data frame
get_eps_precision_recall <- function(orig_distances, emb_distances_list, method_labels, eps_vals = NULL) {
  # if no eps values provided, step over .05 quantiles of the original distances
  if (is.null(eps_vals)) {
    eps_vals <- quantile(orig_distances, probs = seq(0, 1, .05))
  }

  # reformat dist objects
  orig_distances <- as.matrix(orig_distances)

  m <- length(emb_distances_list)
  for (i in 1:m) {
    emb_distances_list[[i]] <- as.matrix(emb_distances_list[[i]])
  }


  precision <- matrix(NA, ncol = m, nrow = length(eps_vals))
  recall <- matrix(NA, ncol = m, nrow = length(eps_vals))


  for (i in 1:length(eps_vals)) {
    eps <- eps_vals[i]
    orig_nbhd <- orig_distances < eps
    orig_counts <- rowSums(orig_nbhd)

    for (j in 1:m) {
      emb_distances <- emb_distances_list[[j]]
      emb_nbhd <- emb_distances < eps
      emb_counts <- rowSums(emb_nbhd)

      overlap <- rowSums(orig_nbhd & emb_nbhd)
      precision[i, j] <- mean(overlap / emb_counts)
      recall[i, j] <- mean(overlap / orig_counts)
    }
  }

  precision_df <- as_tibble(precision)
  colnames(precision_df) <- method_labels
  precision_df$eps <- eps_vals

  recall_df <- as_tibble(recall)
  colnames(recall_df) <- method_labels
  recall_df$eps <- eps_vals

  return(list(precision = precision_df, recall = recall_df))
}


# for each k in k_vals:
# returns the proportion of points in both the embedded
# and original k nbhd of a point (averaged over all points)
# this is the same as the recall of the k nbhd and the precision of the k nbhd
get_knn_overlap <- function(orig_distances, emb_distances_list, method_labels, k_vals) {
  # make sure is dist object
  orig_distances <- as.dist(orig_distances)

  m <- length(emb_distances_list)
  for (i in 1:m) {
    emb_distances_list[[i]] <- as.dist(emb_distances_list[[i]])
  }


  # calculate nearest neighbors
  k <- max(k_vals)
  orig_nn <- kNN(orig_distances, k)
  emb_nn <- list()
  for (i in 1:m) {
    emb_nn[[i]] <- kNN(emb_distances_list[[i]], k)
  }

  knn_overlap <- matrix(NA, nrow = length(k_vals), ncol = m)

  for (i in 1:length(k_vals)) { # for each k
    k <- k_vals[i]
    orig_nbhd <- orig_nn$id
    for (j in 1:m) { # for each embedding method
      emb_nbhd <- emb_nn[[j]]$id
      overlap_total <- 0

      # for each point compare embedding neighborhood to original neighborhood
      for (pt in 1:nrow(orig_nbhd)) {
        overlap_count <- length(intersect(orig_nbhd[pt, 1:k], emb_nbhd[pt, 1:k])) # how many points are in both nbhds
        overlap_total <- overlap_total + overlap_count
      }

      knn_overlap[i, j] <- mean(overlap_count) / k
    }
  }

  knn_df <- as_tibble(knn_overlap)
  colnames(knn_df) <- method_labels
  knn_df$k <- k_vals
  return(knn_df)
}


# -------------------------------------------------------------------------
# Plotting and summary functions ------------------------------------------
# -------------------------------------------------------------------------

# input a list of distances from different embedding methods
# returns plots of precision vs. eps and recall vs. eps with a line for each
# embedding method
precision_recall_plots <- function(orig_distances, emb_distances_list,
                                   method_labels, eps_vals = NULL) {
  scores <- get_eps_precision_recall(orig_distances, emb_distances_list, method_labels, eps_vals)

  precision_df <- scores$precision
  precision_df <- precision_df |> pivot_longer(
    cols = method_labels,
    names_to = "Method",
    values_to = "Precision"
  )

  p1 <- ggplot(precision_df, aes(x = eps, y = Precision, color = Method)) +
    geom_line() +
    geom_point() +
    xlab(expression(epsilon))

  recall_df <- scores$recall
  recall_df <- recall_df |> pivot_longer(
    cols = method_labels,
    names_to = "Method",
    values_to = "Recall"
  )

  p2 <- ggplot(recall_df, aes(x = eps, y = Recall, color = Method)) +
    geom_line() +
    geom_point() +
    xlab(expression(epsilon))


  return(list(precision = p1, recall = p2))
}

# returns plot of knn overlap at each k value for each method
nn_overlap_plot <- function(orig_distances, emb_distances_list, method_labels, k_vals) {
  df <- get_knn_overlap(orig_distances, emb_distances_list, method_labels, k_vals)

  df <- df |> pivot_longer(
    cols = method_labels,
    names_to = "Method",
    values_to = "knn_overlap"
  )

  p <- ggplot(df, aes(x = k, y = knn_overlap, color = Method)) +
    geom_line() +
    geom_point() +
    ylab("kNN Overlap") +
    xlab("k")

  return(p)


  # TO DO: add dashed line if random is one of methods
  # ggplot(k_nbhd_df, aes(x = k, y = knn_overlap, color = method, linetype = method)) +
  #   geom_line() +
  #   scale_linetype_manual(values =  c("ours" = "solid", "tsne" = "solid", "random" = "dashed")) +
  #   scale_color_manual(values = c('ours'= 'red', 'tsne' = 'blue', "random" = "black"))
}


# single function to call all global distance evaluation functions
get_metrics <- function(orig_distances, emb_distances) {
  stress_scores <- get_stress(orig_distances, emb_distances)
  spearman <- get_spearman(orig_distances, emb_distances)

  results <- list(
    stress = stress_scores$stress,
    normalized_stress = stress_scores$normalized_stress,
    spearman = spearman
  )

  return(results)
}

# take in a results df with columns metric-method and rows corresponding to the metrics on different runs
summarize_results <- function(results_df) {
  summary <- results_df |> summarize(across(everything(), list(mean = mean, sd = sd), .names = "{.col}.{.fn}"))

  # reformat data
  temp <- summary |> pivot_longer(
    cols = everything(),
    names_to = c("metric-method", "statistic"),
    names_sep = "\\."
  )
  temp <- temp |> pivot_wider(names_from = "statistic", values_from = "value")
  temp <- temp |> separate(col = `metric-method`, into = c("metric", "method"), sep = "-")

  return(temp)
}
