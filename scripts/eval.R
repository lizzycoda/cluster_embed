library(ggplot2)
library(tidyr)
library(dbscan)

# -------------------------------------------------------------------------
# Distance-based evaluation functions -------------------------------------
# -------------------------------------------------------------------------

# calculates the spearman correlation of the original distances and embedded distances
get_spearman <- function(orig_distances, emb_distances) {
  orig_distances <- as.dist(orig_distances)
  emb_distances <- as.dist(emb_distances)

  return(cor(orig_distances, emb_distances, method = "spearman"))
}

# calculates the overall stress and normalized stress
get_stress <- function(orig_distances, emb_distances) {
  orig_distances <- as.dist(orig_distances)
  emb_distances <- as.dist(emb_distances)

  stress <- sum((orig_distances - emb_distances)^2)
  normalized_stress <- stress / sum(orig_distances^2)
  return(list(
    stress = stress,
    normalized_stress = sqrt(normalized_stress)
  ))
}


get_class_preservation <- function(orig_distances, emb_distances, classes) {
  orig_distances <- as.matrix(orig_distances)
  emb_distances <- as.matrix(emb_distances)

  unique_classes <- unique(classes)
  m <- length(unique(classes))

  # calculate average pairwise distance between classes
  avg_class_distances_orig <- matrix(NA, nrow = m, ncol = m)
  avg_class_distances_emb <- matrix(NA, nrow = m, ncol = m)
  for (i in 1:(m - 1)) {
    idxs1 <- which(classes == unique_classes[i])
    for (j in (i + 1):m) {
      idxs2 <- which(classes == unique_classes[j])
      avg_class_distances_orig[i, j] <- mean(orig_distances[idxs1, idxs2])
      avg_class_distances_orig[j, i] <- avg_class_distances_orig[i, j]

      avg_class_distances_emb[i, j] <- mean(emb_distances[idxs1, idxs2])
      avg_class_distances_emb[j, i] <- avg_class_distances_emb[i, j]
    }
  }

  # for each class calculate spearman correlation between distances to other classes
  correlation_scores <- c()
  for (i in 1:m) {
    orig <- na.omit(avg_class_distances_orig[i, ])
    emb <- na.omit(avg_class_distances_emb[i, ])

    correlation_scores <- c(
      correlation_scores,
      cor(orig, emb, method = "spearman")
    )
  }

  return(mean(correlation_scores))
}


# -------------------------------------------------------------------------
#
# -------------------------------------------------------------------------


# single function to call all global distance evaluation functions
get_global_metrics <- function(orig_distances, emb_distances, classes) {
  stress_scores <- get_stress(orig_distances, emb_distances)
  spearman <- get_spearman(orig_distances, emb_distances)
  class_preservation <- get_class_preservation(orig_distances, emb_distances, classes)

  results <- c(spearman, stress_scores$normalized_stress, class_preservation)

  return(results)
}

# returns average over each class
get_local_metrics <- function(orig_distances, emb_distances, classes) {
  class_labels <- unique(classes)
  n <- length(class_labels)
  spearman_scores <- c()
  stress_scores <- c()

  orig_distances <- as.matrix(orig_distances)
  emb_distances <- as.matrix(emb_distances)


  for (i in class_labels) {
    idxs <- which(classes == i)
    stress <- get_stress(orig_distances[idxs, idxs], emb_distances[idxs, idxs])$normalized_stress
    spearman <- get_spearman(orig_distances[idxs, idxs], emb_distances[idxs, idxs])


    stress_scores <- c(stress_scores, stress)
    spearman_scores <- c(spearman_scores, spearman)
  }
  results <- c(mean(spearman_scores), mean(stress_scores))
  return(results)
}


get_distance_metrics <- function(embedding_paths, embedding_names, distances, labels) {
  n <- length(embedding_paths)
  local_results_df <- matrix(NA, nrow = n, ncol = 2)
  global_results_df <- matrix(NA, nrow = n, ncol = 3)

  for (i in 1:n) {
    run_i <- read_csv(embedding_paths[i])
    distances_i <- dist(run_i[, 1:2])
    local_results_df[i, ] <- get_local_metrics(distances, distances_i, labels)
    global_results_df[i, ] <- get_global_metrics(distances, distances_i, labels)
  }

  local_results_df <- as_tibble(local_results_df)
  global_results_df <- as_tibble(global_results_df)
  colnames(local_results_df) <- c("spearman", "n_stress")
  colnames(global_results_df) <- c("spearman", "n_stress", "class_pres")
  local_results_df["method"] <- embedding_names
  global_results_df["method"] <- embedding_names

  local_summary <- summarize(local_results_df)
  global_summary <- summarize(global_results_df)

  return(list(local = local_summary, global = global_summary))
}

summarize <- function(df){
  summary <- df |> group_by(method) |> 
    summarise(
      across(
        .cols = where(is.numeric),
        .fns = list(mean = ~mean(.x, na.rm = TRUE),
                    sd = ~sd(.x, na.rm = TRUE)),
        .names = "{.col}_{.fn}"
      ),
      .groups = "drop"
    )
  return(summary)
}

# -------------------------------------------------------------------------
# Recall functions --------------------------------------------------------
# -------------------------------------------------------------------------

get_recall_metrics <- function(embedding_paths, embedding_names, orig_distances, k_vals) {
  m <- length(embedding_paths)
  recall_df <- matrix(NA, nrow = m, ncol = length(k_vals))

  orig_distances <- as.dist(orig_distances)
  emb_distances_list <- list()

  for (i in 1:m) {
    run_i <- read_csv(embedding_paths[i])
    distances_i <- dist(run_i[, 1:2])
    emb_distances_list[[i]] <- distances_i
  }

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

      knn_overlap[i, j] <- (overlap_total / nrow(orig_nbhd)) / k
    }
  }



  knn_df <- as_tibble(knn_overlap)
  colnames(knn_df) <- embedding_names
  knn_df$k <- k_vals
  return(knn_df)
}


#--------------------------------------------------------------------------
# Plotting functions ------------------------------------------------------
#--------------------------------------------------------------------------


# plot the embeddings
plot <- function(df, x_max, y_max, stepsize = .5, cmap = "tab10") {
  if (cmap == "shapesI") {
    colors <- c("#1F77B4", "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF")
  } else if (cmap == "shapesII") {
    colors <- c("#17BECF", "#E377C2", "#BCBD22")
  } else if (cmap == "tab10") {
    colors <- c(
      "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF"
    )
  } else if (cmap == "rainbow6") {
    colors <- c("royalblue", "skyblue", "lightgreen", "gold", "tomato", "maroon", "navy")
  } else if (cmap == "rainbow60") {
    colors <- c(
      "#E76BF3", "#EF67EC", "#F564E3", "#FA62DB", "#FD61D1", "#FF61C7", "#FF62BC", "#FF64B0", "#FF67A4", "#FF6A98",
      "#FE6E8A", "#FB727C", "#F47B5C", "#EF7F49", "#EA8331", "#E58700", "#DE8C00", "#D89000", "#D09400", "#FFB600",
      "#FCB336FF", "#F9BD39FF", "#F4C73AFF", "#EECF3AFF", "#E6D839FF", "#FFC800", "#FFDB00", "#FFED00", "#FFFF00", "#39B600",
      "#00B816", "#00BA38", "#00BB4E", "#00BD5F", "#00BE6F", "#00BF7D", "#00C08B", "#00C097", "#00C1A3", "#00C0AF",
      "#00C0BA", "#00BFC4", "#00BECE", "#00BCD8", "#00BAE0", "#00B7E9", "#00B4F0", "#00B0F6", "#00ACFC", "#00A7FF",
      "#35A2FF", "#3D9EFEFF", "#619CFF", "#7E96FF", "#9590FF", "#A889FF", "#B983FF", "#C77CFF", "#D376FF", "#DE70F9"
    )
  }



  p <- ggplot(df, aes(x = V1, y = V2, color = labels)) +
    geom_point(size = 2) +
    scale_x_continuous(
      limits = c(-x_max, x_max),
      breaks = seq(-x_max, x_max, by = stepsize)
    ) +
    scale_y_continuous(
      limits = c(-y_max, y_max),
      breaks = seq(-y_max, y_max, by = stepsize)
    ) +
    coord_fixed() +
    theme_bw() +
    theme(
      axis.text = element_text(size = 14, color = "black"),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "none",
      axis.line = element_line(color = "black"),
      panel.grid = element_blank(),
      plot.title = element_text(size = 20, hjust = 0.5)
    ) +
    scale_color_manual(values = colors)
  return(p)
}


# knn plots
knn_ce_colors <- c("#ff9896", "#ffbb78", "#c7c7c7")
knn_tsne_colors <- c("#98df8a", "#aec7e8")
knn_isomap_color <- "#c5b0d5"
knn_pca_color <- "#f7b6d2"

knn_plot <- function(knn_summary, ce_cols, tsne_cols, ce_labels, tsne_labels) {
  values <- character(0)

  p <- ggplot(knn_summary, aes(x = k))

  for (i in 1:length(ce_cols)) {
    label <- ce_labels[i]

    p <- p + geom_ribbon(aes(
      ymin = .data[[paste0(ce_cols[i], "_mean")]] - .data[[paste0(ce_cols[i], "_sd")]],
      ymax = .data[[paste0(ce_cols[i], "_mean")]] + .data[[paste0(ce_cols[i], "_sd")]],
    ), alpha = 0.3, fill = knn_ce_colors[i]) +
      geom_line(aes(y = .data[[paste0(ce_cols[i], "_mean")]], color = !!label), linewidth = 1.25)
    values[label] <- knn_ce_colors[i]
  }

  for (i in 1:length(tsne_cols)) {
    label <- tsne_labels[i]

    p <- p + geom_ribbon(aes(
      ymin = .data[[paste0(tsne_cols[i], "_mean")]] - .data[[paste0(tsne_cols[i], "_sd")]],
      ymax = .data[[paste0(tsne_cols[i], "_mean")]] + .data[[paste0(tsne_cols[i], "_sd")]],
    ), alpha = 0.3, fill = knn_tsne_colors[i]) +
      geom_line(aes(y = .data[[paste0(tsne_cols[i], "_mean")]], color = !!label), linewidth = 1.25)
    values[label] <- knn_tsne_colors[i]
  }


  p <- p + geom_line(aes(y = isomap, color = "Isomap"), linewidth = 1.25) +
    geom_line(aes(y = pca, color = "PCA"), linewidth = 1.25)

  values["Isomap"] <- knn_isomap_color
  values["PCA"] <- knn_pca_color



  p <- p +
    labs(x = "k", y = "kNN Recall") +
    theme_bw() +
    theme(
      axis.text = element_text(size = 18, color = "black"),
      axis.line = element_line(color = "black"),
      panel.grid = element_blank(),
      plot.title = element_text(size = 20, hjust = 0.5),
      legend.text = element_text(size = 20),
      legend.key.size = unit(2.0, "lines"),
      axis.title.x = element_text(size = 18),
      axis.title.y = element_text(size = 18),
      axis.text.x = element_text(size = 18),
      axis.text.y = element_text(size = 18)
    )


  p <- p + scale_color_manual(
    name = "",
    values = values,
    breaks = c(ce_labels, tsne_labels, "Isomap", "PCA")
  )


  return(p)
}
