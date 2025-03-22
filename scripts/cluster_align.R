library(stats)
library(smacofx)

source("./scripts/utils.R")
source("./scripts/stress_functions.R")

# -------------------------------------------------------------------------
# Main method -------------------------------------------------------------
# -------------------------------------------------------------------------


cluster_align <- function(distances,
                          clusters,
                          embedding_method = c("MDS", "LMDS", "KL"),
                          alignment_method = c("stress", "s-stress", "sammon"),
                          alpha = 1,
                          learn_translations = T,
                          max_iter = 10,
                          max_init_iter = 10,
                          embeddings = NULL,
                          use_grad = F) {
  embedding_method <- match.arg(embedding_method)
  alignment_method <- match.arg(alignment_method)

  orig_labels <- clusters # store original cluster labels
  clusters <- reindex_clusters(clusters) # index from 1:num_clusters
  n_clusters = max(clusters)

  distances <- as.matrix(distances) # will be easier to work with this format

  # embed each cluster by the embedding method unless the embeddings are prespecified
  if (is.null(embeddings)) {
    embedding_obj <- embed_clusters(distances, clusters, embedding_method)
  } else {
    embedding_obj <- list(
      embeddings = embeddings,
      translations = matrix(0, nrow = n_clusters, ncol = 2),
      rotations = matrix(0, nrow = n_clusters, ncol = 2),
      reflections = matrix(0, nrow = n_clusters, ncol = 2),
      clusters = clusters
    )
  }


  if (alpha != 1) {
    distances <- scale_distances(distances, clusters, alpha)
  }

  # initialize translations on each cluster
  embedding_obj$translations <- initialize_translations(distances, clusters, method = "MDS")


  ###### Alignment stage 1 #####

  # Fix the largest cluster and orient the remaining clusters with respect to only this cluster
  cluster_idxs <- 1:max(clusters)
  cluster_counts <- table(clusters)
  fixed_cluster <- as.numeric(names(cluster_counts)[which.max(cluster_counts)]) # fix the largest cluster
  moveable_clusters <- cluster_idxs[cluster_idxs != fixed_cluster]

  thresh = .01
  for (j in moveable_clusters) {
    for (iter in 1:max_init_iter) {
      out <- align_cluster(j,
        anchors = c(fixed_cluster),
        embedding_obj = embedding_obj,
        distances = distances,
        alignment_method = alignment_method,
        learn_translations = learn_translations,
        use_grad = use_grad
      )
      
      embedding_obj <- out$embedding_obj
      curr_loss <- out$loss
    

    if (iter > 1) { # stop updating a cluster once loss stabilizes
      if (abs(curr_loss - prev_loss) < thresh) {
        break
      }
    }
    prev_loss <- curr_loss
    }
  }


  # Alignment stage 2
  # Fix the largest cluster and orient the remaining clusters
  for (iter in 1:max_iter) {
    moveable_clusters <- sample(moveable_clusters) # shuffle order
    for (i in moveable_clusters) {
      out <- align_cluster(i,
        anchors = cluster_idxs[cluster_idxs != i],
        embedding_obj = embedding_obj,
        distances = distances,
        alignment_method = alignment_method,
        learn_translations = learn_translations,
        use_grad = use_grad
      )
      embedding_obj <- out$embedding_obj
    }
  }

  # Alignment stage 3
  # Allow all clusters to update
  for (iter in 1:max_iter) {
    moveable_clusters <- sample(cluster_idxs) # shuffle order of all clusters
    for (i in moveable_clusters) {
      out <- align_cluster(i,
        anchors = cluster_idxs[cluster_idxs != i],
        embedding_obj = embedding_obj,
        distances = distances,
        alignment_method = alignment_method,
        learn_translations = learn_translations,
        use_grad = use_grad
      )
      embedding_obj <- out$embedding_obj
    }
  }

  # reformat final embedding based on learned parameters
  out <- format_all_embeddings(embedding_obj)

  return(list(embeddings = out, clusters = orig_labels))
}




# -------------------------------------------------------------------------
# Helper functions --------------------------------------------------------
# -------------------------------------------------------------------------



# embed individual clusters
# returns object with embedding of each point,  translation (num_clusters x2) of each cluster,
embed_clusters <- function(distances, clusters, embedding_method) {
  n_clusters <- max(clusters)

  if (embedding_method == "MDS") {
    embed_function <- stats::cmdscale
  } else if (embedding_method == "LMDS") {
    embed_function <- function(distances) {
      Y <- smacofx::lmds(distances, k = k, tau = tau)
      return(Y$conf)
    }
  } else if (embedding_method == "KL") {
    embed_function <- function(distances) {
      return(kl_func(distances, k = k))
    }
  }


  embeddings <- matrix(NA, nrow = length(clusters), ncol = 2)
  for (j in 1:n_clusters) {
    idxs <- which(clusters == j)
    embeddings[idxs, ] <- embed_function(distances[idxs, idxs])
  }


  embedding_obj <- list(
    embeddings = embeddings,
    translations = matrix(0, nrow = n_clusters, ncol = 2),
    rotations = matrix(0, nrow = n_clusters, ncol = 2),
    reflections = matrix(0, nrow = n_clusters, ncol = 2),
    clusters = clusters
  )

  return(embedding_obj)
}


# initialize translation on each cluster
initialize_translations <- function(distances, clusters, method = c("MDS", "random")) {
  method <- match.arg(method)

  n_clusters <- max(clusters)

  if (method == "random") {
    scale <- median(distances) # coarse approximation of scale of embedding
    translations <- matrix(runif(2 * n_clusters, -scale, scale), ncol = 2)
  } else if (method == "MDS") {
    mean_cluster_distances <- matrix(0, nrow = n_clusters, ncol = n_clusters)

    for (i in 1:(n_clusters - 1)) {
      cluster_i_idxs <- which(clusters == i)
      for (j in (i + 1):n_clusters) {
        cluster_j_idxs <- which(clusters == j)
        batch_distances <- distances[cluster_i_idxs, cluster_j_idxs]
        mean_cluster_distances[i, j] <- mean(batch_distances)
        mean_cluster_distances[j, i] <- mean_cluster_distances[i, j]
      }
    }

    if (n_clusters == 2) { # cmd scale gives an error in this setting
      translations <- matrix(c(
        mean_cluster_distances[1, 2] / 2,
        -mean_cluster_distances[1, 2] / 2, 0, 0
      ), nrow = 2)
    } else {
      translations <- stats::cmdscale(d = mean_cluster_distances, k = 2)
    }
  }

  return(translations)
}



# find rigid transformation on cluster to minimize the stress
align_cluster <- function(cluster_idx, 
                          anchors,
                          embedding_obj,
                          distances,
                          alignment_method = "stress",
                          learn_translations = T,
                          use_grad = F) {
  rotation <- embedding_obj$rotations[cluster_idx]
  translation <- embedding_obj$translations[cluster_idx, ]


  # find best rotation matrix for cluster
  rot_out <- optimize(rotation_stress,
    interval = c(0, 2 * pi),
    cluster_idx = cluster_idx,
    reflection = 0,
    translation = translation,
    anchors = anchors,
    embedding_obj = embedding_obj,
    distances = distances,
    alignment_method = alignment_method
  )
  rot_loss <- rot_out$objective


  # find best reflection matrix for cluster
  ref_out <- optimize(rotation_stress,
    interval = c(0, 2 * pi),
    cluster_idx = cluster_idx,
    reflection = 1,
    translation = translation,
    anchors = anchors,
    embedding_obj = embedding_obj,
    distances = distances,
    alignment_method = alignment_method
  )
  ref_loss <- ref_out$objective

  # update reflection and rotation parameters
  if (ref_loss < rot_loss) {
    reflection <- 1
    rotation <- ref_out$minimum
  } else {
    reflection <- 0
    rotation <- rot_out$minimum
  }
  embedding_obj$rotations[cluster_idx] <- rotation
  embedding_obj$reflections[cluster_idx] <- reflection
  loss = min(ref_loss, rot_loss)


  # optionally update the translation of the cluster
  if (learn_translations) {
    if (use_grad == T) {
      if (alignment_method != "stress") {
        stop("Gradient method has not been implemented")
      }

      # in practice this seemed to work better than gradient descent which required extensive learning rate tuning
      translation_out <- optim(translation, translation_stress,
        gr = cluster_translation_gradient,
        cluster_idx = cluster_idx,
        rotation = rotation,
        reflection = reflection,
        anchors = anchors,
        embedding_obj = embedding_obj,
        distances = distances,
        alignment_method = alignment_method,
        method = "BFGS",
        control = list(maxit = 10)
      )
    } else {
      translation_out <- optim(translation, translation_stress,
        cluster_idx = cluster_idx,
        rotation = rotation,
        reflection = reflection,
        anchors = anchors,
        embedding_obj = embedding_obj,
        distances = distances,
        alignment_method = alignment_method,
        method = "Nelder-Mead",
        control = list(maxit = 10)
      )
    }
    embedding_obj$translations[cluster_idx, ] <- translation_out$par
    loss = translation_out$value
  }


  return(list(embedding_obj = embedding_obj, loss = loss ))
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




# scale distances between any points not in the same cluster by alpha
scale_distances <- function(distances, clusters, alpha) {
  for (i in 1:max(clusters)) {
    idxs1 <- which(clusters == i)
    idxs2 <- which(clusters != i)
    distances[idxs1, idxs2] <- alpha * distances[idxs1, idxs2]
  }
  return(distances)
}
