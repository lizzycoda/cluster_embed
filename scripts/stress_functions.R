library(proxy)

# -------------------------------------------------------------------------
# Stress functions --------------------------------------------------------
# -------------------------------------------------------------------------


# Calculates the stress (or stress variant) with contributions from points in cluster_idx
# and point in the clusters in the set of anchors
#
cluster_stress <- function(cluster_idx, rotation, reflection, translation,
                           anchors, embedding_obj, distances, alignment_method) {
  stress <- 0
  clusters <- embedding_obj$clusters

  cluster1_embeddings <- get_cluster_embeddings(cluster_idx, embedding_obj, rotation, reflection, translation)
  cluster1_idxs <- which(clusters == cluster_idx)

  for (idx in anchors) {
    cluster2_embeddings <- get_cluster_embeddings(idx, embedding_obj)

    cluster2_idxs <- which(clusters == idx)
    orig_distances <- distances[cluster1_idxs, cluster2_idxs]
    emb_distances <- proxy::dist(cluster1_embeddings, cluster2_embeddings)

    stress <- stress + pairwise_stress(orig_distances, emb_distances, alignment_method)
  }
  return(stress)
}

pairwise_stress <- function(orig_distances, emb_distances, alignment_method = c("stress", "s-stress", "sammon_stress")) {
  if (alignment_method == "stress") {
    return(sum((orig_distances - emb_distances)**2))
  } else if (alignment_method == "s-stress") {
    return(sum((orig_distances^2 - emb_distances^2)**2))
  } else if (alignment_method == "sammon-stress") {
    return(sum(1 / orig_distances * (orig_distances - emb_distances)**2))
  }
}


# wrapper function for use with optimization function
rotation_stress <- function(rotation, cluster_idx, reflection, translation, anchors,
                            embedding_obj, distances, alignment_method) {
  return(cluster_stress(
    cluster_idx, rotation, reflection, translation, anchors,
    embedding_obj, distances, alignment_method
  ))
}


# wrapper function for use with optimization function
translation_stress <- function(translation, cluster_idx, rotation, reflection, anchors,
                               embedding_obj, distances, alignment_method) {
  return(cluster_stress(
    cluster_idx, rotation, reflection, translation, anchors,
    embedding_obj, distances, alignment_method
  ))
}


#  ------------------------------------------------------------------------
# Gradient functions ------------------------------------------------------
#  ------------------------------------------------------------------------
# Currently only implemented to support stress minimization and not variants of the stress


# gradient of stress with respect to translation vector for cluster i
cluster_translation_gradient <- function(translation, cluster_idx, rotation, reflection, anchors,
                                         embedding_obj, distances, alignment_method) {
  grad <- c(0, 0)
  clusters <- embedding_obj$clusters

  cluster1_embeddings <- get_cluster_embeddings(cluster_idx, embedding_obj, rotation, reflection, translation)
  cluster1_idxs <- which(clusters == cluster_idx)

  for (idx in anchors) {
    cluster2_embeddings <- get_cluster_embeddings(idx, embedding_obj)

    cluster2_idxs <- which(clusters == idx)
    orig_distances <- distances[cluster1_idxs, cluster2_idxs]
    emb_distances <- proxy::dist(cluster1_embeddings, cluster2_embeddings)

    grad <- grad + pairwise_grad(cluster1_embeddings, cluster2_embeddings, orig_distances, emb_distances)
  }
  return(grad)
}

pairwise_grad <- function(embeddings1, embeddings2, orig_distances, emb_distances) {
  embedding_x_diff <- outer(embeddings1[, 1], embeddings2[, 1], "-")
  embedding_y_diff <- outer(embeddings1[, 2], embeddings2[, 2], "-")

  dx <- sum((emb_distances - orig_distances) * embedding_x_diff / emb_distances)
  dy <- sum((emb_distances - orig_distances) * embedding_y_diff / emb_distances)

  return(c(dx, dy))
}

# -------------------------------------------------------------------------
# Stress helper functions --------------------------------------------------------
# -------------------------------------------------------------------------
# used in calculation of the stress

get_cluster_embeddings <- function(cluster_idx, embedding_obj,
                                   rotation = NULL, reflection = NULL, translation = NULL) {
  # if parameters are provided then use those, otherwise use the default settings
  if (is.null(rotation)) {
    rotation <- embedding_obj$rotations[cluster_idx]
  }

  if (is.null(reflection)) {
    reflection <- embedding_obj$reflections[cluster_idx]
  }

  if (is.null(translation)) {
    translation <- embedding_obj$translations[cluster_idx, ]
  }


  idxs <- which(embedding_obj$clusters == cluster_idx)
  out <- rigid_rotation(embedding_obj$embeddings[idxs, ], rotation, reflection, translation)
  return(out)
}


rigid_rotation <- function(Z, rotation, reflection, translation) {
  # apply rigid rotation to object Z
  R <- orthogonal_matrix(rotation, reflection)
  Z_transformed <- Z %*% t(R) # Z is nx2
  Z_transformed[, 1] <- Z_transformed[, 1] + translation[1]
  Z_transformed[, 2] <- Z_transformed[, 2] + translation[2]

  return(Z_transformed)
}


orthogonal_matrix <- function(theta, reflection = 0) {
  R <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nrow = 2)
  if (reflection == 1) {
    R <- matrix(c(1, 0, 0, -1), nrow = 2) %*% R
  }
  return(R)
}
