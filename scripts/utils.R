
# -------------------------------------------------------------------------
# Format embeddings -------------------------------------------------------
# -------------------------------------------------------------------------
# Functions to reformat embedding object to aligned embeddings ready for visualization


format_all_embeddings <- function(embedding_obj) {
  clusters <- embedding_obj$clusters
  embeddings <- embedding_obj$embeddings
  for (i in 1:max(clusters)) {
    embeddings[which(clusters == i), ] <- get_cluster_embeddings(i, embedding_obj)
  }
  return(embeddings)
}


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
