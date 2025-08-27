library(dbscan)
library(MASS)
library(stats)


load_data <- function(dataset = c("shapes_2d", "shapes_3d", "shapes_cylinder", "gmm_10d"), n = 1000) {
  dataset <- match.arg(dataset)

  set.seed(321)

  if (dataset %in% c("shapes_2d", "shapes_3d")) {
    data("DS3", package = "dbscan") # load geometric shapes dataset from dbscan package
    X <- DS3

    cl <- hdbscan(X, minPts = 25) # cluster the data and remove noise
    X <- as.matrix(X[cl$cluster > 0, ])
    X <- scale(X, scale = FALSE) # centering
    labels <- cl$cluster[cl$cluster > 0]

    if (dataset == "shapes_3d") {
      X1 <- X[which(labels == 3), ]
      X2 <- X[which(labels == 4), ]
      X3 <- X[which(labels == 5), ]

      # embed into 3d
      X1 <- cbind(X1, rep(0, nrow(X1)))
      X1[, 1] <- X1[, 1] + 100 # shift
      X2 <- cbind(X2[, 1], rep(0, nrow(X2)), X2[, 2])
      X3 <- cbind(rep(0, nrow(X3)), X3)
      X <- rbind(X1, X2, X3)
      X <- scale(X, scale = FALSE) # center the data
      X <- X / max(abs(X))
      labels <- c(rep(1, nrow(X1)), rep(2, nrow(X2)), rep(3, nrow(X3)))
    }
  }


  # project shapes onto a half cylinder
  else if (dataset == "shapes_cylinder") {
    data <- load_data("shapes_2d", n = n) # recursive call
    X <- data$X
    labels <- data$labels

    # normalize
    X <- X / max(abs(X))

    # points on a cylinder
    Z <- matrix(0, nrow = 1000, ncol = 3)
    theta <- (1 - X[, 1]) * pi / 2
    r <- 2 / pi
    Z[, 1] <- r * cos(theta)
    Z[, 2] <- r * sin(theta)
    Z[, 3] <- X[, 2]
    X <- Z
  }


  # generate a gaussian mixture in 10d with equally spaced components
  else if (dataset == "gmm_10d") {
    d <- 10 # dimension
    k <- 10 # number of clusters
    sep <- 1.5 * sqrt(d) # separation between clusters (within a cluster points are about sqrt(d) from center)
    cluster_counts <- as.vector(rmultinom(n = 1, size = n, prob = rep(1, k)))

    X <- rep(0, d)
    labels <- c()
    for (i in 1:k) {
      mu <- rep(0, d)
      mu[i] <- sep

      X1 <- mvrnorm(cluster_counts[i], mu, 1 * diag(d))
      X <- rbind(X, X1)

      labels <- c(labels, rep(i, cluster_counts[i]))
    }
    X <- X[2:(n + 1), ] # drop first row which was initialized as a placeholder
  }

  max_n <- nrow(X)
  idxs <- sample.int(max_n, size = min(n, max_n))
  X <- X[idxs, ]
  labels <- labels[idxs]

  return(list(X = X, labels = labels))
}
