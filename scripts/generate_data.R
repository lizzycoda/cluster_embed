library(dbscan)

load_data <- function(dataset = c("shapes_2d", "shapes_3d"), n = 1000) {
  dataset <- match.arg(dataset)

  if (grepl("shapes", dataset)) {
    # load geometric shapes dataset from dbscan package
    data("DS3", package = "dbscan")
    X <- DS3

    # cluster the data and remove noise
    cl <- hdbscan(X, minPts = 25) #params taken from https://cran.r-project.org/web/packages/dbscan/vignettes/hdbscan.html
    X <- as.matrix(X[cl$cluster > 0, ])
    labels <- cl$cluster[cl$cluster > 0]

    if (dataset == "shapes_3d") {
      X1 <- X[which(labels == 3), ]
      X2 <- X[which(labels == 4), ]
      X3 <- X[which(labels == 5), ]

      # embed into 3d
      X1 <- cbind(X1, rep(0, nrow(X1)))
      X2 <- cbind(X2[, 1], rep(0, nrow(X2)), X2[, 2])
      X3 <- cbind(rep(0, nrow(X3)), X3)
      X <- rbind(X1, X2, X3)
      labels <- c(rep(1, nrow(X1)), rep(2, nrow(X2)), rep(3, nrow(X3)))
    }
  }


  set.seed(321)
  max_n <- nrow(X)
  idxs <- sample.int(max_n, size = min(n, max_n))
  X <- scale(X[idxs, ], scale = FALSE) # center the data
  labels <- labels[idxs]

  return(list(X = X, labels = labels))
}
