library(dbscan)

load_data <- function(dataset = c("shapes_2d", "shapes_3d", "halfsphere", "gmm_10d"), n = 1000) {
  dataset <- match.arg(dataset)
  
  set.seed(321)

  if (grepl("shapes", dataset)) {
    # load geometric shapes dataset from dbscan package
    data("DS3", package = "dbscan")
    X <- DS3

    # cluster the data and remove noise
    cl <- hdbscan(X, minPts = 25) #params taken from https://cran.r-project.org/web/packages/dbscan/vignettes/hdbscan.html
    X <- as.matrix(X[cl$cluster > 0, ])
    X <- scale(X, scale = FALSE) # centering 
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
      X <- scale(X, scale = FALSE) # center the data
      labels <- c(rep(1, nrow(X1)), rep(2, nrow(X2)), rep(3, nrow(X3)))
    }
  }
  
  else if(dataset == 'halfsphere'){
    r = .15 
    
    # centers of 10 clusters 
   centers = matrix(c(0.11315152, -0.2910270, 0.95,
                      -0.47240967, 0.2330861, 0.85,
                      0.63507596, 0.1848744, 0.75,
                      -0.39454163, -0.6494897, 0.65,
                      -0.16243302, 0.8192164, 0.55,
                      0.71978430, -0.5285930, 0.45,
                      -0.93127144, -0.1011608, 0.35,
                      0.63914835, 0.7273166, 0.25,
                      0.02042954, -0.9884749, 0.15,
                      -0.68971880, 0.7223489, 0.05),
                    nrow = 10, ncol = 3, byrow = TRUE)
   
    
    labels = c()
    X = c()
    pts = uniform_sphere(25*n, min_z = -r) # oversample
    for(i in 1:nrow(centers)){
      b = centers[i,]
      dists = dist_on_sphere(pts,b)
      cluster_i = pts[which(dists < r),]
      labels = c(labels, rep(i, nrow(cluster_i)))
      X = rbind(X, cluster_i)
      
    }
    
    
  }
  
  # else if(dataset == 'gmm_10d'){
  #   
  #   
  # }


  max_n <- nrow(X)
  idxs <- sample.int(max_n, size = min(n, max_n))
  X <- X[idxs,]
  labels <- labels[idxs]

  return(list(X = X, labels = labels))
}


# Helper functions for halfsphere dataset ---------------------------------


# sample n uniform points on a sphere
uniform_sphere <- function(n, min_z){
  X = rnorm(n)
  Y = rnorm(n)
  Z = rnorm(n)
  R = sqrt(X^2 + Y^2 + Z^2 )
  sphere = cbind(X/R, Y/R, Z/R)
  halfsphere = sphere[which(sphere[,3] > min_z),]
  
  return(halfsphere)
}

# distance on sphere between every point in A and the vector b 
dist_on_sphere <- function(A,b, r = 1){
  return(r*acos(A %*% b / r^2))
}
