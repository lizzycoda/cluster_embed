source("./scripts/generate_data.R")
source("./scripts/cluster_embed.R")
source("./scripts/eval.R")
source("./scripts/utils.R")
source("./scripts/loe_embedding.R")

#--------------------------------------------------------------------------
# Data loading and preprocessing ------------------------------------------
#--------------------------------------------------------------------------

# load first 5k training sample from MNIST dataset: https://www.kaggle.com/datasets/hojjatk/mnist-dataset
mnist <- read_csv("data/mnist_flat.csv")
labels <- pull(mnist, 785)
X <- as.matrix(mnist[, 1:784]) # flattened version of the data

# calculate geodesic distances
euclidean_distances <- dist(X)
calculate_geodesic_distances(euclidean_distances,
  k = 5, save = TRUE,
  path = "data/MNIST_geodesic_5.csv"
)
# just load after run above function once
distances <- load_geodesic_distances("data/MNIST_geodesic_5.csv", 5000)

# precompute and save clusters
for (i in 1:5) {
  clusters <- get_leiden_clustering(distances, k = 20, resolution = 1 / 2)
  print(rand_index(labels, clusters))
  save_clusters(clusters, paste0("results/mnist/leiden_clusters", i, ".csv"))
}

#--------------------------------------------------------------------------
# Embed the data  ---------------------------------------------------------
#--------------------------------------------------------------------------

# C+E
for (i in 1:5) {
  clusters <- load_clusters(paste0("results/mnist/leiden_clusters", i, ".csv"))
  out <- cluster_align(distances, clusters,
    embedding_method = "PCA", alignment_method = "stress",
    use_grad = T, alpha = 1
  )
  save_embeddings(out$embeddings, labels, cluster_labels = clusters, path = paste0("results/mnist/ce_leiden_run", i, ".csv"))
}

# C+E (alpha = 1.75)
for (i in 1:5) {
  clusters <- load_clusters(paste0("results/mnist/leiden_clusters", i, ".csv"))
  out <- cluster_align(distances, clusters,
    embedding_method = "PCA", alignment_method = "stress",
    use_grad = T, alpha = 1.75
  )
  save_embeddings(out$embeddings, labels, cluster_labels = clusters, path = paste0("results/mnist/ce_leiden_175_run", i, ".csv"))
}

# C+E with LOE

# first only compute the loe embeddings
for (i in 1:5) {
  clusters <- load_clusters(paste0("results/mnist/leiden_clusters", i, ".csv"))
  embeddings <- clusterwise_loe(clusters, distances)
  save_embeddings(embeddings, labels, cluster_labels = clusters, paste0("results/mnist/loe_run", i, ".csv"))
}

# sample of how to determine alpha
get_alpha(distances, dist(embeddings), clusters, .995)

# align cluster embeddings to get global embedding
for (i in 1:5) {
  clusters <- load_clusters(paste0("results/mnist/leiden_clusters", i, ".csv"))
  embeddings <- read_csv(paste0("results/mnist/loe_run", i, ".csv"))
  embeddings <- as.matrix(embeddings[, 1:2])
  out <- cluster_align(distances, clusters,
    embedding_method = "NA", alignment_method = "stress",
    use_grad = T, alpha = 3, embeddings = embeddings
  )
  save_embeddings(out$embeddings, labels, cluster_labels = clusters, path = paste0("results/mnist/ce_leiden_loe_3_run", i, ".csv"))
}


# tuning of tsne parameter
scores <- c()
vals <- seq(50, 1000, by = 50)
for (p in vals) {
  out <- get_tsne(X, perplexity = p)
  spearman <- get_spearman(distances, dist(out))
  scores <- c(scores, spearman)
}
perplexity <- vals[which.max(scores)] # 550

for (i in 1:5) {
  out <- get_tsne(X, perplexity = 550)
  save_embeddings(out, labels, paste0("results/mnist/tsne_550_run", i, ".csv"))
}

# good local behavior
for (i in 1:5) {
  out <- get_tsne(X, perplexity = 20)
  save_embeddings(out, labels, paste0("results/mnist/tsne_20_run", i, ".csv"))
}


pca <- get_pca(euclidean_distances)
save_embeddings(pca, labels, path = "results/mnist/pca.csv")

isomap <- get_isomap(distances, "geodesic")
save_embeddings(isomap, labels, path = "results/mnist/isomap.csv")



#--------------------------------------------------------------------------
# Metrics -----------------------------------------------------------------
#--------------------------------------------------------------------------

embedding_paths <- c(
  paste0("results/mnist/ce_leiden_run", 1:5, ".csv"),
  paste0("results/mnist/ce_leiden_175_run", 1:5, ".csv"),
  paste0("results/mnist/tsne_run", 1:5, ".csv"),
  paste0("results/mnist/tsne_550_run", 1:5, ".csv"),
  paste0("results/mnist/ce_leiden_loe_3_run", 1:5, ".csv"),
  "results/mnist/isomap.csv",
  "results/mnist/pca.csv"
)
embedding_names <- c(
  rep("ce_leiden", 5), rep("ce_leiden2", 5), rep("tsne", 5), rep("tsne_2", 5),
  rep("loe", 5), "isomap", "pca"
)

results <- get_distance_metrics(embedding_paths, embedding_names, distances, labels)
results_euc <- get_distance_metrics(embedding_paths, embedding_names, euclidean_distances, labels)

#--------------------------------------------------------------------------
# Embedding Plots ---------------------------------------------------------
#--------------------------------------------------------------------------

ce1 <- read_csv("results/mnist/ce_leiden_run1.csv") |> mutate(across(3:4, as.factor))
ce2 <- read_csv("results/mnist/ce_leiden_175_run1.csv") |> mutate(across(3:4, as.factor))
ce3 <- read_csv("results/mnist/ce_leiden_loe_3_run1.csv") |> mutate(across(3:4, as.factor))
tsne <- read_csv("results/mnist/tsne_550_run1.csv") |> mutate(across(3, as.factor))
tsne2 <- read_csv("results/mnist/tsne_run1.csv") |> mutate(across(3, as.factor))
pca <- read_csv("results/mnist/pca.csv") |> mutate(across(3, as.factor))
isomap <- read_csv("results/mnist/isomap.csv") |> mutate(across(3, as.factor))

loe_rotated <- procrustes(X = ce2[, 1:2], Y = ce3[, 1:2], scale = F)$Yrot
ce3[, 1:2] <- loe_rotated

tsne_rotated <- procrustes(X = ce2[, 1:2], Y = tsne[, 1:2], scale = F)$Yrot
tsne[, 1:2] <- tsne_rotated

tsne_rotated <- procrustes(X = ce2[, 1:2], Y = tsne2[, 1:2], scale = F)$Yrot
tsne2[, 1:2] <- tsne_rotated

isomap_rotated <- procrustes(X = ce2[, 1:2], Y = isomap[, 1:2], scale = F)$Yrot
isomap[, 1:2] <- isomap_rotated

pca_rotated <- procrustes(X = ce2[, 1:2], Y = pca[, 1:2], scale = F)$Yrot
pca[, 1:2] <- pca_rotated

p1 <- plot(ce1, 60, 60, stepsize = 30)
p2 <- plot(ce2, 100, 100, stepsize = 50)
p3 <- plot(ce3, 150, 150, , stepsize = 75)
p4 <- plot(tsne, 15, 15, , stepsize = 7.5)
p5 <- plot(tsne2, 60, 60, , stepsize = 30)
p6 <- plot(pca, 8, 8, stepsize = 4)
p7 <- plot(isomap, 40, 40, stepsize = 20)

ggsave("figures/mnist_ce.png", plot = p1, width = 6, height = 6, units = "in", dpi = 100)
ggsave("figures/mnist_ce_scaled.png", plot = p2, width = 6, height = 6, units = "in", dpi = 100)
ggsave("figures/mnist_ce_loe.png", plot = p3, width = 6, height = 6, units = "in", dpi = 100)
ggsave("figures/mnist_tsne_550.png", plot = p4, width = 6, height = 6, units = "in", dpi = 100)
ggsave("figures/mnist_tsne_20.png", plot = p5, width = 6, height = 6, units = "in", dpi = 100)
ggsave("figures/mnist_pca.png", plot = p6, width = 6, height = 6, units = "in", dpi = 100)
ggsave("figures/mnist_isomap.png", plot = p7, width = 6, height = 6, units = "in", dpi = 100)
?

#--------------------------------------------------------------------------
# Recall Plots ------------------------------------------------------------
#--------------------------------------------------------------------------

k_vals <- seq(5, 500, 5)
embedding_names <- c(
  paste0("ce1_run", 1:5), paste0("ce175_run", 1:5),
  paste0("tsne20_run", 1:5), paste0("tsne550_run", 1:5),
  paste0("ce_loe_run", 1:5), "isomap", "pca"
)
knn_results <- get_recall_metrics(embedding_paths, embedding_names, distances, k_vals)
write_csv(knn_results, "results/mnist/knn_results.csv")

ce_alpha1_mean <- rowMeans(knn_results[, 1:5])
ce_alpha1_sd <- apply(knn_results[, 1:5], 1, sd)
ce_alpha175_mean <- rowMeans(knn_results[, 6:10])
ce_alpha175_sd <- apply(knn_results[, 6:10], 1, sd)
tsne20_mean <- rowMeans(knn_results[, 11:15])
tsne20_sd <- apply(knn_results[, 11:15], 1, sd)
tsne550_mean <- rowMeans(knn_results[, 16:20])
tsne550_sd <- apply(knn_results[, 16:20], 1, sd)
ce_loe_mean <- rowMeans(knn_results2[, 21:25])
ce_loe_sd <- apply(knn_results2[, 21:25], 1, sd)
isomap <- knn_results[26]
pca <- knn_results[27]

knn_summary <- as_tibble(cbind(
  ce_alpha1_mean, ce_alpha1_sd,
  ce_alpha175_mean, ce_alpha175_sd,
  tsne20_mean, tsne20_sd,
  tsne550_mean, tsne550_sd,
  ce_loe_mean, ce_loe_sd,
  isomap, pca
))

knn_summary$k <- k_vals

ce_cols <- c("ce_alpha1", "ce_alpha175", "ce_loe")
ce_labels <- c("C+E (Isomap, α = 1)", "C+E (Isomap, α = 1.75)", "C+E (LOE, α = 3.00)")
tsne_cols <- c("tsne20", "tsne550")
tsne_labels <- c("t-SNE (u = 20)", "t-SNE (u = 550)")


p <- knn_plot(knn_summary, ce_cols, tsne_cols, ce_labels, tsne_labels)
ggsave("figures/mnist_knn.png", plot = p, width = 10, height = 6, units = "in", dpi = 100)
