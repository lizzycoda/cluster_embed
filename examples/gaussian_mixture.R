source("./scripts/generate_data.R")
source("./scripts/cluster_embed.R")
source("./scripts/eval.R")
source("./scripts/utils.R")

#--------------------------------------------------------------------------
# GMM with equal spacing --------------------------------------------------
#--------------------------------------------------------------------------

# load the data
data <- load_data("gmm_10d", n = 1000)
X <- data$X
labels <- data$labels
distances <- dist(X)

#--------------------------------------------------------------------------
# Embed the data  ---------------------------------------------------------
#--------------------------------------------------------------------------

# our method with no scaling
for (i in 1:5) {
  clusters <- kmeans(X, 10, iter.max = 300)$cluster
  out <- cluster_align(distances, clusters,
    embedding_method = "PCA",
    alignment_method = "stress",
    use_grad = T
  )
  save_embeddings(out$embeddings, labels, cluster_labels = clusters, path = paste0("results/gmm/ce_alpha1_run", i, ".csv"))
}


# our method with scaling (alpha = 2)
for (i in 1:5) {
  clusters <- kmeans(X, 10, iter.max = 300)$cluster

  out <- cluster_align(distances, clusters,
    embedding_method = "PCA",
    alignment_method = "stress",
    use_grad = T, alpha = 2
  )
  save_embeddings(out$embeddings, labels, cluster_labels = clusters, path = paste0("results/gmm/ce_alpha2_run", i, ".csv"))
}

for (i in 1:5) {
  clusters <- kmeans(X, 10, iter.max = 300)$cluster

  out <- cluster_align(distances, clusters,
    embedding_method = "PCA",
    alignment_method = "stress",
    use_grad = T, alpha = 1.75
  )
  save_embeddings(out$embeddings, labels, cluster_labels = clusters, path = paste0("results/gmm/ce_alpha175_run", i, ".csv"))
}

# tsne with perplexity tuned to maximize global spearman score
for (i in 1:5) {
  out <- get_tsne(X, perplexity = 330)
  save_embeddings(out, labels, paste0("results/gmm/tsne_330_run", i, ".csv"))
}

# tsne with perplexity adjusted for `nicer` plots
for (i in 1:5) {
  out <- get_tsne(X, perplexity = 200)
  save_embeddings(out, labels, paste0("results/gmm/tsne_200_run", i, ".csv"))
}

# pca
pca <- get_pca(distances)
save_embeddings(pca, labels, "results/gmm/pca.csv")


# isomap with k tuned to maximize global spearman score
out <- get_isomap(distances, "euclidean", k = 400)
save_embeddings(out, labels, paste0("results/gmm/isomap.csv"))

#--------------------------------------------------------------------------
# Metrics -----------------------------------------------------------------
#--------------------------------------------------------------------------
embedding_paths <- c(
  paste0("results/gmm/ce_alpha1_run", 1:5, ".csv"), paste0("results/gmm/ce_alpha175_run", 1:5, ".csv"),
  paste0("results/gmm/ce_alpha2_run", 1:5, ".csv"), paste0("results/gmm/tsne_200_run", 1:5, ".csv"),
  paste0("results/gmm/tsne_330_run", 1:5, ".csv"), "results/gmm/isomap.csv", "results/gmm/pca.csv"
)

embedding_names <- c(rep("ce", 5), rep("ce175", 5), rep("ce2", 5), rep("tsne200", 5), rep("tsne330", 5), "isomap", "pca")

results <- get_distance_metrics(embedding_paths, embedding_names, distances, labels)


#--------------------------------------------------------------------------
# Embedding Plots ---------------------------------------------------------
#--------------------------------------------------------------------------

ce1 <- read_csv("results/gmm/ce_alpha1_run1.csv") |> mutate(across(3:4, as.factor))
ce2 <- read_csv("results/gmm/ce_alpha2_run1.csv") |> mutate(across(3:4, as.factor))
tsne <- read_csv("results/gmm/tsne_200_run1.csv") |> mutate(across(3, as.factor))
pca <- read_csv("results/gmm/pca.csv") |> mutate(across(3, as.factor))
isomap <- read_csv("results/gmm/isomap.csv") |> mutate(across(3, as.factor))
ce3 <- read_csv("results/gmm/ce_alpha175_run1.csv") |> mutate(across(3:4, as.factor))
tsne2 <- read_csv("results/gmm/tsne_330_run1.csv") |> mutate(across(3, as.factor))

p1 <- plot(ce1, 9, 9, stepsize = 3)
p2 <- plot(ce2, 15, 15, stepsize = 5)
p3 <- plot(tsne, 9, 9, , stepsize = 3)
p4 <- plot(pca, 6, 6, stepsize = 2)
p5 <- plot(isomap, 6, 6, stepsize = 2)
p6 <- plot(ce3, 15, 15, stepsize = 5)
p7 <- plot(tsne2, 6, 6, stepsize = 2)

ggsave("figures/gmm_ce.png", plot = p1, width = 6, height = 6, units = "in", dpi = 100)
ggsave("figures/gmm_ce_scaled2.png", plot = p2, width = 6, height = 6, units = "in", dpi = 100)
ggsave("figures/gmm_tsne.png", plot = p3, width = 6, height = 6, units = "in", dpi = 100)
ggsave("figures/gmm_pca.png", plot = p4, width = 6, height = 6, units = "in", dpi = 100)
ggsave("figures/gmm_isomap.png", plot = p5, width = 6, height = 6, units = "in", dpi = 100)
ggsave("figures/gmm_ce175.png", plot = p6, width = 6, height = 6, units = "in", dpi = 100)
ggsave("figures/gmm_tsne330.png", plot = p7, width = 6, height = 6, units = "in", dpi = 100)


#--------------------------------------------------------------------------
# Recall Plots ------------------------------------------------------------
#--------------------------------------------------------------------------

embedding_names <- c(
  paste0("ce1_run", 1:5), paste0("ce175_run", 1:5), paste0("ce2_run", 1:5),
  paste0("tsne200_run", 1:5), paste0("tsne330_run", 1:5), "isomap", "pca"
)
k_vals <- seq(2, 100, 1)
knn_results <- get_recall_metrics(embedding_paths, embedding_names, distances, k_vals)
# write_csv(knn_results, 'results/gmm/knn_results.csv')

ce_alpha1_mean <- rowMeans(knn_results[, 1:5])
ce_alpha1_sd <- apply(knn_results[, 1:5], 1, sd)
ce_alpha175_mean <- rowMeans(knn_results[, 6:10])
ce_alpha175_sd <- apply(knn_results[, 6:10], 1, sd)
ce_alpha2_mean <- rowMeans(knn_results[, 11:15])
ce_alpha2_sd <- apply(knn_results[, 11:15], 1, sd)
tsne200_mean <- rowMeans(knn_results[, 16:20])
tsne200_sd <- apply(knn_results[, 16:20], 1, sd)
tsne330_mean <- rowMeans(knn_results[, 21:25])
tsne330_sd <- apply(knn_results[, 21:25], 1, sd)
isomap <- knn_results[26]
pca <- knn_results[27]

knn_summary <- as_tibble(cbind(
  ce_alpha1_mean, ce_alpha1_sd,
  ce_alpha175_mean, ce_alpha175_sd,
  ce_alpha2_mean, ce_alpha2_sd,
  tsne200_mean, tsne200_sd,
  tsne330_mean, tsne330_sd, isomap, pca
))
knn_summary$k <- seq(2, 100, 1)

ce_cols <- c("ce_alpha1", "ce_alpha175", "ce_alpha2")
ce_labels <- c("C+E (α = 1)", "C+E (α = 1.75)", "C+E (α = 2)")
tsne_cols <- c("tsne200", "tsne330")
tsne_labels <- c("t-SNE (u = 200)", "t-SNE (u = 330)")

p <- knn_plot(knn_summary, ce_cols, tsne_cols, ce_labels, tsne_labels)

ggsave("figures/gmm_knn.png", plot = p, width = 10, height = 6, units = "in", dpi = 100)
