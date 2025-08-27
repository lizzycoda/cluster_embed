source("./scripts/generate_data.R")
source("./scripts/cluster_embed.R")
source("./scripts/eval.R")
source("./scripts/utils.R")
source("./scripts/loe_embedding.R")

#--------------------------------------------------------------------------
# Data loading and preprocessing ------------------------------------------
#--------------------------------------------------------------------------

data <- read_csv("data/inhNeurons.csv")
X <- as.matrix(data[1:3000])
labels <- pull(data, "labels")

euclidean_distances <- dist(X)
out <- calculate_geodesic_distances(euclidean_distances,
  k = 10, save = TRUE,
  path = "data/neurons_geodesic_10.csv"
)

distances <- load_geodesic_distances("data/neurons_geodesic_10.csv", 5000)


# precompute clusters and remove clusters with only 1 points
set.seed(50)
clusters <- get_leiden_clustering(distances, k = 10, resolution = 1)
idxs <- which(clusters < 22)
X <- X[idxs, ]
distances <- as.matrix(distances)
distances <- distances[idxs, idxs]
distances <- as.dist(distances)
clusters <- clusters[idxs]


#--------------------------------------------------------------------------
# Embed the data  ---------------------------------------------------------
#--------------------------------------------------------------------------

# C+E
for (i in 1:5) {
  out <- cluster_align(distances, clusters,
    embedding_method = "PCA", alignment_method = "stress",
    use_grad = T, alpha = 1
  )
  save_embeddings(out$embeddings, labels, cluster_labels = clusters, path = paste0("results/neurons/ce_leiden_run", i, "b.csv"))
}

for (i in 1:5) {
  out <- cluster_align(distances, clusters,
    embedding_method = "PCA", alignment_method = "stress",
    use_grad = T, alpha = 1.75
  )
  save_embeddings(out$embeddings, labels[idxs], cluster_labels = clusters, path = paste0("results/neurons/ce_leiden_175_run", i, ".csv"))
}


# C+E with LOE
for (i in 1:5) {
  embeddings <- clusterwise_loe(clusters, distances)
  save_embeddings(embeddings, labels, cluster_labels = clusters, paste0("results/neurons/loe_11_run", i, ".csv"))
}

# align cluster embeddings to get global embedding
for (i in 1:5) {
  embeddings <- read_csv(paste0("results/neurons/loe_run", i, ".csv"))
  embeddings <- as.matrix(embeddings[, 1:2])
  out <- cluster_align(distances, clusters,
    embedding_method = "NA", alignment_method = "stress",
    use_grad = T, alpha = 1.1, embeddings = embeddings
  )
  save_embeddings(out$embeddings, labels,
    cluster_labels = clusters,
    path = paste0("results/neurons/loe_11_run", i, ".csv")
  )
}


# tsne
for (i in 1:5) {
  out <- get_tsne(X, perplexity = 400)
  save_embeddings(out, labels, paste0("results/neurons/tsne_400_run", i, ".csv"))
}
for (i in 1:5) {
  out <- get_tsne(X, perplexity = 30)
  save_embeddings(out, labels, paste0("results/neurons/tsne_30_run", 1, ".csv"))
}


pca <- get_pca(euclidean_distances)
save_embeddings(pca, labels, path = "results/neurons/pca.csv")

isomap <- get_isomap(distances, "geodesic")
save_embeddings(isomap, labels, path = "results/neurons/isomap.csv")


#--------------------------------------------------------------------------
# Metrics -----------------------------------------------------------------
#--------------------------------------------------------------------------

embedding_paths <- c(
  paste0("results/neurons/ce_leiden_run", 1:5, ".csv"),
  paste0("results/neurons/ce_leiden_175_run", 1:5, ".csv"),
  paste0("results/neurons/tsne_30_run", 1:5, ".csv"),
  paste0("results/neurons/tsne_400_run", 1:5, ".csv"),
  paste0("results/neurons/loe_11_run", 1:5, ".csv"),
  paste0("results/neurons/pca.csv"),
  "results/neurons/isomap.csv"
)

embedding_names <- c(rep("ce_leiden", 5), rep("ce_leiden2", 5), rep("tsne", 5), rep("tsne_2", 5), rep("loe", 5), "pca", "isomap")

results <- get_distance_metrics(embedding_paths, embedding_names, distances, labels)


#--------------------------------------------------------------------------
# Embedding Plots ---------------------------------------------------------
#--------------------------------------------------------------------------

ce1 <- read_csv("results/neurons/ce_leiden_run1.csv") |> mutate(across(3:4, as.factor))
ce2 <- read_csv("results/neurons/ce_leiden_175_run1.csv") |> mutate(across(3:4, as.factor))
loe <- read_csv("results/neurons/loe_11_run1.csv") |> mutate(across(3:4, as.factor))
tsne1 <- read_csv("results/neurons/tsne_30_run1.csv") |> mutate(across(3, as.factor))
tsne2 <- read_csv("results/neurons/tsne_400_run1.csv") |> mutate(across(3, as.factor))
isomap <- read_csv("results/neurons/isomap.csv") |> mutate(across(3, as.factor))
pca <- read_csv("results/neurons/pca.csv") |> mutate(across(3, as.factor))

ce_rotated <- procrustes(X = ce2[, 1:2], Y = ce1[, 1:2], scale = F)$Yrot
ce1[, 1:2] <- ce_rotated

loe_rotated <- procrustes(X = ce2[, 1:2], Y = loe[, 1:2], scale = F)$Yrot
loe[, 1:2] <- loe_rotated

tsne_rotated <- procrustes(X = ce2[, 1:2], Y = tsne1[, 1:2], scale = F)$Yrot
tsne1[, 1:2] <- tsne_rotated

tsne_rotated <- procrustes(X = ce2[, 1:2], Y = tsne2[, 1:2], scale = F)$Yrot
tsne2[, 1:2] <- tsne_rotated

isomap_rotated <- procrustes(X = ce2[, 1:2], Y = isomap[, 1:2], scale = F)$Yrot
isomap[, 1:2] <- isomap_rotated

pca_rotated <- procrustes(X = ce2[, 1:2], Y = pca[, 1:2], scale = F)$Yrot
pca[, 1:2] <- pca_rotated

p1 <- plot(ce1, 1200, 1200, 400, cmap = "rainbow60")
p2 <- plot(ce2, 1500, 1500, 500, cmap = "rainbow60")
p3 <- plot(loe, 1200, 1200, 400, cmap = "rainbow60")
p4 <- plot(tsne1, 60, 60, 20, cmap = "rainbow60")
p5 <- plot(tsne2, 24, 24, 8, cmap = "rainbow60")
p6 <- plot(pca, 60, 60, 20, cmap = "rainbow60")
p7 <- plot(isomap, 800, 800, 200, cmap = "rainbow60")

ggsave("figures/neurons_ce.png", plot = p1, width = 6, height = 6, units = "in", dpi = 100)
ggsave("figures/neurons_ce_175.png", plot = p2, width = 6, height = 6, units = "in", dpi = 100)
ggsave("figures/neurons_loe.png", plot = p3, width = 6, height = 6, units = "in", dpi = 100)
ggsave("figures/neurons_tsne_20.png", plot = p4, width = 6, height = 6, units = "in", dpi = 100)
ggsave("figures/neurons_tsne_400.png", plot = p5, width = 6, height = 6, units = "in", dpi = 100)
ggsave("figures/neurons_pca.png", plot = p6, width = 6, height = 6, units = "in", dpi = 100)
ggsave("figures/neurons_isomap.png", plot = p7, width = 6, height = 6, units = "in", dpi = 100)


#--------------------------------------------------------------------------
# Recall Plots ------------------------------------------------------------
#--------------------------------------------------------------------------

embedding_names <- c(
  paste0("ce1_run", 1:5), paste0("ce2_run", 1:5),
  paste0("tsne30_run", 1:5), paste0("tsne400_run", 1:5),
  paste0("loe_run", 1:5), "pca", "isomap"
)


k_vals <- seq(5, 500, 5)
knn_results <- get_recall_metrics(embedding_paths, embedding_names, distances, k_vals)

ce_alpha1_mean <- rowMeans(knn_results[, 1:5])
ce_alpha1_sd <- apply(knn_results[, 1:5], 1, sd)
ce_alpha175_mean <- rowMeans(knn_results[, 6:10])
ce_alpha175_sd <- apply(knn_results[, 6:10], 1, sd)
tsne30_mean <- rowMeans(knn_results[, 11:15])
tsne30_sd <- apply(knn_results[, 11:15], 1, sd)
tsne400_mean <- rowMeans(knn_results[, 16:20])
tsne400_sd <- apply(knn_results[, 16:20], 1, sd)
ce_loe_mean <- rowMeans(knn_results[, 21:25])
ce_loe_sd <- apply(knn_results[, 21:25], 1, sd)
pca <- knn_results[26]
isomap <- knn_results[27]


knn_summary <- as_tibble(cbind(
  ce_alpha1_mean, ce_alpha1_sd,
  ce_alpha175_mean, ce_alpha175_sd,
  ce_loe_mean, ce_loe_sd,
  tsne30_mean, tsne30_sd,
  tsne400_mean, tsne400_sd, isomap, pca
))
knn_summary$k <- k_vals
write_csv(knn_summary, "results/neurons/knn.csv")
knn_summary <- read_csv("results/neurons/knn.csv")


ce_cols <- c("ce_alpha1", "ce_alpha175", "ce_loe")
ce_labels <- c("C+E (Isomap, α = 1)", "C+E (Isomap, α = 1.75)", "C+E (LOE, α = 1.10)")
tsne_cols <- c("tsne30", "tsne400")
tsne_labels <- c("t-SNE (u = 30)", "t-SNE (u = 400)")

p <- knn_plot(knn_summary, ce_cols, tsne_cols, ce_labels, tsne_labels)
ggsave("figures/neurons_knn.png", plot = p, width = 11, height = 6, units = "in", dpi = 100)
