source("./scripts/generate_data.R")
source("./scripts/cluster_embed.R")
source("./scripts/eval.R")
source("./scripts/utils.R")
source("./scripts/loe_embedding.R")

#--------------------------------------------------------------------------
# Data loading and preprocessing ------------------------------------------
#--------------------------------------------------------------------------

# load the data
# dataset loaded using this reference:
# https://github.com/berenslab/ne_spectrum_scRNAseq/blob/982f619c0b84a67999e4a14e3c016c86705da270/utils/utils.py#L99
# this dataset was then subsampled from 20272 samples to 5k samples
data <- read_csv("data/bio_data_5k.csv")

X <- as.matrix(data[1:50]) # preprocessing has already been done
labels <- pull(data, "clusters") # 0 through 6 labeled

# calculate geodesic distances
euclidean_distances <- dist(X)
calculate_geodesic_distances(euclidean_distances,
  k = 5, save = TRUE,
  path = "data/bio_geodesic_5.csv"
)

# just load after run above function
distances <- load_geodesic_distances("data/bio_geodesic_5.csv", 5000)

# precompute and save clusters
for (i in 1:5) {
  clusters <- get_leiden_clustering(distances, k = 100, resolution = 1 / 2)
  print(rand_index(labels, clusters))
  save_clusters(clusters, paste0("results/bio/leiden_clusters", i, ".csv"))
}


#--------------------------------------------------------------------------
# Embed the data  ---------------------------------------------------------
#--------------------------------------------------------------------------

# C+E
for (i in 1:5) {
  clusters <- load_clusters(paste0("results/bio/leiden_clusters", i, ".csv"))
  out <- cluster_align(distances, clusters,
    embedding_method = "PCA", alignment_method = "stress",
    use_grad = T, alpha = 1
  )
  save_embeddings(out$embeddings, labels, cluster_labels = clusters, path = paste0("results/bio/ce_leiden_run", i, ".csv"))
}

for (i in 1:5) {
  print(i)
  clusters <- load_clusters(paste0("results/bio/leiden_clusters", i, ".csv"))
  out <- cluster_align(distances, clusters,
    embedding_method = "PCA", alignment_method = "stress",
    use_grad = T, alpha = 1.15
  )
  save_embeddings(out$embeddings, labels, cluster_labels = clusters, path = paste0("results/bio/ce_leiden_115_run", i, ".csv"))
}

# C+E with LOE
for (i in 1:5) {
  clusters <- load_clusters(paste0("results/bio/leiden_clusters", i, ".csv"))
  embeddings <- clusterwise_loe(clusters, distances)
  save_embeddings(embeddings, labels, cluster_labels = clusters, paste0("results/bio/loe_run", i, ".csv"))
}

# align cluster embeddings to get global embedding
for (i in 1:5) {
  clusters <- load_clusters(paste0("results/bio/leiden_clusters", i, ".csv"))
  embeddings <- read_csv(paste0("results/bio/loe_run", i, ".csv"))
  embeddings <- as.matrix(embeddings[, 1:2])
  out <- cluster_align(distances, clusters,
    embedding_method = "NA", alignment_method = "stress",
    use_grad = T, alpha = 1.75, embeddings = embeddings
  )
  save_embeddings(out$embeddings, labels, cluster_labels = clusters, path = paste0("results/bio/ce_leiden_loe_175_run", i, ".csv"))
}



# tuning of tsne parameter
scores <- c()
vals <- seq(50, 1000, by = 50)
for (p in vals) {
  out <- get_tsne(X, perplexity = p)
  spearman <- get_spearman(distances, dist(out))
  scores <- c(scores, spearman)
}
perplexity <- vals[which.max(scores)] # 500

for (i in 1:5) {
  out <- get_tsne(X, perplexity = 500)
  save_embeddings(out, labels, paste0("results/bio/tsne_500_run", i, ".csv"))
}

# good local behavior
for (i in 1:5) {
  out <- get_tsne(X, perplexity = 20)
  save_embeddings(out, labels, paste0("results/bio/tsne_20_run", i, ".csv"))
}


pca <- get_pca(euclidean_distances)
save_embeddings(pca, labels, path = "results/bio/pca.csv")

isomap <- get_isomap(distances, "geodesic")
save_embeddings(isomap, labels, path = "results/bio/isomap.csv")



#--------------------------------------------------------------------------
# Metrics -----------------------------------------------------------------
#--------------------------------------------------------------------------

embedding_paths <- c(
  paste0("results/bio/ce_leiden_run", 1:5, ".csv"),
  paste0("results/bio/ce_leiden_115_run", 1:5, ".csv"),
  paste0("results/bio/tsne_20_run", 1:5, ".csv"),
  paste0("results/bio/tsne_500_run", 1:5, ".csv"),
  paste0("results/bio/ce_leiden_loe_175_run", 1:5, ".csv"),
  paste0("results/bio/pca.csv"),
  "results/bio/isomap.csv"
)

embedding_names <- c(rep("ce_leiden", 5), rep("ce_leiden2", 5), rep("tsne", 5), rep("tsne_2", 5), rep("loe", 5), "pca", "isomap")

results <- get_distance_metrics(embedding_paths, embedding_names, distances, labels)
results2 <- get_distance_metrics(embedding_paths, embedding_names, euclidean_distances, labels)


#--------------------------------------------------------------------------
# Embedding Plots ---------------------------------------------------------
#--------------------------------------------------------------------------

ce1 <- read_csv("results/bio/ce_leiden_run1.csv") |> mutate(across(3:4, as.factor))
ce2 <- read_csv("results/bio/ce_leiden_115_run1.csv") |> mutate(across(3:4, as.factor))
loe <- read_csv("results/bio/ce_leiden_loe_175_run1.csv") |> mutate(across(3:4, as.factor))
tsne1 <- read_csv("results/bio/tsne_20_run1.csv") |> mutate(across(3, as.factor))
tsne2 <- read_csv("results/bio/tsne_500_run1.csv") |> mutate(across(3, as.factor))
pca <- read_csv("results/bio/pca.csv") |> mutate(across(3, as.factor))
isomap <- read_csv("results/bio/isomap.csv") |> mutate(across(3, as.factor))


# alignment of outputs
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

p1 <- plot(ce1, 150, 150, 50, cmap = "rainbow6")
p2 <- plot(ce2, 150, 150, 50, cmap = "rainbow6")
p3 <- plot(loe, 225, 225, 75, cmap = "rainbow6")
p4 <- plot(tsne1, 60, 60, 20, cmap = "rainbow6")
p5 <- plot(tsne2, 21, 21, 7, cmap = "rainbow6")
p6 <- plot(pca, 21, 21, 7, cmap = "rainbow6")
p7 <- plot(isomap, 120, 120, 40, cmap = "rainbow6")

ggsave("figures/bio_ce.png", plot = p1, width = 6, height = 6, units = "in", dpi = 100)
ggsave("figures/bio_ce_115.png", plot = p2, width = 6, height = 6, units = "in", dpi = 100)
ggsave("figures/bio_ce_loe.png", plot = p3, width = 6, height = 6, units = "in", dpi = 100)
ggsave("figures/bio_tsne_20.png", plot = p4, width = 6, height = 6, units = "in", dpi = 100)
ggsave("figures/bio_tsne_500.png", plot = p5, width = 6, height = 6, units = "in", dpi = 100)
ggsave("figures/bio_pca.png", plot = p6, width = 6, height = 6, units = "in", dpi = 100)
ggsave("figures/bio_isomap.png", plot = p7, width = 6, height = 6, units = "in", dpi = 100)

#--------------------------------------------------------------------------
# Recall Plots ------------------------------------------------------------
#--------------------------------------------------------------------------

embedding_names <- c(
  paste0("ce1_run", 1:5), paste0("ce115_run", 1:5),
  paste0("tsne20_run", 1:5), paste0("tsne500_run", 1:5),
  paste0("ce_loe_run", 1:5), "pca", "isomap"
)

k_vals <- seq(5, 500, 5)
knn_results <- get_recall_metrics(embedding_paths, embedding_names, distances, k_vals) 

ce_alpha1_mean <- rowMeans(knn_results[, 1:5])
ce_alpha1_sd <- apply(knn_results[, 1:5], 1, sd)
ce_alpha115_mean <- rowMeans(knn_results[, 6:10])
ce_alpha115_sd <- apply(knn_results[, 6:10], 1, sd)
tsne20_mean <- rowMeans(knn_results[, 11:15])
tsne20_sd <- apply(knn_results[, 11:15], 1, sd)
tsne500_mean <- rowMeans(knn_results[, 16:20])
tsne500_sd <- apply(knn_results[, 16:20], 1, sd)
ce_loe_mean <- rowMeans(knn_results[, 21:25])
ce_loe_sd <- apply(knn_results[, 21:25], 1, sd)
pca <- knn_results[26]
isomap <- knn_results[27]

knn_summary <- as_tibble(cbind(
  ce_alpha1_mean, ce_alpha1_sd,
  ce_alpha115_mean, ce_alpha115_sd,
  ce_loe_mean, ce_loe_sd,
  tsne20_mean, tsne20_sd,
  tsne500_mean, tsne500_sd, pca, isomap
))
knn_summary$k <- k_vals


write_csv(knn_summary, "results/bio/knn_results.csv")
knn_summary <- read_csv("results/bio/knn_results.csv")


ce_cols <- c("ce_alpha1", "ce_alpha115", "ce_loe")
ce_labels <- c("C+E (Isomap, α = 1)", "C+E (Isomap, α = 1.15)", "C+E (LOE, α = 1.75)")
tsne_cols <- c("tsne20", "tsne500")
tsne_labels <- c("t-SNE (u = 20)", "t-SNE (u = 500)")


p <- knn_plot(knn_summary, ce_cols, tsne_cols, ce_labels, tsne_labels)
ggsave("figures/bio_knn.png", plot = p, width = 10, height = 6, units = "in", dpi = 100)
