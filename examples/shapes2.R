source("./scripts/generate_data.R")
source("./scripts/cluster_embed.R")
source("./scripts/eval.R")
source("./scripts/utils.R")
library(scatterplot3d)

#----------------------------------------------------------------------------
# 3D example -------------------------------------------------------
#----------------------------------------------------------------------------
data <- load_data("shapes_3d", n = 1000)
X <- data$X

labels <- data$labels

distances <- dist(X) # use euclidean distances
clusters <- dbscan(X, eps = .1, minPts = 10)$cluster # this clustering method is deterministic
table(clusters)


#--------------------------------------------------------------------------
# Embed the data  ---------------------------------------------------------
#--------------------------------------------------------------------------


# our method with no scaling
for (i in 1:5) {
  out <- cluster_align(distances, clusters,
    embedding_method = "PCA",
    alignment_method = "stress",
    use_grad = T
  )
  save_embeddings(out$embeddings, labels, cluster_labels = clusters, path = paste0("results/shapes2/ce_run", i, ".csv"))
}


for (i in 1:5) {
  out <- cluster_align(distances, clusters,
    embedding_method = "PCA",
    alignment_method = "stress",
    use_grad = T, alpha = 1.42
  )
  save_embeddings(out$embeddings, labels, cluster_labels = clusters, path = paste0("results/shapes2/ce_alpha142_run", i, ".csv"))
}


# tuning of tsne parameter
scores <- c()
vals <- seq(10, 330, by = 10)
for (p in vals) {
  out <- get_tsne(X, perplexity = p)
  spearman <- get_spearman(distances, dist(out))
  scores <- c(scores, spearman)
}
perplexity <- vals[which.max(scores)] # 260

# tsne
for (i in 1:5) {
  out <- get_tsne(X, perplexity = perplexity)
  save_embeddings(out, labels, paste0("results/shapes2/tsne_run", i, ".csv"))
}

# tsne
for (i in 1:5) {
  out <- get_tsne(X, perplexity = 50)
  save_embeddings(out, labels, paste0("results/shapes2/tsne50_run", i, ".csv"))
}


# tuning of isomap parameter
scores <- c()
vals <- seq(175, 500, by = 25)
for (k in vals) {
  out <- get_isomap(distances, "euclidean", k = k)
  spearman <- get_spearman(distances, dist(out))
  scores <- c(scores, spearman)
}
k <- vals[which.max(scores)] # k=500

# isomap (is deterministic)
out <- get_isomap(distances, "euclidean", k = k)
save_embeddings(out, labels, paste0("results/shapes2/isomap.csv"))

# mds (is deterministic)
pca <- get_pca(distances)
save_embeddings(pca, labels, "results/shapes2/pca.csv")


#--------------------------------------------------------------------------
# Metrics -----------------------------------------------------------------
#--------------------------------------------------------------------------

embedding_paths <- c(
  paste0("results/shapes2/ce_run", 1:5, ".csv"),
  paste0("results/shapes2/ce_alpha142_run", 1:5, ".csv"),
  paste0("results/shapes2/tsne_run", 1:5, ".csv"),
  paste0("results/shapes2/tsne50_run", 1:5, ".csv"),
  "results/shapes2/isomap.csv",
  "results/shapes2/pca.csv"
)

embedding_names <- c(rep("ce", 5), rep("ce2", 5), rep("tsne", 5), rep("tsne50", 5), "isomap", "pca")
results <- get_distance_metrics(embedding_paths, embedding_names, distances, labels)



#--------------------------------------------------------------------------
# Embedding Plots ---------------------------------------------------------
#--------------------------------------------------------------------------

ce1 <- read_csv("results/shapes2/ce_run1.csv") |> mutate(across(3:4, as.factor))
ce2 <- read_csv("results/shapes2/ce_alpha142_run1.csv") |> mutate(across(3:4, as.factor))
tsne <- read_csv("results/shapes2/tsne50_run1.csv") |> mutate(across(3, as.factor))
tsne2 <- read_csv("results/shapes2/tsne_run1.csv") |> mutate(across(3, as.factor))
pca <- read_csv("results/shapes2/pca.csv") |> mutate(across(3, as.factor))
isomap <- read_csv("results/shapes2/isomap.csv") |> mutate(across(3, as.factor))

# align outputs for visualization purposes only
ce_rotated <- procrustes(X = X[, 1:2], Y = ce1[, 1:2], scale = F)$Yrot
ce1[, 1:2] <- ce_rotated
tsne_rotated <- procrustes(X = ce1[, 1:2], Y = tsne[, 1:2], scale = F)$Yrot
tsne[, 1:2] <- tsne_rotated
tsne_rotated2 <- procrustes(X = ce1[, 1:2], Y = tsne2[, 1:2], scale = F)$Yrot
tsne2[, 1:2] <- tsne_rotated2
ce_rotated2 <- procrustes(X = ce1[, 1:2], Y = ce2[, 1:2], scale = F)$Yrot
ce2[, 1:2] <- ce_rotated2
isomap_rotated <- procrustes(X = ce1[, 1:2], Y = isomap[, 1:2], scale = F)$Yrot
isomap[, 1:2] <- isomap_rotated
pca_rotated <- procrustes(X = ce1[, 1:2], Y = pca[, 1:2], scale = F)$Yrot
pca[, 1:2] <- pca_rotated

p1 <- plot(ce1, 1.5, 1.5, cmap = "shapesII")
p2 <- plot(tsne, 30, 30, stepsize = 10, cmap = "shapesII")
p3 <- plot(tsne2, 12, 12, stepsize = 4, cmap = "shapesII")
p4 <- plot(ce2, 1.5, 1.5, cmap = "shapesII")
p5 <- plot(pca, 1.5, 1.5, stepsize = .5, cmap = "shapesII")
p6 <- plot(isomap, 1.5, 1.5, stepsize = .5, cmap = "shapesII")


ggsave("figures/shapes2_ce.png", plot = p1, width = 6, height = 6, units = "in", dpi = 100)
ggsave("figures/shapes2_tsne.png", plot = p2, width = 6, height = 6, units = "in", dpi = 100)
ggsave("figures/shapes2_ce_142.png", plot = p4, width = 6, height = 6, units = "in", dpi = 100)
ggsave("figures/shapes2_tsne_260.png", plot = p3, width = 6, height = 6, units = "in", dpi = 100)
ggsave("figures/shapes2_pca.png", plot = p5, width = 6, height = 6, units = "in", dpi = 100)
ggsave("figures/shapes2_isomap.png", plot = p6, width = 6, height = 6, units = "in", dpi = 100)


# plot original data
colors <- c("#17BECF", "#E377C2", "#BCBD22")
colors <- colors[as.numeric(labels)]

png("figures/shapes2_original.png", width = 6, height = 6, units = "in", res = 100)
par(mar = c(0, 0, 0, 0))
p0 <- scatterplot3d(
  x = X[, 1:3], pch = 16, box = FALSE, color = colors, tick.marks = TRUE,
  xlab = "", ylab = "", zlab = "",
  xlim = c(-1, 1),
  ylim = c(-1, 1),
  zlim = c(-1, 1),
  cex.lab = 1,
  cex.axis = 1,
  cex.symbols = 1
)
dev.off()

#--------------------------------------------------------------------------
# Recall Plots ------------------------------------------------------------
#--------------------------------------------------------------------------
embedding_names <- c(
  paste0("ce1_run", 1:5), paste0("ce142_run", 1:5),
  paste0("tsne260_run", 1:5), paste0("tsne50_run", 1:5), "isomap", "pca"
)

k_vals <- seq(2, 200, 5)
knn_results <- get_recall_metrics(embedding_paths, embedding_names, distances, k_vals)
ce_alpha1_mean <- rowMeans(knn_results[, 1:5])
ce_alpha1_sd <- apply(knn_results[, 1:5], 1, sd)
ce_alpha142_mean <- rowMeans(knn_results[, 6:10])
ce_alpha142_sd <- apply(knn_results[, 6:10], 1, sd)
tsne260_mean <- rowMeans(knn_results[, 11:15])
tsne260_sd <- apply(knn_results[, 11:15], 1, sd)
tsne50_mean <- rowMeans(knn_results[, 16:20])
tsne50_sd <- apply(knn_results[, 16:20], 1, sd)
isomap <- knn_results[21]
pca <- knn_results[22]
knn_summary <- as_tibble(cbind(
  ce_alpha1_mean, ce_alpha1_sd,
  ce_alpha142_mean, ce_alpha142_sd,
  tsne260_mean, tsne260_sd,
  tsne50_mean, tsne50_sd, isomap, pca
))
knn_summary$k <- k_vals
write_csv(knn_summary, "results/shapes2/knn_results.csv")

knn_summary <- read_csv("results/shapes2/knn_results.csv")

ce_cols <- c("ce_alpha1", "ce_alpha142")
ce_labels <- c("C+E (α = 1)", "C+E (α = 1.42)")
tsne_cols <- c("tsne50", "tsne260")
tsne_labels <- c("t-SNE (u = 50)", "t-SNE (u = 260)")


p <- knn_plot(knn_summary, ce_cols, tsne_cols, ce_labels, tsne_labels)
ggsave("figures/shapes2_knn.png", plot = p, width = 10, height = 6, units = "in", dpi = 100)
