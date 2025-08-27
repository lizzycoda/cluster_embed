source("./scripts/generate_data.R")
source("./scripts/cluster_embed.R")
source("./scripts/eval.R")
source("./scripts/utils.R")
library(scatterplot3d)

#----------------------------------------------------------------------------
# Cylindrical example -------------------------------------------------------
#----------------------------------------------------------------------------
data <- load_data("shapes_cylinder", n = 1000)
X <- data$X
labels <- data$labels
distances <- dist(X) # use euclidean distances
clusters <- dbscan(X, eps = .065, minPts = 3)$cluster # this clustering method is deterministic
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
  save_embeddings(out$embeddings, labels, cluster_labels = clusters, path = paste0("results/shapes1/ce_run", i, ".csv"))
}


# tuning of tsne parameter
scores <- c()
vals <- seq(10, 330, by = 10)
for (p in vals) {
  out <- get_tsne(X, perplexity = p)
  spearman <- get_spearman(distances, dist(out))
  scores <- c(scores, spearman)
}
perplexity <- vals[which.max(scores)] # 230

# tsne
for (i in 1:5) {
  out <- get_tsne(X, perplexity = perplexity)
  save_embeddings(out, labels, paste0("results/shapes1/tsne_run", i, ".csv"))
}


# tuning of isomap parameter
scores <- c()
vals <- seq(25, 500, by = 25)
for (k in vals) {
  out <- get_isomap(distances, "euclidean", k = k)
  spearman <- get_spearman(distances, dist(out))
  scores <- c(scores, spearman)
}
k <- vals[which.max(scores)] # k=400

# isomap (is deterministic)
out <- get_isomap(distances, "euclidean", k = k)
save_embeddings(out, labels, paste0("results/shapes1/isomap.csv"))

# pca (is deterministic)
pca <- get_pca(distances)
save_embeddings(pca, labels, "results/shapes1/pca.csv")


#--------------------------------------------------------------------------
# Metrics -----------------------------------------------------------------
#--------------------------------------------------------------------------

embedding_paths <- c(
  paste0("results/shapes1/ce_run", 1:5, ".csv"),
  paste0("results/shapes1/tsne_run", 1:5, ".csv"),
  "results/shapes1/isomap.csv",
  "results/shapes1/pca.csv"
)
embedding_names <- c(rep("ce", 5), rep("tsne", 5), "isomap", "pca")
results <- get_distance_metrics(embedding_paths, embedding_names, distances, labels)

#--------------------------------------------------------------------------
# Embedding Plots ---------------------------------------------------------
#--------------------------------------------------------------------------

ce1 <- read_csv("results/shapes1/ce_run1.csv") |> mutate(across(3:4, as.factor))
tsne <- read_csv("results/shapes1/tsne_run1.csv") |> mutate(across(3, as.factor))
pca <- read_csv("results/shapes1/pca.csv") |> mutate(across(3, as.factor))
isomap <- read_csv("results/shapes1/isomap.csv") |> mutate(across(3, as.factor))


# alignment of outputs
tsne_rotated <- procrustes(X = ce1[, 1:2], Y = tsne[, 1:2], scale = F)$Yrot
tsne[, 1:2] <- tsne_rotated

isomap_rotated <- procrustes(X = ce1[, 1:2], Y = isomap[, 1:2], scale = F)$Yrot
isomap[, 1:2] <- isomap_rotated

pca_rotated <- procrustes(X = ce1[, 1:2], Y = pca[, 1:2], scale = F)$Yrot
pca[, 1:2] <- pca_rotated


p1 <- plot(ce1, 1, 1, cmap = "shapesI")
p2 <- plot(tsne, 10, 10, stepsize = 5, cmap = "shapesI")
p3 <- plot(pca, 1, 1, stepsize = .5, cmap = "shapesI")
p4 <- plot(isomap, 1, 1, stepsize = .5, cmap = "shapesI")

ggsave("figures/shapes1_ce.png", plot = p1, width = 6, height = 6, units = "in", dpi = 100)
ggsave("figures/shapes1_tsne.png", plot = p2, width = 6, height = 6, units = "in", dpi = 100)
ggsave("figures/shapes1_pca.png", plot = p3, width = 6, height = 6, units = "in", dpi = 100)
ggsave("figures/shapes1_isomap.png", plot = p4, width = 6, height = 6, units = "in", dpi = 100)



# also plot original data in R3
colors <- c("#1F77B4", "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF")
colors <- colors[as.numeric(labels)]


# plot a cylinder underneath the points a
theta <- seq(0, pi, length.out = 20)
r <- 2 / pi
x_pts <- r * cos(theta)
y_pts <- r * sin(theta)
z_pts <- seq(min(X[, 3]), max(X[, 3]), length.out = 5)

png("figures/shapes1_original.png", width = 6, height = 6, units = "in", res = 100)
par(mar = c(0, 0, 0, 0))
p0 <- scatterplot3d(
  x = -X[, 1], X[, 2], X[, 3], pch = 16, box = FALSE, color = colors, tick.marks = TRUE,
  xlab = "", ylab = "", zlab = "",
  xlim = c(-.8, .8),
  ylim = c(-.8, .8),
  zlim = c(-.8, .8),
  cex.lab = 1,
  cex.axis = 1,
  cex.symbols = 1,
  angle = 45
)

for (i in 1:length(x_pts)) {
  p0$points3d(x_pts, y_pts, rep(z_pts[i], length(x_pts)), type = "l", col = "gray", lwd = .1)
  p0$points3d(rep(x_pts[i], length(z_pts)), rep(y_pts[i], length(z_pts)),
    z_pts,
    type = "l", col = "gray", lwd = .1
  )
}
p0$points3d(x = -X[, 1], X[, 2], X[, 3], pch = 16, col = colors)
dev.off()



#--------------------------------------------------------------------------
# Recall Plots ------------------------------------------------------------
#--------------------------------------------------------------------------


k_vals <- seq(2, 250, 5)
embedding_names <- c(
  paste0("ce1_run", 1:5),
  paste0("tsne_run", 1:5), "isomap", "pca"
)
knn_results <- get_recall_metrics(embedding_paths, embedding_names, distances, k_vals)
# write_csv(knn_results, 'results/shapes1/knn_results.csv')

knn_results <- read_csv("results/shapes1/knn_results.csv")


ce_mean <- rowMeans(knn_results[, 1:5])
ce_sd <- apply(knn_results[, 1:5], 1, sd)
tsne_mean <- rowMeans(knn_results[, 6:10])
tsne_sd <- apply(knn_results[, 6:10], 1, sd)
iso <- knn_results[11]
pca <- knn_results[12]
knn_summary <- as_tibble(cbind(ce_mean, ce_sd, tsne_mean, tsne_sd, iso, pca))
knn_summary$k <- k_vals

ce_cols <- c("ce")
ce_labels <- c("C+E (α = 1)")
tsne_cols <- c("tsne")
tsne_labels <- c("t-SNE (u = 230)")

p <- knn_plot(knn_summary, ce_cols, tsne_cols, ce_labels, tsne_labels)
ggsave("figures/shapes1_knn.png", plot = p, width = 10, height = 6, units = "in", dpi = 100)
