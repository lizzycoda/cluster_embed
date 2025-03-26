library(ggplot2)
library(plotly)
library(Rtsne)
library(tidyr)
library(dplyr)
library(readr)

source("./scripts/generate_data.R")
source("./scripts/cluster_align.R")
source("./scripts/metrics.R")

#--------------------------------------------------------------------------
# Two-dimensional example -------------------------------------------------
#--------------------------------------------------------------------------

# load the data
data <- load_data("shapes_2d", n = 1000)
X <- data$X
labels <- data$labels


# Cluster and embed

# our method
clusters <- dbscan(X, eps = 20, minPts = 8)$cluster # cluster and embed
distances <- dist(X)
out <- cluster_align(distances, clusters,
  embedding_method = "MDS",
  alignment_method = "stress",
  use_grad = T
)

# tsne
tsne <- Rtsne(X, dims = 2, perplexity = 200)$Y


# Plot results ------------------------------------------------------------

df <- as_tibble(cbind(X, out$embeddings, tsne))
colnames(df) <- c("X1", "X2", "V1", "V2", "Tsne1", "Tsne2")
df$clusters <- as.factor(labels)

p1 <- ggplot(df, aes(x = X1, y = X2, color = clusters)) +
  geom_point(size = 2) +
  xlim(-300, 300) +
  ylim(-300, 300) +
  coord_fixed() +
  theme(
    axis.text = element_text(size = 12),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none"
  )

p2 <- ggplot(df, aes(x = V1, y = V2, color = clusters)) +
  geom_point(size = 2) +
  xlim(-300, 300) +
  ylim(-300, 300) +
  coord_fixed() +
  theme(
    axis.text = element_text(size = 12),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none"
  )

p3 <- ggplot(df, aes(x = Tsne1, y = Tsne2, color = clusters)) +
  geom_point(size = 2) +
  xlim(-12, 12) +
  ylim(-12, 12) +
  coord_fixed() +
  theme(
    axis.text = element_text(size = 12),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none"
  )

ggsave("figures/shapes2d_original.png", plot = p1, width = 6, height = 6, units = "in", dpi = 100)
ggsave("figures/shapes2d_ce.png", plot = p2, width = 6, height = 6, units = "in", dpi = 100)
ggsave("figures/shapes2d_tsne.png", plot = p3, width = 6, height = 6, units = "in", dpi = 100)


# Evaluate -----------------------------------------------------------------

# distances of each embedding method
emb_distances <- dist(out$embeddings)
tsne_distances <- dist(tsne)
emb_distances_list <- list(emb_distances, tsne_distances)


# precision-recall plots
out <- precision_recall_plots(distances, emb_distances_list, c("Ours", "t-SNE"))
precision_plot <- out$precision
recall_plot <- out$recall

# k nbhd plot
k_vals <- c(seq(1, 24), seq(25, 1000, by = 25))
knn_plot <- nn_overlap_plot(distances, emb_distances_list, c("Ours", "t-SNE"), k_vals)


# calculate metrics over several runs
results <- matrix(NA, nrow = 10, ncol = 6)
for (i in 1:10) {
  embeddings <- cluster_align(distances, clusters,
    embedding_method = "MDS",
    alignment_method = "stress",
    use_grad = T
  )$embeddings
  tsne_embeddings <- Rtsne(X, dims = 2, perplexity = 200)$Y

  emb_distances <- dist(embeddings)
  tsne_distances <- dist(tsne_embeddings)

  results[i, 1:3] <- unlist(get_metrics(distances, emb_distances))
  results[i, 4:6] <- unlist(get_metrics(distances, tsne_distances))
}

results_df <- as_tibble(results)
colnames(results_df) <- c(
  "stress-ours", "normalized_stress-ours", "spearman-ours",
  "stress-tsne", "normalized_stress-tsne", "spearman-tsne"
)
write_csv(results_df, "shapes2d_distance_metrics.csv") # save results


# look and mean and sd of metrics over several runs
summarize_results(results_df)


#----------------------------------------------------------------------------
# Three-dimensional example -------------------------------------------------
#----------------------------------------------------------------------------

data <- load_data("shapes_3d", n = 1000)
X <- data$X
labels <- data$labels

# Cluster and embed

# our method
clusters <- dbscan(X, eps = 20, minPts = 8)$cluster # cluster and embed
distances <- dist(X)
out <- cluster_align(distances, clusters,
  embedding_method = "MDS",
  alignment_method = "stress",
  use_grad = T
)

# tsne
tsne <- Rtsne(X, dims = 2, perplexity = 50)$Y

# Plot results ------------------------------------------------------------

df <- as_tibble(cbind(X, out$embeddings, tsne))
colnames(df) <- c("X1", "X2", "X3", "V1", "V2", "Tsne1", "Tsne2")
df$clusters <- as.factor(labels) # color by the ground truth clusters 

# save this directly from viewer
gg_colors <- scale_color_discrete()$palette(3) # use same colors as with ggplot
p1 <- plot_ly(df,
  x = ~X1, y = ~X2, z = ~X3, color = ~clusters, size = 1,
  colors = gg_colors, type = "scatter3d", mode = "markers"
) |>
  hide_colorbar()

p2 <- ggplot(df, aes(x = V1, y = V2, color = clusters)) +
  xlim(-100, 250) +
  ylim(-125, 225) +
  geom_point(size = 2) +
  coord_fixed() +
  theme(
    axis.text = element_text(size = 12),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none"
  )

p3 <- ggplot(df, aes(x = Tsne1, y = Tsne2, color = clusters)) +
  geom_point(size = 2) +
  xlim(-30, 30) +
  ylim(-30, 30) +
  coord_fixed() +
  theme(
    axis.text = element_text(size = 12),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none"
  )


ggsave("figures/shapes3d_ce.png", plot = p2, width = 6, height = 6, units = "in", dpi = 100)
ggsave("figures/shapes3d_tsne.png", plot = p3, width = 6, height = 6, units = "in", dpi = 100)


# Evaluate -----------------------------------------------------------------

# distances of each embedding method
emb_distances <- dist(out$embeddings)
tsne_distances <- dist(tsne)
emb_distances_list <- list(emb_distances, tsne_distances)


# precision-recall plots
out <- precision_recall_plots(distances, emb_distances_list, c("Ours", "t-SNE"))
precision_plot <- out$precision
recall_plot <- out$recall

# k nbhd plot
k_vals <- c(seq(1, 24), seq(25, 1000, by = 25))
knn_plot <- nn_overlap_plot(distances, emb_distances_list, c("Ours", "t-SNE"), k_vals)

get_metrics(distances, emb_distances)
get_metrics(distances, tsne_distances)

# TO DO: repeat several runs to get variability of some metrics
