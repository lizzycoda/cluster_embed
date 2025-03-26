library(ggplot2)
library(plotly)
library(Rtsne)
library(tidyr)
library(dplyr)
library(readr)
library(vegan)
library(stats)

source("./scripts/generate_data.R")
source("./scripts/cluster_align.R")
source("./scripts/metrics.R")


data <- load_data("halfsphere", n = 1000)
X <- data$X
labels <- data$labels

df = as_tibble(X)
colnames(df) = c('X1', 'X2', 'X3')
df$labels = as.factor(labels)


# our method
clusters <- dbscan(X, eps = .1, minPts = 5)$cluster
distances <- dist(X)
out <- cluster_align(distances, clusters,
                     embedding_method = "MDS",
                     alignment_method = "stress",
                     use_grad = T
)

# Classical MDS
mds <- stats::cmdscale(distances)

# Isomap 
isomap_emb <- isomap(distances, k = 95, ndim = 2)$points #set k large enough so graph is connected

# tSNE 
tsne <- Rtsne(X, dims = 2, perplexity = 200)$Y

# Plot results ------------------------------------------------------------

df <- as_tibble(cbind(X, out$embeddings, tsne, mds, isomap_emb))
colnames(df) <- c("X1", "X2", "X3", "V1", "V2", "Tsne1", "Tsne2",
                  'MDS1', 'MDS2', 'Isomap1', 'Isomap2')
df$clusters <- as.factor(labels) # color by the ground truth clusters 

#save directly from plot window 
gg_colors <- scale_color_discrete()$palette(10)
p1<- plot_ly(df, x = ~X1, y = ~X2, z = ~X3, color = ~clusters, size = 1 ,
              type = "scatter3d", mode = "markers", colors = gg_colors) |>  
  layout(showlegend = FALSE,
         scene = list(
           xaxis = list(title = ""),
           yaxis = list(title = ""),
           zaxis = list(title = "")))

p2 <- ggplot(df, aes(x = V1, y = V2, color = clusters)) +
  geom_point(size = 2) +
  xlim(-1.5,1.5)+
  ylim(-1.75,1.25) + 
  coord_fixed() +
  theme(
    axis.text = element_text(size = 12),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none"
  )

p3 <- ggplot(df, aes(x = Tsne1, y = Tsne2, color = clusters)) +
  geom_point(size = 2) +
  xlim(-10,10)+
  ylim(-10,10) + 
  coord_fixed() +
  theme(
    axis.text = element_text(size = 12),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none"
  ) 

p4 <- ggplot(df, aes(x = MDS1, y = MDS2, color = clusters)) +
  geom_point(size = 2) +
  xlim(-1.5,1.5)+
  ylim(-1.5,1.5) + 
  coord_fixed() +
  theme(
    axis.text = element_text(size = 12),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none"
  )

p5 <- ggplot(df, aes(x = Isomap1, y = Isomap2, color = clusters)) +
  geom_point(size = 2) +
  xlim(-3,3)+
  ylim(-3,3) + 
  coord_fixed() +
  theme(
    axis.text = element_text(size = 12),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none"
  )

ggsave("figures/halfsphere_ce.png", plot = p2, width = 6, height = 6, units = "in", dpi = 100)
ggsave("figures/halfsphere_tsne.png", plot = p3, width = 6, height = 6, units = "in", dpi = 100)
ggsave("figures/halfsphere_mds.png", plot = p4, width = 6, height = 6, units = "in", dpi = 100)
ggsave("figures/halfsphere_isomap.png", plot = p5, width = 6, height = 6, units = "in", dpi = 100)

# Evaluate -----------------------------------------------------------------

emb_distances <- dist(out$embeddings)
tsne_distances <- dist(tsne)
mds_distances<- dist(mds)
isomap_distances <- dist(isomap_emb)
emb_distances_list <- list(emb_distances, tsne_distances, mds_distances, isomap_distances)


# precision-recall plots
eps_vals = c(seq(.1, 2, by = .1), seq(3,18, by = 1))
out <- precision_recall_plots(distances, emb_distances_list, c("Ours", "t-SNE", 'MDS', 'Isomap'), eps_vals)
precision_plot <- out$precision #because of different scales, tsne precision ends up very high 
recall_plot <- out$recall

# k nbhd plot
k_vals <- c(seq(1, 24), seq(25, 1000, by = 25))
knn_plot <- nn_overlap_plot(distances, emb_distances_list, c("Ours", "t-SNE", 'MDS', 'Isomap'), k_vals)

get_metrics(distances, emb_distances)
get_metrics(distances, tsne_distances)
get_metrics(distances, mds_distances)
get_metrics(distances, isomap_distances)
