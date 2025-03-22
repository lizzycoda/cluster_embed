library(ggplot2)
library(dbscan)
library(Rtsne)
library(tidyr)
library(dplyr)

source("./scripts/generate_data.R")
source("./scripts/cluster_align.R")

#--------------------------------------------------------------------------
# Two-dimensional example -------------------------------------------------
#--------------------------------------------------------------------------

# load the data
data <- load_data("shapes_2d", n = 1000)
X <- data$X
labels <- data$labels


# Cluster and embed methods 


# our method
clusters <- dbscan(X, eps = 20, minPts = 8)$cluster # cluster and embed
distances <- dist(X)
out <- cluster_align(distances, clusters, 
                     embedding_method = "MDS", 
                     alignment_method = "stress",
                     use_grad = T)

# tsne 
tsne <- Rtsne(X, dims = 2, perplexity = 200)$Y






#--------------------------------------------------------------------------
# Plot results ------------------------------------------------------------
#--------------------------------------------------------------------------
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
  geom_point( size  = 2) +
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
ggsave("figures/shapes2d_tsne.png", plot = p3, width = 10, height = 10, units = "in", dpi = 100)


# Evaluate -----------------------------------------------------------------




#----------------------------------------------------------------------------
# Three-dimensional example -------------------------------------------------
#----------------------------------------------------------------------------



