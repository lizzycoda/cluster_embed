library(ggplot2)
library(plotly)
library(Rtsne)
library(tidyr)
library(dplyr)
library(readr)

source("./scripts/generate_data.R")
source("./scripts/cluster_align.R")
source("./scripts/metrics.R")


data <- load_data("halfsphere", n = 1000)
X <- data$X
labels <- data$labels


df = as_tibble(X)
colnames(df) = c('X1', 'X2', 'X3')
df$labels = as.factor(labels)

p1 <- plot_ly(df, x = ~X1, y = ~X2, z = ~X3, color = ~labels, size = 1 ) |>  layout(showlegend = FALSE)
p1


