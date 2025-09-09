# Cluster + Embed 

This repository contains the code to reproduce the experiments and figures in the paper [Cluster and then Embed: A Modular Approach for Visualization](https://arxiv.org/abs/2509.03373). 

<p align="center"><img width="800" alt="Schematic of the cluster + embed approach" src="/figures/ce_schematic.png">

Dimensionality reduction methods such as t-SNE and UMAP are popular methods for visualizing data with a potential (latent) clustered structure. The popularity of these methods can (in part) be explained by their tendency to separate existing clusters well in the visualization. However, while they preserve local information well, they are known to distort global distances. 

We propose a modular approach for visualization of clustered data that consists of first clustering the data, then embedding each cluster, and finally aligning the clusters to obtain a global embedding. 
By embedding each cluster individually, we aim to produce an embedding that suffers from less distortion than if we were to embed all of the data together. 
Additionally, the method includes a tuning parameter to allow for separation between clusters. The method is competitive with existing methods, while offering much more flexibility and transparency. 

Below is an example of the C+E approach on the MNIST dataset (left) and on the developmental human brain organoid data from [Kanton et al. (2019)](https://www.nature.com/articles/s41586-019-1654-9) (right).

<p align="center"><img width="800" alt="Example of the cluster + embed approach on real data" src="/figures/ce_examples.png">


## Example Usage

Load the source files

``` 
source("./scripts/generate_data.R")
source("./scripts/cluster_embed.R")
source("./scripts/utils.R")
``` 

Synthetic data can be loaded using 

```
data <- load_data("shapes_3d", n = 1000)
X <- data$X
labels <- data$labels
euclidean_distances <- dist(X) 
``` 

Any clustering method can then be used to cluster the data. In this example, we use DBSCAN. The embedding is obtained using the `cluster_align` function. In this example each cluster is embedded via PCA and the clusters are aligned to minimize Kruskal's formulation of the stress.

```
clusters <-  dbscan(X, eps = .1, minPts = 10)$cluster
out<- cluster_align(euclidean_distances, clusters,
    embedding_method = "PCA",
    alignment_method = "stress",
    use_grad = T
  )
ce_embedding <- out$embeddings
```

With real world data, Euclidean distance may not capture the intrinsic dissimilarity between two observations well and other notions of dissimilarity may be more appropriate to work with. Below we calcualte the geodesic distance. It is recommended to save these distances when the dataset is large as initial distance calculation can be slow. 

```
geodesic_distances <- calculate_geodesic_distances(euclidean_distances,
  k = 5, save = TRUE,
  path = "path_to_data"
)

geodesic_distances <- load_geodesic_distances("path_to_data", n = 1000)

```

Then, the C+E approach can be applied to these distances. In our real-world examples we find that the Leiden algorithm performs well for clustering the data. We then embed each cluster and align the clusters. Note that when we input the geodesic distances to the `cluster_align` function and set the embedding method to PCA, this is equivalent to embedding each cluster via Isomap.
```
clusters<- get_leiden_clustering(distances, k = 100, resolution = 1 / 2)
out<- cluster_align(geodesic_distances, clusters,
    embedding_method = "PCA",
    alignment_method = "stress",
    use_grad = T
  )
ce_embedding <- out$embeddings
```

## Files

All experiments in the paper and its appendix are in the examples folder, with a separate file for each dataset.

Each file contains code to obtain the C+E embedding and to evaluate the embedding on a variety of different metrics. 

The repository contains the following files:

- `scripts/`
  - `cluster_embed.R` contains the main code for the method.
  - `generate_data.R` generates the data in our synthetic examples.
  - `utils.R` contains useful embedding and clustering functions.
  - `eval.R` contains the code to evaluate embeddings on the metrics used in the paper, as well as produce plots.
  - `loe_embedding.R` embeds data using a gradient descent approach to minimize the Local Ordinal Embedding objective in [Terada and Von Luxburg (2014)](https://proceedings.mlr.press/v32/terada14.pdf).
  - `stress_functions.R` contains helper functions for the main method.
 

- `examples/` contains the examples in the paper, listed below.
  - `gaussian_mixture.R`
  - `shapes1.R`
  - `shapes2.R`
  - `mnist.R`
  - `human-409b2.R`
  - `fmnist.R`
  - `neurons.R`


## Requirements

The code is built in R (4.4.2) and uses the following packages:
- ggplot2 (3.5.1)
- readr (2.5.1)
- dplyr (1.1.4)
- tidyr (1.3.1)
- Rtsne (0.17)
- stats (4.4.2)
- MASS (7.3.61)
- vegan (2.6.10)
- FreeSortR (1.3)
- igraph (2.1.3)
- dbscan (1.2.2)
- proxy (0.4.27)
- torch (0.14.2) (Only used to obtain the LOE embedding)


## Citation 

```
@article{cluster_embed25,
      title={Cluster and then Embed: A Modular Approach for Visualization}, 
      author={Elizabeth Coda and Ery Arias-Castro and Gal Mishne},
      journal={arXiv preprint arXiv: 2509.03373}
      year={2025}
}
```