# Cluster + Embed 

This repository contains the code to reproduce the experiments and figures in the paper [Cluster and then Embed: A Modular Approach for Visualization](https://arxiv.org/abs/2509.03373). 

<p align="center"><img width="800" alt="Schematic of the cluster + embed approach" src="/figures/ce_schematic.png">

Dimensionality reduction methods such as t-SNE and UMAP are popular methods for visualizing data with a potential (latent) clustered structure. The popularity of these methods can (in part) be explained by their tendency to separate existing clusters well in the visualization. However, while they preserve local information well, they are known to distort global distances. 

We propose a modular approach for visualization of clustered data that consists of first clustering the data, then embedding each cluster, and finally aligning the clusters to obtain a global embedding. 
By embedding each cluster individually, we aim to produce an embedding that suffers from less distortion than if we were to embed all of the data together. 
Additionally, the method includes a tuning parameter to allow for separation between clusters. The method is competitive with existing methods, while offering much more flexibility and transparency. 

Below is an example of the C+E approach on the MNIST dataset (left) and on the developmental human brain organoid data from [Kanton et al. (2019)](https://www.nature.com/articles/s41586-019-1654-9) (right).

<p align="center"><img width="800" alt="Example of the cluster + embed approach on real data" src="/figures/ce_examples.png">


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

These can be installed by running
```
source("install_packages.R")
```


The code can be cloned by running the following command
```
git clone https://github.com/lizzycoda/cluster_embed.git
```



## Basic Usage


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

Any clustering method can then be used to cluster the data. In this example, we use DBSCAN. The embedding is obtained using the `cluster_align` function. In this example each cluster is embedded via PCA and the clusters are aligned to minimize Kruskal's formulation of the stress. The parameter `alpha` can also (optionally) be set to increase separation between clusters in the embedding. 

```
clusters <-  dbscan(X, eps = .1, minPts = 10)$cluster
out<- cluster_align(euclidean_distances, clusters,
    embedding_method = "PCA",
    alignment_method = "stress",
    use_grad = T,
    alpha = 1.42
  )
ce_embedding <- out$embeddings
```

For a more detailed description of the main functionality see this [tutorial](https://github.com/lizzycoda/cluster_embed/blob/main/tutorial/tutorial.md). 

## Organization

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
  
  
## Data

All synthetic datasets used in the paper can be loaded using `load_data(dataset = c("shapes_3d", "shapes_cylinder", "gmm_10d"), n = 1000)`. 

The human brain organoid data is from [Kanton et al. (2019)](https://www.nature.com/articles/s41586-019-1654-9). We downloaded the [metadata](https://www.ebi.ac.uk/biostudies/files/E-MTAB-7552/metadata_human_cells.tsv) and [counts data]("https://www.ebi.ac.uk/biostudies/files/E-MTAB-7552/human_cell_counts_consensus.mtx"), and used the preprocessing of [Damrich et al.](https://www.biorxiv.org/content/10.1101/2024.04.26.590867v1.abstract), using [this code](https://github.com/berenslab/ne_spectrum_scRNAseq/tree/main).

The mouse cortex data is from [Tasic et al. (2018)](https://www.nature.com/articles/s41586-018-0654-5) downloaded from [here](http://celltypes.brain-map.org/api/v2/well_known_file_download/694413985). We then used the preprocessing in [Kobak & Berens, 2019](https://www.nature.com/articles/s41467-019-13056-x) with code available [here](https://github.com/berenslab/rna-seq-tsne/blob/master/tasic-et-al.ipynb). We then filtered for inhibitory neurons.

The MNIST dataset was downloaded using the [torchvision API](https://docs.pytorch.org/vision/stable/generated/torchvision.datasets.MNIST.html) and the Fashion MNIST dataset was downloaded from [Kaggle](https://www.kaggle.com/datasets/zalando-research/fashionmnist).

For all non-synthetic datasets we randomly subsampled 5000 examples.


## Citation 

```
@article{cluster_embed25,
      title={Cluster and then Embed: A Modular Approach for Visualization}, 
      author={Elizabeth Coda and Ery Arias-Castro and Gal Mishne},
      journal={arXiv preprint arXiv: 2509.03373}
      year={2025}
}
```