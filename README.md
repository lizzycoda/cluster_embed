# Cluster + Embed 

This repository contains the code to reproduce the experiments and figures in the paper [Cluster and then Embed: A Modular Approach for Visualization](todo). 

<p align="center"><img width="800" alt="Schematic of the cluster + embed approach" src="/figures/ce_schematic.png">

Dimensionality reduction methods such as t-SNE and UMAP are popular methods for visualizing data with a potential (latent) clustered structure. The popularity of these methods can (in part) be explained by their tendency to separate existing clusters well in the visualization. However, while they preserve local information well, they are known to distort global distances. 

We propose a modular approach for visualization of clustered data that consists of first clustering the data, then embedding each cluster, and finally aligning the clusters to obtain a global embedding. 
By embedding each cluster individually, we aim to produce an embedding that suffers from less distortion than if we were to embed all of the data together. 
Additionally, the method includes a tuning parameter to allow for separation between clusters. The method is competitive with existing methods, while offering much more flexibility and transparency. 

All experiments in the paper are in the examples folder, with a separate file for each dataset.

Each file contains code to obtain the C+E embedding and to evaluate the embedding on a variety of different metrics. 

Below is an example of the C+E approach on the MNIST dataset (left) and on the developmental human brain organoid data from [Kanton et al. (2019)](https://www.nature.com/articles/s41586-019-1654-9).

<p align="center"><img width="800" alt="Example of the cluster + embed approach on real data" src="/figures/ce_examples.png">

## Files 

The repository contains the following files:

- `scripts/`
  - `cluster_embed.R` contains the main code for the method.
  - `generate_data.R` generates the data in our synthetic examples.
  - `utils.R` contains useful embedding and clustering functions.
  - `eval.R` contains the code to evaluate embeddings on the metrics used in the paper, as well as produce plots.
  - `loe_embedding.R` embeds data using a gradient descent approach to minimize the Local Ordinal Embedding objective in [Terada and Von Luxburg (2014)](https://proceedings.mlr.press/v32/terada14.pdf).
  - `stress_functions.R` contains helper functions for the main method.
 


- `examples/` contains the examples in the paper, listed below in the order they appear in the paper.
  - `gaussian_mixture.R`
  - `shapes1.R`
  - `shapes2.R`
  - `mnist.R`
  - `human-409b2.R`
  - `fmnist.R`
  - `neurons.R`

