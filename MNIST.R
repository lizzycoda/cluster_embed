library(dbplyr)
library(ggplot2)
library(vegan)

#### Data loading and preprocessing #### 

#load first 5k training sample from MNIST dataset: https://www.kaggle.com/datasets/hojjatk/mnist-dataset


#calculate geodesic distances 


#### STEP 1: Cluster the data #### 





#### STEP 2-3: Embed each cluster and align the clusters #### 

out = some_function(distances, clusters, embedding_method = 'MDS' )
out = some_function(distances, clusters, embedding_method = 'MDS' , alpha = 3)
out = some_function(distances, clusters, embedding_method = 'kl_div' , k = 10 , scale = 6, alpha = 3)
out = some_function(distances, clusters, embedding_method = 'local_MDS', alpha = 3)

#generate figures
p1 <- ggplot()
#save plots



#### METRICS ####

N = 10
results = matrix(NA, nrow = N, ncol = 6)
# run each method N times and evaluate

for(i in 1:N){
  # fit each method and evaluate each time 
 
  # just one line of code to do all of the evaluation in  
}

# average the results for reporting 



#### BASELINES FOR COMPARISON ####

#isomap


#tsne 


# evaluation of the methods for comparison 



#### VISUALIZE DIGITS ####

# look at the digits from comparisons 
