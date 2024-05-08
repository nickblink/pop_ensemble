library(tmaptools)
library(igraph)
library(dplyr)
library(rstan)
library(ggplot2)

rstan_options(auto_write = TRUE)

# set working directory
if(file.exists('C:/Users/Admin-Dell')){
  root_dir = 'C:/Users/Admin-Dell'
}else{
  root_dir = 'C:/Users/nickl'
}
setwd(sprintf('%s/Documents/github_projects/pop_ensemble/', root_dir))

# load extra functions
source('code/extra_functions_CAR.R')

## pull in the data
D2010 = read.csv('data/merged_wp_census_data2_081122.csv')
county_adj = read.csv('data/countyadj2.csv', row.names = 1)

## subset data by state
NY_lst <- subset_data_by_state(D2010, county_adj, 'New York', 'NY')

A <- NY_lst$adjacency



g <- make_ring(30, directed=TRUE)

al <- as_adj_list(g, mode="out")

mc <- map_coloring(x = NY_lst$adjacency)

table(mc)


test <- list('A' = 2, '2' = c(1,3,4), '3' = c(2,4), '4' = c(2,3))
map_coloring(x = test)

Adj_list <- list()
for(i in 1:nrow(A)){
  Adj_list[[rownames(A)[i]]] <- which(A[i,] == 1)
}

tt <- map_coloring(Adj_list, ncols = 62)

### Make an adjacency list out of an adjacency matrix.
# adj_mat: An adjacency matrix of 0's and 1's
make_adjacency_list <- function(adj_mat){
  adj_list <- list()
  
  # cycle through each row
  for(i in 1:nrow(adj_mat)){
    adj_list[[rownames(adj_mat)[i]]] <- which(adj_mat[i,] == 1)
  }
  
  return(adj_list)
}

### Make K non-neighboring folds of data from an adjacency matrix.
# adj_mat: An adjacency matrix of 0's and 1's
# K: Number of folds. Should be at least 5.
make_data_folds <- function(adj_mat, K){
  
  if(K < 5){
    stop('too few folds.')
  }
  
  # make the adjacency list
  adj_list <- make_adjacency_list(adj_mat)
  
  # make the map
  groups <- tmaptools::map_coloring(adj_list, ncols = K)
  
  print(table(groups))
  
  return(groups)
}

