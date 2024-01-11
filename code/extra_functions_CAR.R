
subset_data_by_state <- function(data, adjacency, state, abbrev = NULL){
  # get the indices corresponding to the state
  ind1 <- grep(state, data$NAME)
  ind2 <- grep(state, rownames(adjacency))
  
  # check that the indices match between data and adjacency matrix
  if(!identical(ind1, ind2)){
    stop('indices of data names and adjacency names do not match')
  }
  
  # check abbreviated adjacency column names
  if(!is.null(abbrev)){
    ind3 <- grep(abbrev, colnames(adjacency))
    if(!identical(ind1, ind3)){
      stop('indices of data names and adjacency columnss do not match')
    }
  }
  
  # subset the data
  data_subset <- data[ind1, ]
  adjacency_subset <- adjacency[ind1, ind1]
  
  # return it!
  return(list(data = data_subset, adjacency = adjacency_subset))
}

### simulate the numbers in the data according to a normal distribution
# data: the input data
# models: the models to simulate data for
# means: the means of the normals for each model
# variances: the variances of the normals for each model
simulate_numbers <- function(data, models, means, variances){
  # check that the lengths of the inputs match
  if(length(models) != length(means) | length(models) != length(variances)){
    stop('length of models, means, and variances, not matching up')
  }
  
  # sample data for each model
  for(i in 1:length(models)){
    m = models[i]
    data[,m] <- rnorm(nrow(data), means[i], sd = sqrt(variances[i]))
  }
  
  return(data)
}
  

simulate_data <- function(data, adjacency, models, sim_numbers = F, ...){
  
  if(sim_numbers){
    # hmm I don't like this. I want a simple way for the original data to be 
    # transformed and then... what? There are many options to change the data.
  }
}
  