plot_embeddings2 <- function(obj, membership_vec, data_1 = T, data_2 = F,
                             nn, frnn_approx = 0, radius_quantile = 0.9,
                             verbose = T){
  embedding <- .prepare_embeddings(obj, data_1 = data_1, data_2 = data_2, 
                                   add_noise = add_noise)
  n <- nrow(embedding[1])
  
  # extract frnn's 
  list_g <- construct_frnn(c_embedding = embedding[1], d_embedding = embedding[2], 
                           e_embedding = embedding[3], 
                           cell_subidx = 1:n, nn = nn, 
                           frnn_approx = frnn_approx, radius_quantile = radius_quantile,
                           bool_matrix = F, verbose = verbose)
  
  # symmetrize
  
  # convert into kernels
  
  # pass into python
}