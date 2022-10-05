#' Fine tune the common tilts, one for each latent dimension
#' 
#' Tune each latent dimension (representing a pair of canonical score vectors)
#' for an appropriate tilt of the common component
#'
#' @param input_obj \code{multiSVD} class, after creation via \code{tiltedCCA()} 
#' @param max_iter maximum number of epochs (i.e., cycling through all latent dimensions)
#' @param fix_tilt_perc scalar between 0 and 1, or \code{NA}, for setting all the tilts in each latent dimension
#' to the scalar if not \code{NA}
#' @param temp_path filepath for saving temporary progress to
#' @param tol small positive number
#' @param verbose non-negative integer             
#'
#' @return updated \code{multiSVD} object
#' @export
fine_tuning <- function(input_obj, 
                        max_iter = 5,
                        fix_tilt_perc = NA,
                        temp_path = NULL,
                        tol = 1e-5,
                        verbose = 0){
  
  if(verbose >= 1) print(paste0(Sys.time(),": Tilted-CCA: Gathering relevant objects"))
  
  n <- nrow(.get_SVD(input_obj)$u)
  metacell_clustering_list <- .get_metacell(input_obj,
                                       resolution = "cell", 
                                       type = "list", 
                                       what = "metacell_clustering")
  if(!all(is.null(metacell_clustering_list))){
    averaging_mat <- .generate_averaging_matrix(metacell_clustering_list = metacell_clustering_list,
                                                n = n)
  } else {
    averaging_mat <- NULL
  }
  
  input_obj <- .set_defaultAssay(input_obj, assay = 1)
  svd_1 <- .get_SVD(input_obj)
  score_1 <- .get_tCCAobj(input_obj, apply_postDimred = F, what = "score")
  distinct_score_1 <- .get_tCCAobj(input_obj, apply_postDimred = F, what = "distinct_score")
  
  input_obj <- .set_defaultAssay(input_obj, assay = 2)
  svd_2 <- .get_SVD(input_obj)
  score_2 <- .get_tCCAobj(input_obj, apply_postDimred = F, what = "score")
  distinct_score_2 <- .get_tCCAobj(input_obj, apply_postDimred = F, what = "distinct_score")
  
  target_dimred <- .get_Laplacian(input_obj, bool_common = T)
  common_score <- .get_tCCAobj(input_obj, apply_postDimred = F, what = "common_score")
  rank_c <- ncol(common_score)
  n <- nrow(svd_1$u)
  param <- .get_param(input_obj)
  snn_bool_cosine <- param$snn_bool_cosine
  snn_bool_intersect <- param$snn_bool_intersect
  snn_k <- param$snn_latent_k
  snn_min_deg <- param$snn_min_deg
  snn_num_neigh <- param$snn_num_neigh
  
  stopifnot(all(is.na(fix_tilt_perc)) || length(fix_tilt_perc) == rank_c)
  stopifnot(all(sapply(1:rank_c, function(k){
    val <- score_1[,k] %*% score_2[,k]; val >= 0 
  }))) # ensures score matrices contain pair of acute vectors
  
  basis_list <- lapply(1:rank_c, function(k){
    .representation_2d(score_1[,k], score_2[,k])
  })
  
  circle_list <- lapply(1:rank_c, function(k){
    vec1 <- basis_list[[k]]$rep1
    vec2 <- basis_list[[k]]$rep2
    .construct_circle(vec1, vec2)
  })
  
  if(verbose >= 1) print(paste0(Sys.time(),": Tilted-CCA: Starting optimization"))
  iter <- 1
  percentage_grid <- seq(0, 1, length.out = param$tcca_discretization_gridsize)
  common_score_prev <- common_score
  if(all(is.na(fix_tilt_perc))) {
    tilt_perc <- rep(input_obj$tcca_obj$tilt_perc[1], rank_c)
  } else {
    tilt_perc <- fix_tilt_perc
  }
  
  while(iter <= max_iter){
    for(k in 1:rank_c){
      if(all(is.na(fix_tilt_perc))){
        tmp <- percentage_grid
      } else {
        tmp <- tilt_perc[k]
      }
      
      res <- .fine_tuning_dimension(averaging_mat = averaging_mat,
                                    basis_list = basis_list,
                                    circle_list = circle_list,
                                    common_score = common_score,
                                    latent_dim = k,
                                    percentage_grid = tmp,
                                    score_1 = score_1,
                                    score_2 = score_2,
                                    snn_bool_cosine = snn_bool_cosine,
                                    snn_bool_intersect = snn_bool_intersect,
                                    snn_k = snn_k,
                                    snn_min_deg = snn_min_deg,
                                    snn_num_neigh = snn_num_neigh,
                                    svd_1 = svd_1, 
                                    svd_2 = svd_2,
                                    target_dimred = target_dimred,
                                    verbose = verbose)
      if(verbose >= 1) {
        print(paste0("On iteration ", iter, " for latent dimension ", k))
        print(res$df)
        print("=======")
      }
      tilt_perc[k] <- res$percentage
      common_score <- res$common_score
    }
    
    if(iter > 2 && sum(abs(common_score - common_score_prev)) <= tol) break()
    if(!is.null(temp_path) && is.character(temp_path)){
      tmp <- .compute_distinct_score(score_1, score_2, common_score)
      distinct_score_1 <- tmp$distinct_score_1; distinct_score_2 <- tmp$distinct_score_2
      res <- list(common_score = common_score, 
                  distinct_score_1 = distinct_score_1,
                  distinct_score_2 = distinct_score_2,
                  score_1 = score_1, score_2 = score_2, 
                  svd_1 = svd_1, svd_2 = svd_2, 
                  tilt_perc = tilt_perc)
      save(res, file = temp_path)
    }
    
    iter <- iter + 1
    common_score_prev <- common_score
  }
  
  tmp <- .compute_distinct_score(score_1, score_2, common_score)
  distinct_score_1 <- tmp$distinct_score_1; distinct_score_2 <- tmp$distinct_score_2
  
  tcca_obj <- .create_tcca_obj(common_basis = res$common_basis,
                               common_score = res$common_score, 
                               distinct_score_1 = distinct_score_1,
                               distinct_score_2 = distinct_score_2,
                               df_percentage = NA,
                               tilt_perc = tilt_perc)
  input_obj$tcca_obj <- tcca_obj
  input_obj$param$ft_max_tier <- max_iter
  
  input_obj
}

##############################

.fine_tuning_dimension <- function(averaging_mat,
                                   basis_list,
                                   circle_list,
                                   common_score,
                                   latent_dim,
                                   percentage_grid,
                                   score_1,
                                   score_2,
                                   snn_bool_cosine,
                                   snn_bool_intersect,
                                   snn_k,
                                   snn_min_deg,
                                   snn_num_neigh,
                                   svd_1, 
                                   svd_2,
                                   target_dimred,
                                   verbose = 0){
  r <- length(basis_list)
  
  value_list <- lapply(percentage_grid, function(percentage){
    if(verbose >= 2) print(paste0("Working on percentage ", percentage))
    radian_val <- .compute_radian(circle = circle_list[[latent_dim]],
                                  enforce_boundary = T,
                                  percentage_val = percentage, 
                                  vec1 = basis_list[[latent_dim]]$rep1,
                                  vec2 = basis_list[[latent_dim]]$rep2)
    
    common_representation_new <- .position_from_circle(circle_list[[latent_dim]], radian_val)
    common_score[,latent_dim] <- basis_list[[latent_dim]]$basis_mat %*% common_representation_new
    
    common_mat <- .convert_common_score_to_mat(common_score,
                                               score_1,
                                               score_2,
                                               svd_1,
                                               svd_2)
    if(!all(is.null(averaging_mat))){
      avg_common_mat <- as.matrix(averaging_mat %*% common_mat)
    } else {
      avg_common_mat <- common_mat
    }
    
    snn_mat <- .form_snn_mat(bool_cosine = snn_bool_cosine,
                             bool_intersect = snn_bool_intersect,
                             mat = avg_common_mat, 
                             min_deg = snn_min_deg,
                             num_neigh = snn_num_neigh,
                             verbose = verbose)
    common_basis <- .compute_laplacian_basis(snn_mat, 
                                             latent_k = snn_k,
                                             verbose = verbose)
    
    value <- .grassmann_distance(orthonormal_1 = common_basis, 
                                 orthonormal_2 = target_dimred)
    
    list(value = value, 
         common_vec = common_score[,latent_dim],
         common_basis = common_basis)
  })
  names(value_list) <- percentage_grid
  value_vec <- sapply(value_list, function(x){x$value})
  
  idx_min <- .select_minimum(minimum = T,
                             x_val = percentage_grid,
                             y_val = value_vec)
  df <- data.frame(percentage = percentage_grid,
                   ratio_val = value_vec)
  prev_vec <- common_score[,latent_dim]
  common_score[,latent_dim] <- value_list[[idx_min]]$common_vec
  common_basis <- value_list[[idx_min]]$common_basis
  
  list(common_basis = common_basis,
       common_score = common_score,
       df = df, 
       percentage = percentage_grid[idx_min])
}