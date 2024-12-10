#' Filter coefficient
#'
#' @param res.pls A list.
#' @param exp_mat_norm A data frame.
#' @param fet_mat_norm A data frame.
#' @param p_value_thresh A number.
#'
#' @importFrom pls plsr
#'
sigCoef <- function(res.pls,
                    exp_mat_norm,
                    fet_mat_norm,
                    p_value_thresh = 0.05){

  ## Prepare variable
  res_pls <- res.pls$res_pls
  opt_comp_num <- res.pls$opt_comp_num

  x_score_mat <- res_pls$scores
  y_score_mat <- res_pls$Yscores

  ## Calculate sig_coef_mat

  # Declare sig_coef_mat
  sig_coef_mat <- matrix(0, nrow = nrow(res_pls$coefficient[,,1]),
                         ncol = ncol(res_pls$coefficients[,,1]))
  rownames(sig_coef_mat) <- rownames(res_pls$coefficient[,,1])
  colnames(sig_coef_mat) <- colnames(res_pls$coefficient[,,1])
  sig_coef_mat_raw <- sig_coef_mat

  # Add each component
  for (calc_comp_ind in 1:opt_comp_num){

    if (calc_comp_ind > 1){
      sig_coef_mat_comp_raw <- res_pls$coefficients[,,calc_comp_ind] - res_pls$coefficients[,,calc_comp_ind - 1]
    } else {
      sig_coef_mat_comp_raw <- res_pls$coefficients[,,calc_comp_ind]
    }

    # Extract non-significant features
    del_fet_ind_vec <- returnDelFetInd(fet_mat_norm = fet_mat_norm,
                                       x_score_mat = x_score_mat,
                                       calc_comp_num = calc_comp_ind,
                                       p_value_thresh = p_value_thresh)

    # Extract non-significant genes
    del_gene_ind_vec <- returnDelGeneInd(exp_mat_norm = exp_mat_norm,
                                         y_score_mat = y_score_mat,
                                         calc_comp_num = calc_comp_ind,
                                         p_value_thresh = p_value_thresh)

    # Replace non-significant value with 0
    sig_coef_mat_comp <- sig_coef_mat_comp_raw
    sig_coef_mat_comp[del_fet_ind_vec, ] <- 0
    sig_coef_mat_comp[, del_gene_ind_vec] <- 0

    sig_coef_mat_raw <- sig_coef_mat_raw + sig_coef_mat_comp_raw
    sig_coef_mat <- sig_coef_mat + sig_coef_mat_comp

  }

  return(list(sig_coef_mat = sig_coef_mat,
              sig_coef_mat_raw = sig_coef_mat_raw))

}

#' Get index of feature for deleting
#'
#' @param fet_mat_norm A data frame.
#' @param x_score_mat A data frame.
#' @param calc_comp_num A number.
#' @param p_value_thresh A number.
#'
#' @importFrom stats cor.test
#'
returnDelFetInd <- function(fet_mat_norm,
                            x_score_mat,
                            calc_comp_num,
                            p_value_thresh){

  # Remove self
  self_col_ind <- grep("self", colnames(fet_mat_norm))
  fet_mat_norm2 <- fet_mat_norm[, -self_col_ind]

  fet_p_vec <- c()
  del_fet_ind_vec <- c()

  fet_num <- ncol(fet_mat_norm2)

  for (fet_ind in 1:fet_num){
    fet_p_value <- stats::cor.test(fet_mat_norm2[, fet_ind], x_score_mat[, calc_comp_num])$p.value
    fet_p_vec <- append(fet_p_vec, fet_p_value)
  }

  fet_q_vec <- stats::p.adjust(fet_p_vec, method = "BH")
  del_fet_bin <- rep(1, fet_num) * (fet_q_vec >= p_value_thresh)

  for (fet_ind in 1:fet_num){
    if (del_fet_bin[fet_ind] == 1)
    del_fet_ind_vec <- c(del_fet_ind_vec, fet_ind)
  }

  return(del_fet_ind_vec)

}

#' Get index of gene for deleting
#'
#' @param exp_mat_norm A data frame.
#' @param y_score_mat A data frame.
#' @param calc_comp_num A number.
#' @param p_value_thresh A number.
#'
#' @importFrom stats cor.test
#'
returnDelGeneInd <- function(exp_mat_norm,
                             y_score_mat,
                             calc_comp_num,
                             p_value_thresh){
  gene_p_vec <- c()
  del_gene_ind_vec <- c()

  gene_num <- ncol(exp_mat_norm)

  for (gene_ind in 1:gene_num){
    gene_p_value <- stats::cor.test(exp_mat_norm[, gene_ind], y_score_mat[, calc_comp_num])$p.value
    gene_p_vec <- append(gene_p_vec, gene_p_value)
  }

  gene_q_vec <- stats::p.adjust(gene_p_vec, method = "BH")
  del_gene_bin <- rep(1, gene_num) * (gene_q_vec >= p_value_thresh)

  for (gene_ind in 1:gene_num){
    if (del_gene_bin[gene_ind] == 1)
      del_gene_ind_vec <- c(del_gene_ind_vec, gene_ind)
  }

  return(del_gene_ind_vec)

}

#' Filter coefficient
#'
#' @param res.estimate A list.
#' @param res.sep.mat A list.
#' @param HVG_extract_num A number.
#' @param dev_opt A string.
#'
#' @importFrom purrr map_dbl
#' @importFrom cluster pam
#' @importFrom stats ecdf
#' @importFrom stats p.adjust
#' @importFrom stats kmeans
#'
cellCellRegSelVar <- function(res.estimate = res.estimate,
                              res.sep.mat = res.sep.mat,
                              HVG_extract_num = HVG_extract_num,
                              dev_opt = dev_opt){
  
  
  start_time <- Sys.time()
  print(paste0("=== cellCellRegSelVar started... ", Sys.time(), " ==="))

  set.seed(123)

  cell_type_num <- length(res.estimate$cell_type_list)

  res_xmeans_list <- vector("list", length = cell_type_num)

  sig_coef_mat_non_zero_list <- vector("list", length = cell_type_num)
  sig_coef_mat_with_zero_list <- vector("list", length = cell_type_num)

  sig_coef_mat_list <- vector("list", length = cell_type_num)
  ave_cluster_coef_mat_list <- vector("list", length = cell_type_num)

  sig_coef_mat_bin_2_list <- vector("list", length = cell_type_num)
  gene_cluster_vec_list <- vector("list", length = cell_type_num)

  all_zero_flag_list <- vector("list", length = cell_type_num)

  for (cell_type_ind in 1:cell_type_num){

    # Judge model was built or not.
    if (res.estimate$CCPLS_result[[cell_type_ind]][[1]] != "No model was built."){

      # Get significant coefficient
      res.sig.coef <- sigCoef(res.pls = res.estimate$CCPLS_result[[cell_type_ind]][[2]],
                             exp_mat_norm = scale(res.sep.mat$exp_mat_sep_list[[cell_type_ind]]),
                             fet_mat_norm = res.sep.mat$fet_mat_sep_list[[cell_type_ind]])

      sig_coef_mat <- res.sig.coef$sig_coef_mat
      sig_coef_mat_raw <- res.sig.coef$sig_coef_mat_raw

      non_sig_gene_col <- colnames(sig_coef_mat)[apply(sig_coef_mat, 2, sum) == 0]
      sig_gene_col <- colnames(sig_coef_mat)[apply(sig_coef_mat, 2, sum) != 0]
      non_sig_fet_row <- rownames(sig_coef_mat)[apply(sig_coef_mat, 1, sum) == 0]
      sig_fet_row <- rownames(sig_coef_mat)[apply(sig_coef_mat, 1, sum) != 0]

      sig_coef_mat_bin <- matrix(0, nrow(sig_coef_mat), ncol(sig_coef_mat))
      rownames(sig_coef_mat_bin) <- rownames(sig_coef_mat)
      colnames(sig_coef_mat_bin) <- colnames(sig_coef_mat)

      for (fet_id in 1:nrow(sig_coef_mat)){

        if (length(non_sig_gene_col) != 0){
          null_dist_fet <- stats::ecdf(sig_coef_mat_raw[fet_id, non_sig_gene_col])
          fet_p_vec <- 1 - null_dist_fet(sig_coef_mat_raw[fet_id,])
          names(fet_p_vec) <- colnames(sig_coef_mat_raw)

          fet_q_vec <- stats::p.adjust(fet_p_vec, method = "BH")

          sig_coef_vec_fet <- sig_coef_mat[fet_id,]
          sig_plus_flag <- fet_q_vec < 0.05 & sig_coef_vec_fet > 0
          sig_minus_flag <- fet_q_vec < 0.05 & sig_coef_vec_fet < 0
          sig_coef_mat_bin[fet_id, sig_plus_flag] <- 1
          sig_coef_mat_bin[fet_id, sig_minus_flag] <- -1

        }

      }

      sig_coef_mat_bin_2 <- sig_coef_mat_bin
      sig_coef_mat_bin_2[, non_sig_gene_col] <- 0
      sig_coef_mat_bin_2[non_sig_fet_row, ] <- 0

      sig_coef_mat_bin_2_non_zero <- sig_coef_mat_bin_2[, colSums(sig_coef_mat_bin_2) != 0]

      sig_coef_mat[sig_coef_mat_bin_2 == 0] <- 0
      sig_coef_mat_with_zero_list[[cell_type_ind]] <- sig_coef_mat
      sig_coef_mat_non_zero <- sig_coef_mat[, colSums(sig_coef_mat_bin_2) != 0]
      # zero_gene_num <- ncol(sig_coef_mat) - ncol(sig_coef_mat_non_zero)

      sig_coef_mat_non_zero_list[[cell_type_ind]] <- sig_coef_mat_non_zero

      all_zero_flag <- sum(sig_coef_mat_bin_2) == 0
      all_zero_flag_list[[cell_type_ind]] <- all_zero_flag

      if (!all_zero_flag){

        if (dev_opt == "kmeans"){

          set.seed(123)

          if (!is.null(nrow(sig_coef_mat_non_zero[rowSums(sig_coef_mat_non_zero) != 0,]))){
            mat <- scale(t(sig_coef_mat_non_zero[rowSums(sig_coef_mat_non_zero) != 0,]))
          } else {
            mat <- t(sig_coef_mat_non_zero[rowSums(sig_coef_mat_non_zero) != 0,])
            colnames(mat) <- rownames(sig_coef_mat_non_zero)[rowSums(sig_coef_mat_non_zero) != 0]
            mat <- scale(mat)
          }
                 
          k.max <- 15
          if ((nrow(mat) - 1) < 15){
            k.max <- nrow(mat) - 1
          }

          if (k.max == 1){
            k.max <- 2
          }
          
          sil_width <- purrr::map_dbl(2:k.max,  function(k){
            model <- cluster::pam(x = mat, k = k)
            model$silinfo$avg.width
          })

          cluster_num <- which.max(sil_width) + 1
          res.kmeans <- stats::kmeans(mat, cluster_num)
          gene_cluster_vec <- res.kmeans$cluster
          gene_cluster_vec_list[[cell_type_ind]] <- gene_cluster_vec

        } else if (dev_opt == "binary"){

          gene_bin_vec_list <- vector("list", ncol(sig_coef_mat_bin_2_non_zero))
          for (gene_id in 1:ncol(sig_coef_mat_bin_2_non_zero)){
            gene_bin_vec_list[[gene_id]] <- as.character(sig_coef_mat_bin_2_non_zero[,gene_id])
          }

          cluster_list <- unique(gene_bin_vec_list)
          cluster_num <- length(cluster_list)

          gene_cluster_vec <- c()
          for (cluster_id in 1:cluster_num){
            for (gene_id in 1:length(gene_bin_vec_list)){
              judge_flag <- nrow(sig_coef_mat_bin_2) == sum(gene_bin_vec_list[[gene_id]] == cluster_list[[cluster_id]])
              if (judge_flag){
                gene_cluster_vec <- append(gene_cluster_vec, cluster_id)
              }
            }
          }

          names(gene_cluster_vec) <- colnames(sig_coef_mat_bin_2_non_zero)
          gene_cluster_vec_list[[cell_type_ind]] <- gene_cluster_vec

        }

        cluster_name <- c()
        ave_cluster_coef_mat_non_zero <- matrix(0, nrow = nrow(sig_coef_mat_non_zero_list[[cell_type_ind]]),
                                                ncol = cluster_num)

        for (cluster_ind in 1:cluster_num){

          cluster_gene <- names(gene_cluster_vec[gene_cluster_vec == cluster_ind])
          cluster_coef_mat_non_zero <- sig_coef_mat_non_zero_list[[cell_type_ind]][, cluster_gene]

          if (!is.null(ncol(cluster_coef_mat_non_zero))){ # judge column is 1 or not
            ave_cluster_coef_mat_non_zero[, cluster_ind] <- apply(cluster_coef_mat_non_zero,
                                                                  1, mean)
          } else {
            ave_cluster_coef_mat_non_zero[, cluster_ind] <- cluster_coef_mat_non_zero
          }

          cluster_name <- append(cluster_name, paste0("SRGs Cluster ", cluster_ind))

        }

        rownames(ave_cluster_coef_mat_non_zero) <- rownames(sig_coef_mat_non_zero)
        colnames(ave_cluster_coef_mat_non_zero) <- cluster_name

        if (ncol(sig_coef_mat_non_zero) <= HVG_extract_num){
          non_DEG_vec <- matrix(0, nrow = nrow(ave_cluster_coef_mat_non_zero), ncol = 1)
          colnames(non_DEG_vec) <- "non-SRGs"
          ave_cluster_coef_mat <- cbind(ave_cluster_coef_mat_non_zero,
                                        non_DEG_vec)
          ave_cluster_coef_mat_list[[cell_type_ind]] <- ave_cluster_coef_mat
        }

        sig_coef_mat_non_zero_ordered <- sig_coef_mat_non_zero_list[[cell_type_ind]][, order(gene_cluster_vec)]

        # Remove zero feature
        if (sum(rowSums(sig_coef_mat_non_zero_ordered) != 0) > 1){
          sig_coef_mat <- sig_coef_mat_non_zero_ordered[rowSums(sig_coef_mat_non_zero_ordered) != 0, ]
        } else {
          sig_coef_mat <- t(sig_coef_mat_non_zero_ordered[rowSums(sig_coef_mat_non_zero_ordered) != 0, ])
          rownames(sig_coef_mat) <- rownames(sig_coef_mat_non_zero_ordered)[rowSums(sig_coef_mat_non_zero_ordered) != 0]
        }

        row_name_rm <- c()
        for (row_ind in 1:nrow(sig_coef_mat)){
          row_name_rm <- append(row_name_rm, strsplit(rownames(sig_coef_mat)[row_ind], split = "neig_")[[1]][2])
        }
        rownames(sig_coef_mat) <- row_name_rm

        sig_coef_mat_list[[cell_type_ind]] <- sig_coef_mat
        sig_coef_mat_bin_2_list[[cell_type_ind]] <- sig_coef_mat_bin_2

      } else { # in case of all zero

        res_xmeans_list [[cell_type_ind]] <- "NULL"
        gene_cluster_vec_list[[cell_type_ind]] <- "NULL"
        ave_cluster_coef_mat_list[[cell_type_ind]] <- "NULL"
        sig_coef_mat_list[[cell_type_ind]] <- "NULL"
        sig_coef_mat_bin_2_list[[cell_type_ind]] <- "NULL"

      }

    }

  }


  finish_time <- Sys.time()
  print(paste0("=== cellCellRegSelVar finished... ", Sys.time(), " ==="))
  print(paste0("=========="))

  return(list(res_xmeans_list = res_xmeans_list,
              ave_cluster_coef_mat_list = ave_cluster_coef_mat_list,
              sig_coef_mat_list = sig_coef_mat_list,
              sig_coef_mat_with_zero_list = sig_coef_mat_with_zero_list,
              sig_coef_mat_bin_2_list = sig_coef_mat_bin_2_list,
              gene_cluster_vec_list = gene_cluster_vec_list,
              all_zero_flag_list = all_zero_flag_list,
              cell_type_list = res.estimate$cell_type_list,
              cell_type_model_SG = res.estimate$cell_type_model_SG))

}
