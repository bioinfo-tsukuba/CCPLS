#' Perform Seurat
#'
#' @param exp_mat A data frame.
#' @param HVG_extract_num A number.
#'
#' @importFrom Seurat CreateSeuratObject
#' @importFrom Seurat NormalizeData
#' @importFrom Seurat FindVariableFeatures
#' @importFrom Seurat VariableFeatures
#' @importFrom Seurat ScaleData
#' @importFrom dplyr %>%
#' @importFrom Matrix Matrix
#'
performSeurat <- function(exp_mat, HVG_extract_num){

  to_dgC <- function(m) {
    if (!inherits(m, "dgCMatrix")) Matrix::Matrix(m, sparse = TRUE) else m
  }

  ### Convert to Seurat Object
  seur_obj <- Seurat::CreateSeuratObject(counts = to_dgC(t(exp_mat)))

  ### Normalize
  seur_obj <- Seurat::NormalizeData(seur_obj)

  ### Identify Highly Variable Genes
  seur_obj <- Seurat::FindVariableFeatures(seur_obj, selection.method = "vst", nfeatures = HVG_extract_num)

  HVG_list <- Seurat::VariableFeatures(seur_obj)

  ### Scaling the data
  seur_obj <- Seurat::ScaleData(seur_obj)

  ### Prepare expression matrix of HVG without using GetAssayData
  exp_mat_seur <- as.matrix(seur_obj@assays$RNA$data[HVG_list,]) %>%
    t()

  return(list(exp_mat_seur = exp_mat_seur,
              HVG_list = HVG_list))
}
#' Perform preprocess
#'
#' @param exp_mat A data frame.
#' @param fet_mat A data frame.
#' @param annot_mat A data frame.
#' @param annot_score_mat A data frame.
#' @param HVG_opt A logical.
#' @param HVG_extract_num A number.
#'
cellCellRegSepMat <- function(exp_mat = exp_mat,
                              fet_mat = fet_mat,
                              annot_mat = annot_mat,
                              annot_score_mat = annot_score_mat,
                              HVG_opt = HVG_opt,
                              HVG_extract_num = HVG_extract_num){


  print(paste0("=== cellCellRegSepMat started... ", Sys.time(), " ==="))
  start_time <- Sys.time()

  exp_mat <- as.matrix(exp_mat)

  cell_row_full_list <- rownames(exp_mat)
  rownames(fet_mat) <- cell_row_full_list

  cell_type_list <- colnames(annot_score_mat)

  row_list <- rownames(exp_mat)
  col_list <- colnames(exp_mat)

  col_list_fet <- colnames(fet_mat)

  exp_mat_sep_list <- vector("list", length = ncol(annot_score_mat))
  exp_mat_sep_each <- matrix(0, nrow = nrow(exp_mat), ncol = ncol(exp_mat))
  fet_mat_sep_list <- vector("list", length = ncol(annot_score_mat))
  exp_mat_orig_sep_list <- vector("list", length = ncol(annot_score_mat))
  HVG_sep_list <- vector("list", length = ncol(annot_score_mat))

  for (cell_type_index in 1:ncol(annot_score_mat)){

    exp_mat_sep_each <- exp_mat[annot_mat[,2] == cell_type_list[cell_type_index],]

    null_flag <- is.null(nrow(exp_mat_sep_each))

    if (null_flag){

      exp_mat_sep_each_del2 <- matrix(0, nrow = 2, ncol = ncol(exp_mat))
      rownames(exp_mat_sep_each_del2) <- paste0("s", seq(1,2))
      colnames(exp_mat_sep_each_del2) <- col_list

      fet_mat_sep_each <- matrix(0, nrow = 2, ncol = ncol(fet_mat))
      colnames(fet_mat_sep_each) <- col_list_fet
      rownames(fet_mat_sep_each) <- paste0("s", seq(1,2))

      exp_mat_sep_each_orig <- exp_mat_sep_each_del2

    } else {

      # Remove cells if all gene expression values are zero
      exp_mat_sep_each_del <- exp_mat_sep_each[apply(exp_mat_sep_each, 1, sum) != 0,]

      # Remove genes if variance is zero
      exp_mat_sep_each_del2 <- exp_mat_sep_each_del[,apply(exp_mat_sep_each_del, 2, var) != 0]

      cell_row_list_del2 <- rownames(exp_mat_sep_each_del2)
      fet_mat_sep_each_unnorm <- fet_mat[cell_row_list_del2,]
      fet_mat_sep_each <- scale(fet_mat_sep_each_unnorm, center = TRUE, scale = TRUE)

    }


    if (HVG_opt == TRUE){

      if (null_flag){

        null_row_names <- rownames(exp_mat_sep_each_del2)
        exp_mat_sep_each_del2 <- matrix(0, nrow = nrow(exp_mat_sep_each_del2), ncol = 1)
        rownames(exp_mat_sep_each_del2) <- null_row_names
        colnames(exp_mat_sep_each_del2) <- "NULL"

        loading_mat <- matrix(0, nrow = nrow(exp_mat_sep_each_del2), ncol = 1)
        rownames(loading_mat) <- null_row_names
        colnames(loading_mat) <- "NULL"

        exp_mat_sep_each_orig <- exp_mat_sep_each_del2

        f.loading_mat <- "NULL"
        q.val.f.loading_mat <- "NULL"

        HVG_list <- "NULL"

      } else {

        exp_mat_sep_each_orig <- exp_mat_sep_each_del

        res.perform.Seurat <- performSeurat(exp_mat = exp_mat_sep_each_del2,
                                            HVG_extract_num = HVG_extract_num)

        exp_mat_sep_each_del2 <- res.perform.Seurat$exp_mat_seur
        HVG_list <- res.perform.Seurat$HVG_list

      }

    }

    fet_mat_sep_list[[cell_type_index]] <- fet_mat_sep_each
    exp_mat_sep_list[[cell_type_index]] <- exp_mat_sep_each_del2
    exp_mat_orig_sep_list[[cell_type_index]] <- exp_mat_sep_each_orig
    HVG_sep_list[[cell_type_index]] <- HVG_list

  }

  finish_time <- Sys.time()
  print(paste0("=== cellCellRegSepMat finished... ", Sys.time(), " ==="))
  print(paste0("=========="))


  return(list(exp_mat_sep_list = exp_mat_sep_list,
              fet_mat_sep_list = fet_mat_sep_list,
              exp_mat_orig_sep_list = exp_mat_orig_sep_list,
              HVG_sep_list = HVG_sep_list,
              cell_type_list = cell_type_list,
              start_time = start_time,
              finish_time = finish_time))

}
