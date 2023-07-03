#' Report heat map
#'
#' @param res.sel.var A list.
#' @param output_dir A string.
#'
#' @importFrom ComplexHeatmap HeatmapAnnotation
#' @importFrom grid gpar
#' @importFrom circlize colorRamp2
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom ComplexHeatmap draw
#'
cellCellRegHeatmap <- function(res.sel.var = res.sel.var,
                               output_dir = output_dir){

  start_time <- Sys.time()
  print(paste0("=== cellCellHeatmap started... ", Sys.time(), " ==="))

  cell_type_list <- res.sel.var$cell_type_list
  cell_type_num <- length(cell_type_list)
  sig_coef_mat_2_list <- vector("list", length = cell_type_num)
  
  for(cell_type_ind in 1:cell_type_num){

    # Judge model was built or not.
    if (!is.null(res.sel.var$all_zero_flag_list[[cell_type_ind]])){
      if (!res.sel.var$all_zero_flag_list[[cell_type_ind]]){

        sig_coef_mat <- res.sel.var$sig_coef_mat_list[[cell_type_ind]]

        gene_cluster_vec <- res.sel.var$gene_cluster_vec_list[[cell_type_ind]]

        gene_cluster_vec_2 <- gene_cluster_vec[order(gene_cluster_vec)]

        if (nrow(sig_coef_mat) == 1){
          sig_coef_mat_2 <- t(sig_coef_mat[, names(gene_cluster_vec_2)])
          rownames(sig_coef_mat_2) <- rownames(sig_coef_mat)
        } else {
          sig_coef_mat_2 <- sig_coef_mat[, names(gene_cluster_vec_2)]
        }

        max_sig <- max(sig_coef_mat_2)
        min_sig <- min(sig_coef_mat_2)

        row_num_heat <- nrow(sig_coef_mat_2)

        if (row_num_heat > 1){
          row_ft_sz <- 100 / row_num_heat
        } else {
          row_ft_sz <- 100 / 3
        }

        set.seed("124")
        column_ha = ComplexHeatmap::HeatmapAnnotation(Cluster = as.character(gene_cluster_vec_2),
                                                      annotation_name_gp = grid::gpar(fontsize = 20))
        if (min_sig < 0 & max_sig > 0){
          col_fun = circlize::colorRamp2(c(min_sig, 0, max_sig), c("blue", "gray", "red"))
        } else if (min_sig >= 0 & max_sig >= 0){
          col_fun = circlize::colorRamp2(c(min_sig, max_sig), c("gray", "red"))
        } else if (min_sig <= 0 & max_sig <= 0){
          col_fun = circlize::colorRamp2(c(min_sig, max_sig), c("blue", "gray"))
        }

        coef_p <- ComplexHeatmap::Heatmap(sig_coef_mat_2,
                                          name = "Coefficient",
                                          cluster_columns = FALSE,
                                          show_column_names=FALSE,
                                          column_title = paste0("HVGs in ",
                                                                cell_type_list[[cell_type_ind]]),
                                          column_title_side = "top",
                                          column_title_gp = grid::gpar(fontsize = 20),
                                          top_annotation = column_ha,
                                          row_names_gp = grid::gpar(fontsize = row_ft_sz),
                                          col = col_fun,
                                          cluster_rows = FALSE)

        png(paste0(output_dir, "/cell_type_", cell_type_ind, "_heatmap.png"))
        ComplexHeatmap::draw(coef_p)
        dev.off()

        pdf(paste0(output_dir, "/cell_type_", cell_type_ind, "_heatmap.pdf"))
        ComplexHeatmap::draw(coef_p)
        dev.off()

        sig_coef_mat_2_list[[cell_type_ind]] <- sig_coef_mat_2
        
      } else {

        sig_coef_mat_2_list[[cell_type_ind]] <- NULL
        
        sink(paste0(output_dir, "/cell_type_", cell_type_ind, "_heatmap.txt"))
        print(paste0("No significant combinations of genes and neighboring cell types in ", res.sel.var$cell_type_list[[cell_type_ind]]))
        sink()

      }

    } else {

      sig_coef_mat_2_list[[cell_type_ind]] <- NULL
      
      sink(paste0(output_dir, "/cell_type_", cell_type_ind, "_heatmap.txt"))
      print(paste0("NULL returned in ", res.sel.var$cell_type_list[[cell_type_ind]]))
      sink()

    }
  }

  finish_time <- Sys.time()
  print(paste0("=== cellCellHeatmap finished... ", Sys.time(), " ==="))
  print(paste0("=========="))

  return(list(sig_coef_mat_2_list = sig_coef_mat_2_list))

}
