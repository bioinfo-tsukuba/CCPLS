#' Report bipartite graph
#'
#' @param res.sel.var A list.
#' @param output_dir A string.
#'
#' @importFrom dplyr %>%
#' @importFrom igraph graph_from_incidence_matrix
#' @importFrom igraph V
#' @importFrom igraph E
#' @importFrom grDevices png
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf
#'
cellCellRegGraph <- function(res.sel.var = res.sel.var,
                             output_dir = output_dir){

  start_time <- Sys.time()
  print(paste0("=== cellCellGraph started... ", Sys.time(), " ==="))

  ### Bipartite graph

  cell_type_list <- res.sel.var$cell_type_list
  cell_type_num <- length(cell_type_list)

  bi_graph_list <- vector("list", cell_type_num)

  for (cell_type_ind in 1:cell_type_num){

    target_cell_type_name <- cell_type_list[cell_type_ind]

    # Judge model was built or not.
    if (!is.null(res.sel.var$all_zero_flag_list[[cell_type_ind]])){
      if (!res.sel.var$all_zero_flag_list[[cell_type_ind]]){

        ave_cluster_coef_mat <- res.sel.var$ave_cluster_coef_mat_list[[cell_type_ind]]

        plot_coef_mat <- ave_cluster_coef_mat[, colnames(ave_cluster_coef_mat) != "non-SRGs"] * 10 # FIXME：定数
        g_row_num <- nrow(plot_coef_mat)

        right_name_rm <- c()

        for (col_ind in 1:ncol(plot_coef_mat)){
          right_name_rm <- append(right_name_rm, strsplit(colnames(plot_coef_mat),  " ")[[col_ind]][[3]] )
        }
        plot_coef_mat_2 <- plot_coef_mat
        colnames(plot_coef_mat_2) <- right_name_rm

        right_name <- rev(colnames(plot_coef_mat_2))
        right_name_num <- length(right_name)

        left_name <- c()
        for (row_ind in 1:nrow(plot_coef_mat)){
          left_name <- append(left_name, strsplit(rownames(plot_coef_mat)[row_ind], split = "neig_")[[1]][2])
        }
        rownames(plot_coef_mat_2) <- left_name

        g <- igraph::graph_from_incidence_matrix(plot_coef_mat_2, weighted = T)

        # Set coordinate
        igraph::V(g)$x <- ifelse(names(igraph::V(g)) %in% right_name, 2, 0)
        igraph::V(g)$y <- ifelse(igraph::V(g)$x == 0, rev(1:g_row_num), -1)

        for (right_ind in 1:right_name_num){
          igraph::V(g)$y[names(igraph::V(g)) == right_name[right_ind]] <- ((g_row_num - 1) / (right_name_num - 1)) * (right_ind - 1) + 1
        }

        # Set color for positive or negative
        igraph::E(g)$color <- ifelse(igraph::E(g)$weight > 0, "red", "blue")
        igraph::E(g)$width <- abs(igraph::E(g)$weight)

        bi_graph_list[[cell_type_ind]] <- g

        vsz_1 <- 4 * max(nchar(left_name)) + 20
        vsz_2 <- 80 / nrow(plot_coef_mat_2)
        vlc <- 10 / max(nrow(plot_coef_mat_2), ncol(plot_coef_mat_2))

        grDevices::png(paste0(output_dir, "/cell_type_", cell_type_ind, "_bipartite_graph.png"))
        plot(g,
             vertex.color = "white",
             vertex.frame.color = "black",
             vertex.label.color = "black",
             vertex.shape = "rectangle",
             vertex.size = vsz_1,
             vertex.size2 = vsz_2,
             vertex.label.cex = vlc)
        grDevices::dev.off()

        grDevices::pdf(paste0(output_dir, "/cell_type_", cell_type_ind, "_bipartite_graph.pdf"))
        plot(g,
             vertex.color = "white",
             vertex.frame.color = "black",
             vertex.label.color = "black",
             vertex.shape = "rectangle",
             vertex.size = vsz_1,
             vertex.size2 = vsz_2,
             vertex.label.cex = vlc)
        grDevices::dev.off()

      } else {

        sink(paste0(output_dir, "/cell_type_", cell_type_ind, "_bipartite_graph.txt"))
        print(paste0("No significant combinations of genes and neighboring cell types in ", res.sel.var$cell_type_list[[cell_type_ind]]))
        sink()

      }

    } else {

      sink(paste0(output_dir, "/cell_type_", cell_type_ind, "_bipartite_graph.txt"))
      print(paste0("NULL returned in ", res.sel.var$cell_type_list[[cell_type_ind]]))
      sink()

    }

  }

  finish_time <- Sys.time()
  print(paste0("=== cellCellGraph finished... ", Sys.time(), " ==="))
  print(paste0("=========="))

  return(list(bi_graph_list = bi_graph_list))

}
