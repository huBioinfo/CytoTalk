#' Main CytoTalk Pipeline
#'
#' @param lst_scrna List containing scRNA-seq data; for example, lists returned
#'   from `read_matrix_folder` or `read_matrix_with_meta`
#'
#' @param cell_type_a Name of cell type A that matches scRNA-seq file; for
#'   example, `"Fibroblasts"`
#'
#' @param cell_type_b Name of cell type B that matches scRNA-seq file; for
#'   example, `"LuminalEpithelialCells"`
#'
#' @param dir_out Folder used for output; if not specified, a "CytoTalk-output"
#'   folder will be generated
#'
#' @param cutoff_a Proportional threshold for lowly expressed genes in cell
#'   type A (range of \[0-1\]); for example, 0.2 means genes with some
#'   expression in at least 20% of cells are retained
#'
#' @param cutoff_b Proportional expression threshold for cell type B (range of
#'   \[0-1\])
#'
#' @param pcg A character vector, contains the names of protein coding genes;
#'   by default, uses the `pcg_human` data. This package also includes
#'   `pcg_mouse`, but you can also use your own data
#'
#' @param lrp A dataframe or matrix object with two columns, ligands names and
#'   the names of their receptors; by default, uses the `lrp_human` data. This
#'   package also includes `lrp_mouse`, but you can also use your own data
#'
#' @param beta_max Upper limit of the test values of the PCSF objective
#'   function parameter $I^2$, which is inversely proportional to the total
#'   number of genes in a given cell-type pair; suggested to be 100 (default)
#'   if the total number of genes in a given cell-type pair is above 10,000; if
#'   the total number of genes is below 5,000, increase to 500
#'
#' @param omega_min Start point of omega range; omega represents the edge cost
#'   of the artificial network, but has been found to be less significant than
#'   beta. Recommended minimum of `0.5`
#'
#' @param omega_max End point of range between `omega_min` and `omega_max`,
#'   step size of `0.1`. Recommended maximum of `1.5`
#'
#' @param depth Starting at each ligand-receptor pair in the resultant network,
#'   how many steps out from that pair should be taken to generate each
#'   neighborhood?
#'
#' @param ntrial How many random network subsets shall be created to get an
#'   empirical p-value for node prize and edge cost?
#'
#' @param cores How many cores to use for parallel processing?
#'
#' @param echo Should update messages be printed?
#'
#' @examples {
#' cell_type_a <- "Macrophages"
#' cell_type_b <- "LuminalEpithelialCells"
#' cutoff_a <- 0.6
#' cutoff_b <- 0.6
#' # result <- CytoTalk::run_cytotalk(CytoTalk::scrna_cyto,
#' #                                  cell_type_a, cell_type_b,
#' #                                  cutoff_a, cutoff_b,
#' #                                  cores = 2)
#' }
#'
#' @return A list containing model parameters, prefential expression measure,
#' the integrated co-expression network, the results of the PCST, and resulting
#' pathways from the final extracted network
#'
#' @export
run_cytotalk <- function(
    lst_scrna, cell_type_a, cell_type_b,
    cutoff_a=0.2, cutoff_b=0.2,
    pcg=CytoTalk::pcg_mouse, lrp=CytoTalk::lrp_mouse,
    beta_max=100, omega_min=0.5, omega_max=0.5,
    depth=3, ntrial=1000,
    cores=NULL, echo=TRUE, dir_out=NULL,
    use_cache = TRUE, cache_dir = "C:/Users/Administrator/Desktop/temp2_Data") {

    # save numeric parameters
    params <- list(
        cell_type_a = cell_type_a, cell_type_b = cell_type_b,
        cutoff_a = cutoff_a, cutoff_b = cutoff_b,
        beta_max = beta_max, omega_min = omega_min, omega_max = omega_max,
        depth = depth, ntrial = ntrial
    )

    # register parallel backend
    if (is.null(cores) || 1 < cores) {
        unregister_parallel()
        register_parallel(cores)
    }

    # create directory
    if (!is.null(dir_out) && !dir.exists(dir_out)) {
        dir.create(dir_out, recursive = TRUE)
    }

    # --- 新增：处理缓存目录 ---
    if (use_cache) {
      if (!dir.exists(cache_dir)) {
        dir.create(cache_dir, recursive = TRUE)
      }
      # 为不同的细胞类型组合创建唯一的缓存文件名，避免冲突
      cache_suffix <- paste0(cell_type_a, "_vs_", cell_type_b)
      # 合并后的预处理缓存，包含 mat_pem, mat_a, mat_b, vec_nst, mat_filt, mat_disc
      pre_cache_file <- file.path(cache_dir, paste0("preprocessing_", cache_suffix, ".rda"))
      # MI 矩阵缓存
      mi_a_cache_file <- file.path(cache_dir, paste0("mat_mi_a_", cache_suffix, ".rda"))
      mi_b_cache_file <- file.path(cache_dir, paste0("mat_mi_b_", cache_suffix, ".rda"))
    }
    
    if (echo) {
      tick(1, "Preprocessing...")
    }
    
    # --- 修改：缓存 PEM, NST, Filtered, 和 Discretized matrices ---
    if (use_cache && file.exists(pre_cache_file)) {
      if (echo) cat("Loading preprocessing results from cache...\n")
      load(pre_cache_file)
    } else {
      # 计算 PEM
      mat_pem <- pem(lst_scrna)
      # 提取细胞类型矩阵
      mat_a <- extract_group(cell_type_a, lst_scrna)
      mat_b <- extract_group(cell_type_b, lst_scrna)
      errorifnot(ncol(mat_a)>0 & ncol(mat_b)>0, "No corresponding data for the celltype")
      # 计算 NST
      vec_nst_a <- nonselftalk(mat_a, lrp)
      vec_nst_b <- nonselftalk(mat_b, lrp)
      # 过滤基因
      mat_filt_a <- subset_rownames(subset_non_zero_old(mat_a, cutoff_a), pcg)
      mat_filt_b <- subset_rownames(subset_non_zero_old(mat_b, cutoff_b), pcg)
      # 离散化
      mat_disc_a <- discretize_sparse(Matrix::t(mat_filt_a))
      mat_disc_b <- discretize_sparse(Matrix::t(mat_filt_b))
      
      if (use_cache) {
        if (echo) cat("Saving preprocessing results to cache...\n")
        save(mat_pem, mat_a, mat_b, vec_nst_a, vec_nst_b, mat_filt_a, mat_filt_b, mat_disc_a, mat_disc_b, file = pre_cache_file)
      }
    }
    

    # write out PEM matrix
    if (!is.null(dir_out)) {
        fpath <- file.path(dir_out, "PEM.txt")
        vroom_write_silent(mat_pem, fpath, rownames = TRUE)
    }

    
    if (echo) {
      tick(2, "Mutual information matrix...")
    }
    
    # --- 修改：缓存 mat_mi_a, mat_mi_b ---
    # mat_disc_a 和 mat_disc_b 已经从上面的缓存中加载
    if (use_cache && file.exists(mi_a_cache_file) && file.exists(mi_b_cache_file)) {
      if (echo) cat("Loading Mutual Information matrices from cache...\n")
      load(mi_a_cache_file)
      load(mi_b_cache_file)
    } else {
      mat_mi_a <- mi_mat_parallel(mat_disc_a, method = "mm")
      mat_mi_b <- mi_mat_parallel(mat_disc_b, method = "mm")
      
      if (use_cache) {
        if (echo) cat("Saving Mutual Information matrices to cache...\n")
        save(mat_mi_a, file = mi_a_cache_file)
        save(mat_mi_b, file = mi_b_cache_file)
      }
    }
    
    dimnames(mat_mi_a) <- list(colnames(mat_disc_a), colnames(mat_disc_a))
    dimnames(mat_mi_b) <- list(colnames(mat_disc_b), colnames(mat_disc_b))
    
    # --- 新增：定义网络构建后的缓存文件 ---
    if (use_cache) {
      # 这个缓存文件将保存 aracne 过滤和归一化后的矩阵
      network_cache_file <- file.path(cache_dir, paste0("network_", cache_suffix, ".rda"))
    }
    

    if (echo) {
      tick(3, "Indirect edge-filtered network...")
    }
    
    # --- 修改：缓存 aracne 过滤和归一化后的矩阵 ---
    if (use_cache && file.exists(network_cache_file)) {
      if (echo) cat("Loading network construction results from cache...\n")
      load(network_cache_file)
    } else {
      # 如果缓存不存在，则执行计算
      mat_intra_a <- Matrix::Matrix(parmigene::aracne.m(zero_diag(mat_mi_a)))
      mat_intra_b <- Matrix::Matrix(parmigene::aracne.m(zero_diag(mat_mi_b)))
      
      # remove some independent genes
      inde_gene_id_a = which(Matrix::rowSums(mat_intra_a)<=0)
      if(length(inde_gene_id_a)>0){
        mat_intra_a = mat_intra_a[-inde_gene_id_a,-inde_gene_id_a]
      }
      inde_gene_id_b = which(Matrix::rowSums(mat_intra_b)<=0)
      if(length(inde_gene_id_b)>0){
        mat_intra_b = mat_intra_b[-inde_gene_id_b,-inde_gene_id_b]
      }
      
      # 如果启用了缓存，则保存计算结果
      if (use_cache) {
        if (echo) cat("Saving network construction results to cache...\n")
        save(mat_intra_a, mat_intra_b, file = network_cache_file)
      }
    }
    
    
    if (echo) {
        tick(4, "Integrate network...")
    }

    lst_net <- integrate_network(
        vec_nst_a, vec_nst_b, mat_intra_a, mat_intra_b,
        cell_type_a, cell_type_b, mat_pem, mat_a, lrp
    )

    # write out integrated nodes and edges
    if (!is.null(dir_out)) {
        fpath <- file.path(dir_out, "IntegratedNodes.txt")
        vroom_write_silent(lst_net$nodes, fpath)
        fpath <- file.path(dir_out, "IntegratedEdges.txt")
        vroom_write_silent(lst_net$edges, fpath)
    }

    if (echo) {
        tick(5, "PCSF...")
    }

    lst_pcst <- run_pcst(lst_net, beta_max, omega_min, omega_max)

    # write out PCST nodes and edges
    if (!is.null(dir_out)) {
        fpath <- file.path(dir_out, "PCSTNodeOccurance.txt")
        vroom_write_silent(lst_pcst$nodes, fpath)
        fpath <- file.path(dir_out, "PCSTEdgeOccurance.txt")
        vroom_write_silent(lst_pcst$edges, fpath)
    }

    if (echo) {
        tick(6, "Determine best signaling network...")
    }

    df_test <- ks_test_pcst(lst_pcst)
    index <- order(as.numeric(df_test[, "pval"]))[1]
    beta <- df_test[index, "beta"]
    omega <- df_test[index, "omega"]

    # write out PCST scores
    if (!is.null(dir_out)) {
        vroom_write_silent(df_test, file.path(dir_out, "PCSTScores.txt"))
    }

    if (echo) {
        tick(7, "Generate network output...")
    }

    df_net <- extract_network(lst_net, lst_pcst, mat_pem, beta, omega)
    lst_path <- extract_pathways(df_net, cell_type_a, depth)
    lst_graph <- lapply(lst_path, graph_pathway)

    # no pathways found
    if (is.null(lst_path)) {
        result <- list(
            params = params,
            pem = mat_pem,
            integrated_net = lst_net,
            pcst = list(
                occurances = lst_pcst,
                ks_test_pval = df_test,
                final_network = df_net
            ),
            pathways = NULL
        )

        # unregister parallel backend
        if (is.null(cores) || 1 < cores) {
            unregister_parallel()
        }

        tick(8, "NOTE: No pathways found, analysis skipped!")
        return(result)
    }

    # write out pathways
    if (!is.null(dir_out)) {
        dir_path <- file.path(dir_out, "pathways")
        if (!dir.exists(dir_path)) {
            dir.create(dir_path, recursive = TRUE)
        }
        fnames <- names(lst_path)
        for (fn in fnames) {
            fpath <- file.path(dir_path, sprintf("%s.txt", fn))
            vroom_write_silent(lst_path[[fn]], fpath)
        }

        dir_gv <- file.path(dir_out, "graphviz")
        if (!dir.exists(dir_gv)) {
            dir.create(dir_gv, recursive = TRUE)
        }
        fnames <- names(lst_graph)
        for (fn in fnames) {
            fpath <- file.path(dir_gv, sprintf("%s.svg", fn))
            content <- DiagrammeRsvg::export_svg(lst_graph[[fn]])
            write(content, fpath)
        }
    }

    # write out final network
    if (!is.null(dir_out)) {
        write_network_sif(df_net, cell_type_a, dir_out)
        vroom_write_silent(df_net, file.path(dir_out, "FinalNetwork.txt"))
    }

    if (echo) {
        tick(8, "Analyze pathways...")
    }

    lst_pval <- lapply(
        lst_path, analyze_pathway, lst_net,
        cell_type_a, cell_type_b, beta, ntrial
    )

    # format the ligand and receptor cell types
    nodes <- do.call(rbind, strsplit(names(lst_path), "--"))
    df_pval <- do.call(rbind, apply(nodes, 1, function(x) {
        t <- ifelse(endsWith(x, "_"), cell_type_b, cell_type_a)
        x <- gsub("_$", "", x)
        data.frame(
            ligand = x[1], receptor = x[2],
            ligand_type = t[1], receptor_type = t[2]
        )
    }))

    # combine with scores and sort
    df_pval <- cbind(df_pval, do.call(rbind, lst_pval))
    df_pval <- df_pval[order(as.numeric(df_pval$pval_potential)), ]

    # write out analysis
    if (!is.null(dir_out)) {
        fpath <- file.path(dir_out, "PathwayAnalysis.txt")
        vroom_write_silent(df_pval, fpath)
    }

    # return out
    result <- list(
        params = params,
        pem = mat_pem,
        integrated_net = lst_net,
        pcst = list(
            occurances = lst_pcst,
            ks_test_pval = df_test,
            final_network = df_net
        ),
        pathways = list(
            raw = lst_path,
            graphs = lst_graph,
            df_pval = df_pval
        )
    )

    # write out session
    if (!is.null(dir_out)) {
        fpath <- file.path(dir_out, "CytoTalkSession.rda")
        save(result, file = fpath, version = 2)
    }

    # unregister parallel backend
    if (is.null(cores) || 1 < cores) {
        unregister_parallel()
    }

    result
}
