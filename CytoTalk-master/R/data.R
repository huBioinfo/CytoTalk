#' Human ligand-receptor pairs
#'
#' Contains the names of human ligand-receptor pairs
#'
#' @format A data frame with two variables: \code{ligand} and \code{receptor}
"lrp_human"

#' Mouse ligand-receptor pairs
#'
#' Contains the names of mouse ligand-receptor pairs
#'
#' @format A data frame with two variables: \code{ligand} and \code{receptor}
"lrp_mouse"

#' Human protein-coding genes (PCG)
#'
#' Contains the names of human protein-coding genes
#'
#' @format A large character vector
"pcg_human"

#' Mouse protein-coding genes (PCG)
#'
#' Contains the names of mouse protein-coding genes
#'
#' @format A large character vector
"pcg_mouse"

#' Protetta Stone
#'
#' @format A dataframe
"protetta"

#' Example scRNAseq Data (CellphoneDB)
#'
#' @format A list
"scrna_cpdb"

#' Example scRNAseq Data (CytoTalk)
#'
#' @format A list
"scrna_cyto"

#' Example Pipeline Result (CytoTalk)
#'
#' @format A list
"result_cyto"

#' Graphing with Cytoscape and Graphviz
#'
#' @param df_net Final network; for example, the output of the
#'   `extract_network` function
#'
#' @param df_net_sub A subset of the final network (pathway); for example, the
#'   output of the `extract_pathways` function
#'
#' @param cell_type_a Name of cell type A that matches scRNA-seq file; for
#'   example, `"Fibroblasts"`
#'
#' @param dir_out Folder used for output
#'
#' @examples {
#' pathways <- CytoTalk::result_cyto$pathways$raw
#' result <- graph_pathway(pathways[[1]])
#' }
#'
#' @return Graphs which represent the subset pathways of the final network
#'
#' @name doc_graphing
NULL
#> NULL

#' Prize-Collecting Steiner Tree and Pathway Extraction
#'
#' @param lst_net Integrated network
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
#' @param lst_pcst PCST output
#'
#' @param mat_pem PEM output
#'
#' @param beta A single beta value, see `beta_max` paramter for more detail
#'
#' @param omega A single omega value, see `omega_min` for more detail
#'
#' @param df_net Final network; for example, the output of the
#'   `extract_network` function
#'
#' @param cell_type_a Name of cell type A that matches scRNA-seq file; for
#'   example, `"Fibroblasts"`
#'
#' @param depth Starting at each ligand-receptor pair in the resultant network,
#'   how many steps out from that pair should be taken to generate each
#'   neighborhood?
#'
#' @examples {
#' lst_net <- CytoTalk::result_cyto$integrated_net
#' beta_max <- 100
#' omega_min <- 0.5
#' omega_max <- 0.5
#' # result <- run_pcst(lst_net, beta_max, omega_min, omega_max)
#' }
#'
#' @return Scores for the integrated network
#'
#' @name doc_pcst
NULL
#> NULL

#' Integrated Co-Expression Network
#'
#' @param mat_type A matrix of a single cell type
#'
#' @param lrp A dataframe or matrix object with two columns, ligands names and
#'   the names of their receptors; by default, uses the `lrp_human` data. This
#'   package also includes `lrp_mouse`, but you can also use your own data
#'
#' @param mat_intra_a Intracellular network for cell type A; for example, using
#'   `parmigene::aracne.m` on the result of `mi_mat_parallel` will filter out
#'   any indirect edges per the data processing inequality
#'
#' @param mat_intra_b Intracellular network for cell type B
#'
#' @param mat_intra Intracellular network of any type
#'
#' @param mat_pem PEM output
#'
#' @param cell_type_a Name of cell type A that matches scRNA-seq file; for
#'   example, `"Fibroblasts"`
#'
#' @param cell_type_b Name of cell type B that matches scRNA-seq file; for
#'   example, `"LuminalEpithelialCells"`
#'
#' @param cell_type Name of any cell type
#'
#' @param vec_relev Vector containing gene relevance scores
#'
#' @param vec_nst_a Vector containing nonselftalk-scores for cell type A

#' @param vec_nst_b Vector containing nonselftalk-scores for cell type B
#'
#' @examples {
#' lst_scrna <- CytoTalk::scrna_cyto
#' cell_type_a <- "Macrophages"
#' lrp <- CytoTalk::lrp_human
#' mat_a <- extract_group(cell_type_a, lst_scrna)
#' result <- nonselftalk(mat_a, lrp)
#' }
#'
#' @return Various outputs to build the integrated network
#'
#' @name doc_integrated
NULL
#> NULL

#' Information Theory
#'
#' @param x,y Any discrete numerical vector
#'
#' @param mat Any discrete numerical matrix, processed column-wise
#'
#' @param disc Name of discretization method. Options include "equalfreq",
#'   "equalwidth", and "globalequalwidth". See `infotheo::discretize` for more
#'   detail
#'
#' @param nbins The number of bins by which the continuous data should be
#'   discretized; by default the cube root of the number of samples
#'
#' @param method Name of entropy estimator. Options include "emp", "mm",
#'   "shrink", and "sg". See `infotheo::mutinformation` for more detail
#'
#' @param normalize Should the entropies be normalized to fall within the range
#'   \[0, 1\]?
#'
#' @examples {
#' lst_scrna <- CytoTalk::scrna_cyto
#' pcg <- CytoTalk::pcg_human
#' cell_type_a <- "Macrophages"
#' cutoff_a <- 0.8
#' mat_a <- extract_group(cell_type_a, lst_scrna)
#' mat_filt_a <- subset_rownames(subset_non_zero(mat_a, cutoff_a), pcg)
#' mat_disc_a <- discretize_sparse(Matrix::t(mat_filt_a))
#' result <- mi_mat_parallel(mat_disc_a, method = "mm")
#' }
#'
#' @return Mutual information calculations for filtered matrices
#'
#' @name doc_mutinfo
NULL
#> NULL

#' Registering Parallel Backend
#'
#' @param cores How many cores to use for parallel processing?
#'
#' @examples {
#' register_parallel(cores=2)
#' unregister_parallel()
#' }
#'
#' @return Nothing
#'
#' @name doc_parallel
NULL
#> NULL

#' Preferential Expression Measure
#'
#' @param mat_scrna Matrix containing scRNA-seq data, along with `labels`,
#'   defines the transformed count data and the associated cell types
#'
#' @param labels Associated cell types, column-wise, of `mat_scrna`
#'
#' @param lst_scrna List containing scRNA-seq data; for example, lists returned
#'   from `read_matrix_folder` or `read_matrix_with_meta`
#'
#' @examples {
#' lst_scrna <- CytoTalk::scrna_cyto
#' result <- pem(lst_scrna)
#' }
#'
#' @return Preferential expression measure scores
#'
#' @name doc_pem
NULL
#> NULL
