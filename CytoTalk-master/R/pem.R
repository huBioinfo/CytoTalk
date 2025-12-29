#' @rdname doc_pem
#' @export
pem_basic <- function(mat_scrna, labels) {
    # exponentiate the normalized data
    mat_scrna <- expm1(mat_scrna)
    # seperate out the cell types,
    # rowmeans per scRNA file
    lst_means <- lapply(group_meta_basic(mat_scrna, labels), Matrix::rowMeans)
    # sums of rowmeans, per scRNA file
    lst_sums <- lapply(lst_means, sum)
    # sums of rowmeans, per gene (combined scRNA files)
    vec_gene_sums <- rowSums(do.call(cbind, lst_means))
    # overall sum of rowmeans
    total_sum <- sum(vec_gene_sums)

    # for all cell_types
    mat_pem <- list()
    for (i in seq_len(length(lst_means))) {
        mat_pem[[i]] <- vector()

        # what proportion of this cell type's rowmean sum
        # accounts for the whole?
        cell_type_prop <- lst_sums[[i]] / total_sum

        # for all genes
        for (j in seq_len(length(vec_gene_sums))) {

            # scale gene sum to cell type proportion
            gene_prop <- vec_gene_sums[j] * cell_type_prop
            # scale cell type rowmean to gene proportion
            mat_pem[[i]][j] <- log10(lst_means[[i]][j] / gene_prop)
        }
    }

    # join columns to dataframe
    mat_pem <- as.data.frame(do.call(cbind, mat_pem))

    # extract valid type names, copy over rownames
    colnames(mat_pem) <- names(lst_means)
    rownames(mat_pem) <- rownames(mat_scrna)

    # return pem matrix
    as.matrix(mat_pem)
}

#' @rdname doc_pem
#' @export
pem <- function(lst_scrna) {
    pem_basic(lst_scrna[[1]], lst_scrna[[2]])
}
