#' @noRd
is_integer <- function(x) {
    identical(x, round(x))
}

#' @noRd
zero_diag <- function(mat) {
    diag(mat) <- 0
    mat
}

#' @noRd
img <- function(mat, ..., scl=0) {
    mat <- as.matrix(mat)
    diag(mat) <- diag(mat) * scl
    graphics::image(mat, ...)
}

#' @noRd
proportion_non_zero <- function(mat) {
    row_props <- Matrix::rowSums(mat != 0) / ncol(mat)
    return(row_props)
}

#' Subset Rows on Proportion Non-Zero
#'
#' @param mat A numerical matrix
#'
#' @param cutoff Threshold value in range \[0, 1\]; proportion of each row
#'   required to be non-zero
#'
#' @examples {
#' lst_scrna <- CytoTalk::scrna_cyto
#' cell_type_a <- "Macrophages"
#' cutoff_a <- 0.8
#' mat_a <- extract_group(cell_type_a, lst_scrna)
#' result <- subset_non_zero(mat_a, cutoff_a)
#' }
#'
#' @return A matrix
#'
#' @export
subset_non_zero <- function(mat, cutoff) {
    index <- (cutoff <= proportion_non_zero(mat))
    mat[index, ]
}

#' @noRd
subset_non_zero_old <- function(mat, cutoff) {
    thresh <- floor(ncol(mat) * cutoff)
    index <- thresh <= Matrix::rowSums(mat != 0)
    mat[index, ]
}

#' Subset Rows on Rownames
#'
#' @param mat A numerical matrix
#'
#' @param labels A subset of the matrix's rownames to filter by
#'
#' @examples {
#' lst_scrna <- CytoTalk::scrna_cyto
#' pcg <- CytoTalk::pcg_human
#' cell_type_a <- "Macrophages"
#' mat_a <- extract_group(cell_type_a, lst_scrna)
#' result <- subset_rownames(mat_a, pcg)
#' }
#'
#' @return A matrix
#'
#' @export
subset_rownames <- function(mat, labels) {
    index <- which(!is.na(cmatch(rownames(mat), labels)))
    mat[index, ]
}

#' @noRd
add_noise <- function(mat) {
    n <- nrow(mat)
    m <- ncol(mat)
    dummy <- c(1e-20, rep(0, m - 1))
    mat <- mat + matrix(replicate(n, sample(dummy)), n, byrow = TRUE)
    Matrix::Matrix(mat)
}

#' @noRd
extract_group_basic <- function(group, mat, labels) {
    mat[, labels == group, drop = FALSE]
}

#' Extract Cell Type from Meta List
#'
#' @param group Character vector with a single string, name of cell type to
#'   extract
#'
#' @param lst A meta matrix list (outputted from `doc_fileio` methods)
#'
#' @examples {
#' lst_scrna <- CytoTalk::scrna_cyto
#' cell_type_a <- "Macrophages"
#' result <- extract_group(cell_type_a, lst_scrna)
#' }
#'
#' @return A matrix
#'
#' @export
extract_group <- function(group, lst) {
    extract_group_basic(group, lst[[1]], lst[[2]])
}

#' @noRd
group_meta_basic <- function(mat, labels) {
    groups <- sort(unique(labels))
    lst <- lapply(groups, extract_group_basic, mat, labels)
    names(lst) <- groups
    lst
}

#' Meta List to Named List
#'
#' @param lst A meta matrix list (outputted from `doc_fileio` methods)
#'
#' @examples {
#' lst_scrna <- CytoTalk::scrna_cyto
#' result <- group_meta(lst_scrna)
#' }
#'
#' @return A list of matrices
#'
#' @export
group_meta <- function(lst) {
    group_meta_basic(lst[[1]], lst[[2]])
}

#' Named List to Meta List
#'
#' @param lst A meta matrix list (outputted from `doc_fileio` methods)
#'
#' @examples {
#' lst_scrna <- CytoTalk::scrna_cyto
#' lst_group <- group_meta(lst_scrna)
#' result <- ungroup_meta(lst_group)
#' }
#'
#' @return A list of a matrix and a meta vector
#'
#' @export
ungroup_meta <- function(lst) {
    ncols <- vapply(lst, ncol, numeric(1))
    cell_types <- unlist(lapply(seq_len(length(ncols)), function(i) {
        rep(names(ncols)[i], ncols[i])
    }))
    new_named_list(do.call(cbind, lst), cell_types)
}

#' Find Ligand Receptor Pairs in Matrix
#'
#' @param mat A numerical matrix
#'
#' @param lrp Ligand-receptor pair data.frame; see `CytoTalk::lrp_human` for
#'   an example
#'
#' @examples {
#' lst_scrna <- CytoTalk::scrna_cyto
#' lrp <- CytoTalk::lrp_human
#' cell_type_a <- "Macrophages"
#' mat_a <- extract_group(cell_type_a, lst_scrna)
#' result <- match_lr_pairs(mat_a, lrp)
#' }
#'
#' @return A matrix with two columns (ligand-receptor)
#'
#' @export
match_lr_pairs <- function(mat, lrp) {
    # save a copy of the rownames
    hold <- rownames(mat)

    # match lr_pairs to mat rownames
    index <- data.frame(apply(lrp, 2, function(x) {
        match(toupper(x), toupper(hold))
    }))
    index <- index[rowSums(is.na(index)) == 0, ]

    # return out
    data.frame(
        ligand = hold[index[, 1]],
        receptor = hold[index[, 2]]
    )
}

#' scRNAseq Count Normalization
#'
#' @param mat An integer matrix
#'
#' @param scale.factor A single number constant by which to scale the
#'   transformed data by
#'
#' @examples {
#' lst_scrna <- CytoTalk::scrna_cyto
#' cell_type_a <- "Macrophages"
#' mat_a <- extract_group(cell_type_a, lst_scrna)
#' result <- normalize_sparse(mat_a)
#' }
#'
#' @return A matrix
#'
#' @export
normalize_sparse <- function(mat, scale.factor=10000) {
    log1p(Matrix::t(Matrix::t(mat) / Matrix::colSums(mat) * scale.factor))
}

#' Transform Count Data if Detected
#'
#' @param mat A numerical matrix
#'
#' @param auto_transform Should the data be transformed if counts are detected?
#'
#' @examples {
#' lst_scrna <- CytoTalk::scrna_cyto
#' cell_type_a <- "Macrophages"
#' mat_a <- extract_group(cell_type_a, lst_scrna)
#' result <- check_count_data(mat_a)
#' }
#'
#' @return A matrix
#'
#' @export
check_count_data <- function(mat, auto_transform=TRUE) {
    check <- is_integer(mat)
    if (check) {
        if (auto_transform) {
            mat <- normalize_sparse(mat)
            # warn that normalization was performed
            msg <- paste(
                "count data detected;",
                "auto-transformed (see `?check_count_data`)"
            )
            warnifnot(!check, msg)
        } else {
            # warn if counts detected
            msg <- "count data detected, make sure to transform it"
            warnifnot(!check, msg)
        }
    }
    # return matrix
    mat
}
