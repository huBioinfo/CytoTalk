#' @rdname doc_mutinfo
#' @export
discretize_sparse <- function(mat, disc="equalfreq", nbins=nrow(mat)^(1/3)) {
    mat_disc <- infotheo::discretize(as.matrix(mat), disc, max(2, nbins))
    mat_disc <- t(do.call(rbind, mat_disc)) - 1
    Matrix::Matrix(mat_disc, sparse = TRUE)
}

#' @rdname doc_mutinfo
#' @export
mutinfo_xy <- function(x, y, method="emp", normalize=FALSE) {
    if (!normalize) {
        return(infotheo::mutinformation(x, y, method))
    }

    h1 <- infotheo::entropy(x, method)
    h2 <- infotheo::entropy(y, method)
    h12 <- infotheo::entropy(data.frame(x, y), method)

    hm <- min(h1, h2)
    if (hm == 0) {
        return(0)
    }

    ((h1 + h2) - h12) / hm
}

#' @rdname doc_mutinfo
#' @export
mi_mat_parallel <- function(mat, method="emp", normalize=FALSE) {
    mat <- as.matrix(mat)
    n <- ncol(mat)
    ent <- apply(mat, 2, infotheo::entropy, method)

    i <- NULL
    res <- foreach::`%dopar%`(foreach::foreach(i = seq_len(n)), {
        h1 <- ent[i]
        x <- as.numeric(mat[, i])
        vapply(i:n, function(j) {
            h2 <- ent[j]
            y <- as.numeric(mat[, j])
            hm <- min(h1, h2)
            h12 <- infotheo::entropy(data.frame(x, y), method)
            if (hm == 0) {
                0
            } else if (normalize) {
                ((h1 + h2) - h12) / hm
            } else {
                ((h1 + h2) - h12)
            }
        }, numeric(1))
    })

    out <- matrix(0, n, n)
    out[lower.tri(out, diag = TRUE)] <- unlist(res)
    out + t(out * lower.tri(out))
}
