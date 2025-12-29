#' @noRd
cmatch <- function(x, y) {
    match(toupper(x), toupper(y))
}

#' @noRd
all_identical <- function(lst) {
    all(vapply(lst[-1], FUN = identical, logical(1), lst[[1]]))
}

#' @noRd
format_message <- function(msg) {
    sprintf("%s\n", trimws(msg))
}

#' @noRd
errorifnot <- function(expr, msg) {
    if (!expr) stop(format_message(msg))
}

#' @noRd
warnifnot <- function(expr, msg) {
    if (!expr) warning(format_message(msg))
}

#' @noRd
is_integer <- function(mat) {
    identical(mat, round(mat))
}

#' @noRd
zero_diag <- function(mat) {
    diag(mat) <- 0
    mat
}

#' @noRd
dir_full <- function(folder, ...) {
    fnames <- dir(folder, ...)
    file.path(folder, fnames)
}

#' @noRd
totitle <- function(x) {
  x <- as.character(x)
  is_all_upper <- x == toupper(x)
  
  out <- x
  idx <- !is_all_upper
  out[idx] <- gsub("^([A-Za-z])", "\\U\\1", tolower(out[idx]), perl = TRUE)
  
  out
}

#' @noRd
add_suffix <- function(x, type) {
    sprintf("%s__%s", totitle(x), type)
}

#' @noRd
minmax <- function(x, ...) {
    (x - min(x, ...)) / (max(x, ...) - min(x, ...))
}

#' @noRd
zero_na_neg <- function(x) {
    ifelse(x < 0 | is.na(x), 0, x)
}

#' @noRd
now <- function() {
    format(Sys.time(), "%H:%M:%S")
}

#' @noRd
tick <- function(step, msg) {
    cat("[", step, " / 8] (", now(), ") ", msg, "\n", sep = "")
}
