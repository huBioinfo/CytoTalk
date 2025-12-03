#' @noRd
score_subnetwork <- function(prizes, costs, beta) {
    potential <- (beta * sum(prizes) - sum(costs))
    c(mean(prizes), mean(costs), potential)
}

#' @noRd
subsample_network_shuffle <- function(prizes, costs, n_nodes, n_edges) {
    list(prize = sample(prizes, n_nodes), cost = sample(costs, n_edges))
}

#' @noRd
gamma_fit <- function(x, y) {
    x <- x - min(x) + y
    shape <- mean(x)^2 / stats::var(x)
    scale <- stats::var(x) / mean(x)
    list(
        params = c(shape, scale),
        kstest = stats::ks.test(x, "pgamma", shape = shape, scale = scale)
    )
}

#' @noRd
gamma_score <- function(x, inverse=FALSE) {
    if (inverse) { x <- (-x) }

    v <- x[-1]
    y <- stats::optimize(function(y) {
        -gamma_fit(v, y)$kstest$p.value
    }, c(0, 1000))$minimum

    fit <- gamma_fit(v, y)
    shape <- fit$params[1]
    scale <- fit$params[2]

    x1 <- x[1] - min(v) + y
    stats::pgamma(x1, shape = shape, scale = scale, lower.tail = FALSE)
}

#' @noRd
convert_names <- function(x, cell_type_a, cell_type_b) {
    suffix_a <- sprintf("__%s", cell_type_a)
    suffix_b <- sprintf("__%s", cell_type_b)
    gsub(suffix_a, "", gsub(suffix_b, "_", x))
}

#' Final Network Pathway Analysis
#'
#' @param df_net_sub A subset of the final network (pathway); for example, the
#'   output of the `extract_pathways` function
#'
#' @param lst_net Integrated network
#'
#' @param cell_type_a Name of cell type A that matches scRNA-seq file; for
#'   example, `"Fibroblasts"`
#'
#' @param cell_type_b Name of cell type B that matches scRNA-seq file; for
#'   example, `"LuminalEpithelialCells"`
#'
#' @param beta Upper limit of the test values of the PCSF objective function
#'   parameter $I^2$, which is inversely proportional to the total number of
#'   genes in a given cell-type pair; suggested to be 100 (default) if the
#'   total number of genes in a given cell-type pair is above 10,000; if the
#'   total number of genes is below 5,000, increase to 500
#'
#' @param ntrial How many empirical simulations to run? (Sample used to form
#'   theoretical Gamma distribution)
#'
#' @examples {
#' df_net_sub <- result_cyto$pathways$raw[[1]]
#' lst_net <- result_cyto$integrated_net
#' cell_type_a <- "Macrophages"
#' cell_type_b <- "LuminalEpithelialCells"
#' beta <- 20
#' ntrial <- 1000
#' result <- analyze_pathway(df_net_sub, lst_net, cell_type_a, cell_type_b, beta, ntrial)
#' }
#'
#' @return A data-frame containing information relating to pathway size, mean
#' node prize, mean edge cost, potential scores, and p-values from a fitted
#' Gamma distribution
#'
#' @export
analyze_pathway <- function(
    df_net_sub, lst_net, cell_type_a, cell_type_b, beta, ntrial) {

    # extract node and tables
    df_net_nodes <- lst_net$nodes[, c(1, 2)]
    df_net_edges <- lst_net$edges[, c(1, 2, 3)]

    # set columns numeric
    df_net_nodes$prize <- as.numeric(df_net_nodes$prize)
    df_net_edges$cost <- as.numeric(df_net_edges$cost)

    # convert names
    df_net_nodes$node <- convert_names(
        df_net_nodes$node, cell_type_a, cell_type_b)
    df_net_edges$node1 <- convert_names(
        df_net_edges$node1, cell_type_a, cell_type_b)
    df_net_edges$node2 <- convert_names(
        df_net_edges$node2, cell_type_a, cell_type_b)

    # prepare nodes
    df_node <- data.frame(
        node = c(df_net_sub$node1, df_net_sub$node2),
        prize = c(df_net_sub$node1_prize, df_net_sub$node2_prize)
    )
    df_node <- df_node[!duplicated(df_node[, 1]), ]
    n_nodes <- nrow(df_node)

    # prepare edges
    df_edge <- df_net_sub[, c("node1", "node2", "cost")]
    n_edges <- nrow(df_edge)

    # begin scores
    scores <- score_subnetwork(df_node[, 2], df_edge[, 3], beta)

    # simulate random subsets
    for (i in seq_len(ntrial)) {
        lst <- subsample_network_shuffle(
            df_net_nodes[, 2], df_net_edges[, 3], n_nodes, n_edges)
        scores <- rbind(scores, score_subnetwork(lst[[1]], lst[[2]], beta))
    }

    # new row of data
    data.frame(
        num_edges = n_edges, num_nodes = nrow(df_node),
        mean_prize = scores[1, 1], mean_cost = scores[1, 2],
        potential = scores[1, 3],
        pval_prize = gamma_score(scores[, 1]),
        pval_cost = gamma_score(scores[, 2], inverse = TRUE),
        pval_potential = gamma_score(scores[, 3])
    )
}
