# CONSTANTS


HEX <- c(0:9, LETTERS[seq_len(6)])
URL_GENECARDS <- "https://www.genecards.org/cgi-bin/carddisp.pl?gene="
URL_WIKIPI <- "https://hagrid.dbmi.pitt.edu/wiki-pi/index.php/search?q="
FORM_NODE <- paste0(
    "\"%s\" [label = \"%s\" href = \"",
    URL_GENECARDS,
    "%s\" width = %s height = %s fillcolor = \"%s\"]"
)
FORM_EDGE <- paste0(
  "\"%s\" -> \"%s\" [href=\"",
  URL_WIKIPI,
  "%s+%s\" penwidth=%s style=%s color=\"%s\"]"
)
FORM_GV <- trimws("
digraph {\n
pad=0.25
layout=dot
labeljust=l
splines=true
rankdir=TB
ranksep=1
nodesep=0.1
compound=true
outputorder=\"edgesfirst\"\n
node [shape=oval fixedsize=true target=\"_blank\"]
node [fontname=\"Arial\" fontsize=9 style=filled ordering=out]
edge [arrowhead=none target=\"_blank\"]\n
subgraph cluster0 {\n
margin=20
color=none
style=filled
fillcolor=\"#EEEEEE\"\n
%s\n
}\n
subgraph cluster1 {\n
margin=20
color=none
style=filled
fillcolor=\"#EEEEEE\"\n
%s\n
}\n
// cluster external horizontal order
%s
// cluster external
edge [color=limegreen arrowhead=normal arrowtail=normal dir=back]
%s\n
}")


# FUNCTIONS


#' @rdname doc_graphing
#' @export
write_network_sif <- function(df_net, cell_type_a, dir_out) {
    # format dir path
    dir_out_cs <- file.path(dir_out, "cytoscape")

    # make sure it exists
    if (!dir.exists(dir_out_cs)) {
        dir.create(dir_out_cs, recursive = TRUE)
    }

    # format filepaths
    fpath_node <- file.path(dir_out_cs, "CytoscapeNodes.txt")
    fpath_edge <- file.path(dir_out_cs, "CytoscapeEdges.txt")
    fpath_sif <- file.path(dir_out_cs, "CytoscapeNetwork.sif")

    # make names unique
    df_net$node1 <- paste0(
        df_net$node1, ifelse(df_net$node1_type == cell_type_a, "", "_"))
    df_net$node2 <- paste0(
        df_net$node2, ifelse(df_net$node2_type == cell_type_a, "", "_"))

    # format edges
    edge_types <- ifelse(df_net$is_ct_edge, "pr", "pp")
    edge_names <- sprintf("%s (%s) %s", df_net$node1, edge_types, df_net$node2)
    edges <- gsub("[()]", "", edge_names)

    # create node table
    index1 <- c("node1", "node1_type", "node1_prize", "node1_pem")
    index2 <- c("node2", "node2_type", "node2_prize", "node2_pem")
    df_node <- data.frame(rbind(
        as.matrix(df_net[, index1]),
        as.matrix(df_net[, index2])
    ))

    # naming and type conversion
    names(df_node) <- c("node", "type", "prize", "pem")
    df_node <- utils::type.convert(df_node, as.is = TRUE)

    # create edge table
    df_edge <- cbind(edge_names, df_net[, c("cost", "is_ct_edge")])
    names(df_edge) <- c("edge", "cost", "is_ct_edge")

    # write out sif
    writeLines(paste(edges, collapse = "\n"), fpath_sif)

    # write out node and edge tables
    vroom_write_silent(df_node, fpath_node)
    vroom_write_silent(df_edge, fpath_edge)
    NULL
}

#' @rdname doc_graphing
#' @export
graph_pathway <- function(df_net_sub) {
    # reorder nodes
    df_cpy <- df_net_sub
    index_lt <- (df_net_sub$node1_pem < df_net_sub$node2_pem)
    index1 <- c("node1", "node1_type", "node1_prize", "node1_pem")
    index2 <- c("node2", "node2_type", "node2_prize", "node2_pem")
    for (i in seq_len(length(index_lt))) {
        if (index_lt[i]) {
            df_net_sub[i, index2] <- df_cpy[i, index1]
            df_net_sub[i, index1] <- df_cpy[i, index2]
        }
    }

    # reorder edges
    df_net_sub <- df_net_sub[order(as.numeric(df_net_sub$cost)), ]

    # normalize PEM
    index <- c("node1_pem", "node2_pem")
    if(nrow(df_net_sub)>1){
      df_net_sub[index] <- apply(df_net_sub[index], 2, function(x) {
        minmax(ifelse(x < 0 | is.na(x), min(x, na.rm = TRUE), x))
      })
    }

    # prepare nodes
    df_node <- data.frame(rbind(
        as.matrix(df_net_sub[, index1]),
        as.matrix(df_net_sub[, index2])
    ))

    df_node <- df_node[!duplicated(df_node), ]
    df_node <- utils::type.convert(df_node, as.is = TRUE)
    names(df_node) <- c("node", "type", "prize", "pem")

    # string format nodes
    ew <- function(x) endsWith(x, "_")
    index_nodes <- ew(df_node$node)
    clean <- trimws(df_node$node, whitespace = "_")
    size <- 0.5 + 3.5 * df_node$prize^2.5
    color <- grDevices::hsv(
        ifelse(ew(df_node$node), 0.02, 0.55), (df_node$pem), 1
    )
    nodes <- sprintf(
        FORM_NODE, df_node$node, clean, clean, size, size, color)

    # string format edges
    index_edges <- ew(df_net_sub$node1) + ew(df_net_sub$node2)
    c1 <- trimws(df_net_sub$node1, whitespace = "_")
    c2 <- trimws(df_net_sub$node2, whitespace = "_")
    size <-  1.25 + 3.75 * (1 - df_net_sub$cost)^2
    # size <- ifelse(df_net_sub$is_ct_edge, 2, 1) * size
    style <- "solid"
    edge_color <- ifelse(df_net_sub$is_ct_edge, "#46E568", "black")
    edges <- sprintf(
        FORM_EDGE, df_net_sub$node1, df_net_sub$node2, c1, c2, size, style, edge_color
    )

    # cell type names
    type_a <- df_node[!ew(df_node$node), "type"][1]
    type_b <- df_node[ew(df_node$node), "type"][1]

    # string format graph
    graph <- sprintf(FORM_GV,
        paste0(
            sprintf("label=\"%s\"\ntooltip=\"%s\"\n", type_a, type_a),
            paste0(nodes[!index_nodes], collapse = "\n"), "\n",
            paste0(edges[index_edges == 0], collapse = "\n"),
            collapse = "\n"
        ),
        paste0(
            sprintf("label=\"%s\"\ntooltip=\"%s\"\n", type_b, type_b),
            paste0(nodes[index_nodes], collapse = "\n"), "\n",
            paste0(edges[index_edges == 2], collapse = "\n"),
            collapse = "\n"
        ),
        sprintf(
            "%s [style=invis, constraint=true]\n",
            gsub("\\s+\\[.+", "", edges[index_edges == 1][1])
        ),
        paste0(
            edges[index_edges == 1],
            collapse = "\n"
        )
    )

    # return out
    DiagrammeR::grViz(graph)
}
