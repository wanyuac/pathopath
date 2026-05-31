#' @title Create a network as an igraph object and identify connected components in the network based on contact lengths
#'
#' @description This function creates an igraph-class contact network from a node table and an edge table, optionally
#' pruned to remove edges with contact lengths greater than a maximum contact length, and then identifies connect
#' components in the network. Users can prune the node and edge tables in other ways before calling this function.
#'
#' @param E A mandatory edge table, which can be the data frame E in function pathopath's output network object. This
#' data frame must contain columns named "Pathway_1" and "Pathway_2".
#' @param V An optional node table, which can be the data frame V in function pathopath's output network object. This
#' data frame must contain a column named Pathway. When provided, V is passed through to the output unchanged; nodes
#' absent from the edge table are treated as singletons (isolated vertices) in the network.
#' @param l_max An optional maximum length of contacts for pruning the input network. Default: Inf (namely, no pruning).
#'
#' @return A Components object with four slots: V (the input node table, unchanged), E (edges of the possibly pruned
#' network), membership (per-node component assignments), and component_size (per-component node counts).
#'
#' @author Yu Wan, \email{yu.wan@liverpool.ac.uk}
#' @author Mohammad Saiful Islam Sajib, \email{saiful.sajib@chrfbd.org}
#'
#' @import igraph

#' @export contact_clustering
#
#  Copyright (C) 2025-2026 Yu Wan <yu.wan@liverpool.ac.uk>, Mohammad Saiful Islam Sajib <saiful.sajib@chrfbd.org>
#  Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
#  Creation: 31 May 2026; the latest update: 31 May 2026

contact_clustering <- function(E = NULL, V = NULL, l_max = Inf) {
    # Validate inputs ###############
    if (is.null(E) || !is.data.frame(E)) {
        stop("E must be a non-NULL data frame of network edges.")
    }
    if (!all(c("Pathway_1", "Pathway_2", "Contact_len") %in% colnames(E))) {
        stop("Edge table E must contain columns 'Pathway_1', 'Pathway_2', and 'Contact_len'.")
    }
    with_V <- !is.null(V)
    if (with_V) {
        if (!is.data.frame(V)) {
            stop("V must be a data frame of network nodes.")
        }
        if (!("Pathway" %in% colnames(V))) {
            stop("Node table V must contain a column named 'Pathway'.")
        }
    }

    # Optionally prune edges by maximum contact length ###############
    # Pruning is skipped when l_max is Inf or non-positive (e.g. 0 or negative).
    if (is.finite(l_max) && l_max > 0) {
        E <- E[E$Contact_len <= l_max, , drop = FALSE]
    }

    # Build the vertex list for the igraph object ###############
    # When V is provided, all its nodes are used as vertices (singletons included).
    # Otherwise, derive the node list solely from the edge table.
    if (!with_V) {
        V <- data.frame(Pathway = sort(unique(c(E$Pathway_1, E$Pathway_2))), stringsAsFactors = FALSE)
    }

    # Build an undirected network as an igraph object ###############
    g <- igraph::graph_from_data_frame(d = E[, c("Pathway_1", "Pathway_2", "Contact_len")],
                                       directed = FALSE,
                                       vertices = V$Pathway)

    # Identify connected components ###############
    comps <- igraph::components(g)  # Since g is an undirected network, the "mode" parameter of the components function has no effect on the result.

    # membership: one row per node, recording its component ID
    membership_df <- data.frame(
        Pathway = names(comps$membership),
        component = as.integer(comps$membership),
        stringsAsFactors = FALSE
    )

    # When V is provided, restore any additional node-level metadata from V
    if (with_V) {
        membership_df <- merge(membership_df, V, by = "Pathway", all.x = TRUE)
    }

    # Preserve the original node ordering
    membership_df <- membership_df[order(membership_df$component, membership_df$Pathway), ]
    rownames(membership_df) <- NULL

    # component_size: one row per component, recording the number of member nodes
    # Component identifiers are sequential integers.
    component_size_df <- data.frame(
        component = as.integer(seq_along(comps$csize)),
        size = as.integer(comps$csize),
        stringsAsFactors = FALSE
    )

    # Assemble and return a Components object ###############
    c <- new(
        "Components",
        V = V,
        E = E,
        membership = membership_df,
        component_size = component_size_df
    )

    return(c)
}
