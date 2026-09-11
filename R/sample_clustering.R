#' @title Clustering of microbiological samples by Hamming distances
#'
#' @description This function takes as input a matrix of Hamming distances, which can be a SNP-distance matrix generated
#' from a sequence alignment by software snp-dists for pathogens, and a vector of distance thresholds. It clusters samples
#' using either complete-linkage hierarchical clustering or a distance-network-based approach. In the latter, The distance
#' matrix is treated as an equivalent all-to-all bidirectional distance network. Under each distance threshold, this function
#' prunes edges that exceeds the threshold in the network and identifies components in the remaining network as clusters of
#' nodes representing samples. Comparing to the complete-linkage hierarchical clustering, this approach walks through a
#' contact network to identify maximal clusters of nodes based on a given maximum stepwise distance.
#'
#' After clustering, the function partitions samples into clusters under each distance threshold. This function always determines
#' clusters of genetically identical samples (namely, with the threshold of zero), generating column "Cluster_0" in the output
#' data frame of cluster identifiers.
#'
#' @param mat Path to the input distance matrix in the TSV format.
#' @param method Clustering method: "hc", complete-linkage hierarchical clustering (default); "dn", distance-network clustering.
#' @param thresholds Positive integers for distance thresholds. Default: c(0L, 5L, 10L, 15L, 20L, 50L).
#' @param excl_samples An optional vector of sample names to be excluded from the clustering analysis.
#'
#' @return A Clusters object comprising the following slots: (1) data frame "clusters" consisting of columns "Sample",
#' "Cluster_0", ..., "Cluster_n", where n represents the largest distance threshold in the vector "thresholds"; (2) data
#' frame "cluster_counts" of two columns "Threshold", "Cluster_count"; (3) named list "tree", which consists of a dendrogram
#' "tree" (rooted by the clustering algorithm) and a character field "method" (complete linkage) when method = "hc", otherwise,
#' an empty list object; and (4) symmetric matrix "distances", storing the input distance matrix; if the excl_samples parameter is
#' specified, certain samples are excluded from this matrix.
#'
#' @importFrom ape as.phylo
#' @importFrom dplyr n_distinct
#' @importFrom fs file_exists
#' @importFrom igraph graph_from_adjacency_matrix components
#'
#' @author Yu Wan, \email{yu.wan@liverpool.ac.uk}
#'
#' @export sample_clustering
#
#  Copyright (C) 2025-2026 Yu Wan <yu.wan@liverpool.ac.uk>, Mohammad Saiful Islam Sajib <saiful.sajib@chrfbd.org>
#  Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
#  Creation: 11 September 2026; the latest update: 12 September 2026

sample_clustering <- function(mat = NULL, method = "hc", thresholds = c(0L, 5L, 10L, 15L, 20L, 50L), excl_samples = "") {
    # require(readr); require(fs); require(dplyr); require(tibble); require(ape)
    # Validation of parameters ###############
    # Validate the matrix argument
    if (is.null(mat)) {
        stop("Error: Parameter 'mat' is not specified.")
    } else if (!file_exists(mat)) {
        stop("Error: the input matrix file was not found.")
    }

    # Validate the method argument
    if (!method %in% c("hc", "dn")) {
        stop("Error: Parameter 'method' must be either \"hc\" (hierarchical clustering) or \"dn\" (distance-network clustering).")
    }

    # Validate and process the thresholds vector
    if (length(thresholds) > 0) {
        thresholds <- as.integer(thresholds)
        if (any(thresholds < 0)) {
            stop("Error: the vector of distance thresholds contains negative values.")
        }
        if (! 0 %in% thresholds) {
            thresholds <- append(thresholds, 0, after = 0)
        }
        thresholds <- sort(thresholds)  # Sorts thresholds in an ascending order
    } else {
        warning("Since thresholds are not specified, I will take a single threshold zero.")
        thresholds <- 0L
    }

    # Import the distance matrix ###############
    dist_mat <- as.matrix(read.delim(file = mat, check.names = FALSE, row.names = 1, stringsAsFactors = FALSE))

    # Validate the matrix
    if (!isSymmetric(dist_mat)) {
        stop("Error: the input matrix is not symmetric.")
    } else if (any(as.numeric(diag(dist_mat)) != 0)) {
        stop("Error: the diagonal contains at least one non-zero value.")  # This suggests a possible mistake in the generation of this matrix.
    }

    # Exclude specific samples from the matrix
    if (!identical(excl_samples, "")) {
        indices_to_keep <- !(rownames(dist_mat) %in% excl_samples)
        if (sum(indices_to_keep) < 2) {
            stop("Error: Less than two samples remained after excluding specified samples.")
        } else {
            dist_mat <- dist_mat[indices_to_keep, indices_to_keep, drop = FALSE]  # Excluded specified samples from the distance matrix
        }
    }

    # Perform clustering analysis of samples
    if (method == "hc") {
        m <- "Hierarchical clustering"
        clsters <- .sample_hc(dist_mat, thresholds)
    } else {
        m <- "Distance-network clustering"
        clsters <- .sample_dnc(dist_mat, thresholds)
    }

    # Summarise the number of distinct clusters under each threshold
    cluster_counts <- do.call(rbind, lapply(thresholds, function(t) {
        col <- paste0("Cluster_", as.character(t))
        data.frame(
            Threshold = t,
            Cluster_count = n_distinct(clsters$cls[[col]])
        )
    }))
    rownames(cluster_counts) <- NULL

    ## Output ###################################
    # The dendrogram is rooted by the hierarchical clustering algorithm and will be exported as a "phylogenetic" tree.
    return(new("Clusters",
               clusters = clsters$cls,
               cluster_counts = cluster_counts,
               method = m,
               tree = clsters$tr,
               distances = dist_mat))
}


.sample_hc <- function(mat, thresholds) {
    # Cluster samples using complete-linkage clustering ###############
    hc <- hclust(d = as.dist(mat), method = "complete")

    # Determine clusters under each distance threshold ###############
    memberships <- cutree(hc, h = thresholds)  # Returns a matrix with the number of rows equals the number of samples, and columns matching thresholds
    colnames(memberships) <- paste0("Cluster_", as.character(thresholds))
    cls <- as.data.frame(cbind(Sample = rownames(memberships), memberships))
    cls[, -1] <- lapply(cls[, -1], as.integer)  # Restores integer types on columns of Cluster IDs, since the "cbind" function coerces them to character
    rownames(cls) <- NULL  # Removes row names
    return(list(cls = cls, tr = list(tree = as.phylo(hc), method = "complete linkage")))
}


.sample_dnc <- function(mat, thresholds) {
    # Identify components under each distance threshold ###############
    cls <- NULL  # Initiate the result data frame
    for (t in thresholds) {
        # Network pruning
        adj <- mat <= t  # Creates a Boolean adjacency matrix
        adj[is.na(adj)] <- FALSE
        diag(adj) <- FALSE  # Avoids self-connections

        # Component identification
        g <- graph_from_adjacency_matrix(adj, mode = "undirected")
        comp <- components(g)$membership

        # Build this threshold's cluster-membership column
        col_name <- paste0("Cluster_", as.character(t))
        cls_t <- data.frame(Sample = names(comp), Cluster = as.integer(comp), stringsAsFactors = FALSE)
        names(cls_t)[names(cls_t) == "Cluster"] <- col_name
        if (is.null(cls)) {
            cls <- cls_t  # The first and second columns in the result data frame.
        } else {
            cls <- cbind(cls, cls_t[, -1, drop = FALSE])  # Append subsequent columns to the data frame.
        }
    }
    rownames(cls) <- NULL  # Removes row names
    return(list(cls = cls, tr = list()))
}
