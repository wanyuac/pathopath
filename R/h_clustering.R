#' @title Complete-linkage hierarchical clustering of subjects by Hamming distances
#'
#' @description This function takes as input a matrix of Hamming distances, which can be a SNP-distance matrix generated
#' from a sequence alignment by software snp-dists, and a vector of distance thresholds. It clusters subjects using
#' complete-linkage hierarchical clustering and then partitions subjects into clusters under each distance threshold. The
#' prefix "h" in the function name stands for "hierarchical".
#'
#' Users do not need to include zero in the threshold vector, since this function always determines clusters of genetically
#' identical subjects, generating column "Cluster_0" in the output tibble of cluster identifiers.
#'
#' @param mat Path to the input distance matrix in the TSV format.
#' @param threshold Positive integers for distance thresholds. Default: c(5, 10, 15, 20, 50).
#' @param excl_subjects An optional vector of subject names for exclusion from the clustering analysis.
#'
#' @return A Clusters object comprising the following slots: (1) data frame "clusters" consisting of columns "Subject",
#' "Cluster_0", ..., "Cluster_n", where n represents the largest distance threshold in the vector "thresholds"; (2) data
#' frame "cluster_counts" of two columns "Threshold", "Cluster_count"; (3) named list "tree", which consists of a dendrogram
#' "tree" (rooted by the clustering algorithm) and a character field "method" (complete linkage); and (4) symmetric matrix
#' "distances", storing the input distance matrix; if the excl_subjects parameter is specified, certain subjects are
#' excluded from this matrix.
#'
#' @importFrom ape as.phylo
#' @importFrom fs file_exists
#' @importFrom dplyr n_distinct
#'
#' @author Yu Wan, \email{yu.wan@liverpool.ac.uk}
#'
#' @export h_clustering
#
#  Copyright (C) 2025-2026 Yu Wan <yu.wan@liverpool.ac.uk>, Mohammad Saiful Islam Sajib <saiful.sajib@chrfbd.org>
#  Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
#  Creation: 9 May 2026; the latest update: 11 May 2026

h_clustering <- function(mat = NULL, thresholds = c(5L, 10L, 15L, 20L, 50L), excl_subjects = "") {
    # require(readr); require(fs); require(dplyr); require(tibble); require(ape)
    # Validation of parameters ###############
    if (is.null(mat)) {
        stop("Error: Parameter 'mat' is not specified.")
    } else if (!file_exists(mat)) {
        stop("Error: the input matrix file was not found.")
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

    # Exclude specific subjects from the matrix
    if (!identical(excl_subjects, "")) {
        indices_to_keep <- !(rownames(dist_mat) %in% excl_subjects)
        if (sum(indices_to_keep) < 2) {
            stop("Error: Less than two subjects remained after excluding specified subjects.")
        } else {
            dist_mat <- dist_mat[indices_to_keep, indices_to_keep, drop = FALSE]  # Excluded specified subjects from the distance matrix
        }
    }

    # Cluster subjects using complete-linkage clustering ###############
    hc <- hclust(d = as.dist(dist_mat), method = "complete")

    # Determine clusters under each distance threshold ###############
    memberships <- cutree(hc, h = thresholds)  # Returns a matrix with the number of rows equals the number of subjects, and columns matching thresholds
    colnames(memberships) <- paste0("Cluster_", as.character(thresholds))
    cls <- as.data.frame(cbind(Subject = rownames(memberships), memberships))
    cls[, -1] <- lapply(cls[, -1], as.integer)  # Restores integer types on columns of Cluster IDs, since the "cbind" function coerces them to character
    rownames(cls) <- NULL  # Removes row names

    # Summarise the number of distinct clusters under each threshold
    cluster_counts <- do.call(rbind, lapply(thresholds, function(t) {
        col <- paste0("Cluster_", as.character(t))
        data.frame(
            Threshold = t,
            Cluster_count = n_distinct(cls[[col]])
        )
    }))
    rownames(cluster_counts) <- NULL

    ## Output ###################################
    # The dendrogram is rooted by the hierarchical clustering algorithm and will be exported as a "phylogenetic" tree.
    return(new("Clusters",
               clusters = cls,
               cluster_counts = cluster_counts,
               tree = list(tree = as.phylo(hc),
                           method = "complete linkage"),
               distances = dist_mat))
}
