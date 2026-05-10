#' @title Clustering subjects by pairwise genetic distances
#' @description This function takes as input a matrix of genetic distances, which can be a SNP-distance matrix (Hamming
#' distances) generated from a sequence alignment by software snp-dists, and clusters subjects based on a distance threshold.
#'
#' @param matrix_tsv Path to the input distance matrix in the TSV format.
#' @param threshold Threshold of SNP distances. Integer. Default: 15.
#' @param excl_subjects An optional vector of subject names for exclusion from the clustering analysis.
#'
#' @return A tibble comprising three columns: Subject, Threshold, and Genetic_cluster
#'
#' @import readr
#' @importFrom fs file_exists
#' @importFrom dplyr group_by n_distinct summarise filter
#' @importFrom tibble tibble
#' @importFrom ape as.phylo
#'
#' @author Yu Wan, \email{yu.wan@liverpool.ac.uk}
#'
#' @export clustering_by_genetic_dists
#
#  Copyright (C) 2025-2026 Yu Wan <yu.wan@liverpool.ac.uk>, Mohammad Saiful Islam Sajib <saiful.sajib@chrfbd.org>
#  Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
#  Creation: 9 May 2026; the latest update: 10 May 2026

clustering_by_genetic_dists <- function(matrix_tsv = NULL, threshold = 15, excl_subjects = "") {
    # require(readr); require(fs); require(dplyr); require(tibble); require(ape)
    # Sanity check ###############
    if (is.null(matrix_tsv)) {
        stop("Error: Parameter 'matrix_tsv' is not specified.")
    } else if (!file_exists(matrix_tsv)) {
        stop("Error: the input matrix file was not found.")
    }

    # Import the distance matrix ###############
    dist_mat <- as.matrix(read.delim(file = matrix_tsv, check.names = FALSE, row.names = 1, stringsAsFactors = FALSE))
    if (!isSymmetric(dist_mat)) {
        stop("Error: the input matrix is not symmetric.")
    } else if (any(as.numeric(diag(dist_mat)) != 0)) {
        stop("Error: the diagonal contains at least one non-zero value.")  # This suggests a prossible mistake in the generation of this matrix.
    }

    # Exclude specific subjects ###############
    if (excl_subjects != "") {
        indices_to_keep <- !(rownames(dist_mat) %in% excl_subjects)
        if (sum(indices_to_keep) < 2) {
            stop("Error: Less than two subjects remained after exclusing specified subjects.")
        } else {
            dist_mat <- dist_mat[indices_to_keep, indices_to_keep, drop = FALSE]  # Excluded specified subjects from the distance matrix
        }
    }

    # Cluster subjects using complete-linkage clustering ###############
    hc <- hclust(d = as.dist(dist_mat), method = "complete")

    # Cut at the SNP threshold
    threshold <- as.integer(threshold)
    membership <- cutree(hc, h = threshold)  # Cut height = threshold; returns a vector of integers named with isolate names
    cls <- tibble(Isolate = names(membership), Cluster = as.integer(membership))

    ## Output ###################################
    results <- list(
        clusters = cls,
        cluster_num = max(cls$Cluster),
        hclust_output = hc,
        dendrogram = as.phylo(hc),
        dist_mat = dist_mat,
        input = matrix_tsv,
        threshold = threshold,
        excl_subjects = excl_subjects
    )

    return(results)
}
