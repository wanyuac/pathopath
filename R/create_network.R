#' @title Create a undirected Cytoscape-compatible network files from a contact summary
#'
#' @description This function creates a node table and edge table that can be visualised in Cytoscape or other compatible
#' software for network analysis.
#'
#' @param contact_summary The tibble or data frame of contact summaries in pathopath's output (slot "summary").
#' @param genotypes Optional path to a tab-delimited spreadsheet of isolates' genotypical data, with three mandatory column
#' names Subject, Pathway, and Sample. The genotypical data will be added to node attributes for network visualisation and analysis.
#'
#' @return A network object with two tibble slots: V for nodes and E for edges.
#'
#' @author Yu Wan, \email{yu.wan@liverpool.ac.uk}
#' @author Mohammad Saiful Islam Sajib, \email{saiful.sajib@chrfbd.org}
#'
#' @importFrom tibble tibble
#' @importFrom readr read_tsv
#' @importFrom fs file_exists
#' @importFrom tidyr pivot_longer
#' @importFrom tidyselect all_of
#' @importFrom dplyr select mutate rename bind_rows distinct
#' @export create_network
#
#  Copyright (C) 2025 Yu Wan <yu.wan@liverpool.ac.uk>, Mohammad Saiful Islam Sajib <saiful.sajib@chrfbd.org>
#  Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
#  Creation: 16 November 2025; the latest update: 18 November 2025

create_network <- function(contact_summary, genotypes = NULL) {
    # Create an undirected network of contacts between pathways, assuming pathway accessions are unique across the input tibble.
    # To-do: genotypical data can be added by this function as node attributes in the future.
    if (nrow(contact_summary) > 0) {
        # Create an edge table ###############
        if (any(contact_summary$Direct_contact_num > 0)) {
            edges_direct <- contact_summary |>
                dplyr::filter(Direct_contact_num > 0) |>
                select(Pathway_1, Pathway_2, Direct_contact_len, Subject_1, Subject_2,
                       Direct_contact_num, Direct_contact_loc) |>
                mutate(Contact_type = "Direct", .before = Subject_1) |>
                rename(Contact_len = Direct_contact_len,
                       Contact_num = Direct_contact_num,
                       Location_num = Direct_contact_loc)
        } else {
            edges_direct <- tibble()
        }
        if (any(contact_summary$Indirect_contact_num > 0)) {
            edges_indirect <- contact_summary |>
                dplyr::filter(Indirect_contact_num > 0) |>
                select(Pathway_1, Pathway_2, Indirect_contact_len, Subject_1, Subject_2,
                       Indirect_contact_num, Indirect_contact_loc) |>
                mutate(Contact_type = "Indirect", .before = Subject_1) |>
                rename(Contact_len = Indirect_contact_len,
                       Contact_num = Indirect_contact_num,
                       Location_num = Indirect_contact_loc)
        } else {
            edges_indirect <- tibble()
        }
        edges <- bind_rows(edges_direct, edges_indirect)

        # Create a node table by reshaping the edge table ###############
        nodes <- edges |>
            select(Pathway_1, Subject_1, Pathway_2, Subject_2) |>
            pivot_longer(cols = c(Pathway_1, Pathway_2, Subject_1, Subject_2),
                         names_to = c(".value", "pair"),
                         names_pattern = "(Pathway|Subject)_(.)") |>
            distinct(Pathway, Subject)

        # Incorporate user's optional genotype data ###############
        genotypes <- .read_genotypes(genotypes)
        if (! is.null(genotypes)) {
            genotypes_unique <- genotypes |> distinct(Subject, Pathway, .keep_all = TRUE)  # Only select the first row for duplicated combinations of Subject and Pathway
            nodes <- nodes |>
                left_join(genotypes_unique, by = c("Subject", "Pathway")) |>
                arrange(Pathway, Subject)
        }
    } else {
        edges <- tibble()
        nodes <- tibble()
    }
    return(new("Network", V = nodes, E = edges))
}

.read_genotypes <- function(genotype_data = NULL) {
    if (! is.null(genotype_data)) {
        if (file_exists(genotype_data)) {
            ESSENTIAL_COLUMNS <- c("Sample", "Subject", "Pathway")
            genotype_table <- read_tsv(file = genotype_data, show_col_types = FALSE, progress = FALSE)
            if (nrow(genotype_table) > 0 & ncol(genotype_table) >3) {
                if (all(ESSENTIAL_COLUMNS %in% names(genotype_table))) {
                    genotypes <- genotype_table[, c(ESSENTIAL_COLUMNS, setdiff(names(genotype_table), ESSENTIAL_COLUMNS))]
                } else {
                    genotypes <- NULL  # Not all essential columns are present in the input spreadsheet
                }
            } else {
                genotypes <- NULL
            }
        } else {
            genotypes <- NULL  # The input file does not exist.
        }
    } else {
        genotypes <- NULL
    }
    return(genotypes)
}
