#' @title Create undirected, Cytoscape-compatible network files from a contact summary
#'
#' @description This function creates a node table and edge table that can be visualised in Cytoscape or other compatible
#' software for network analysis at the level of movement pathways.
#'
#' @param contact_summary The tibble or data frame of contact summaries in pathopath's output (slot "summary").
#' @param samples Optional path to a tab-delimited spreadsheet of microbiological/pathological data, with three mandatory columns
#' Sample, Subject, and Pathway, followed by optional variable columns such as sample metadata, genotypes, or phenotypes. The
#' sample data will be added to node attributes for network visualisation and analysis.
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
#  Creation: 16 November 2025; the latest update: 7 January 2026

create_network <- function(contact_summary, samples = NULL) {
    # Create an undirected network of contacts between pathways, assuming pathway accessions are unique across
    # the input tibble contact_summary.
    # Create an edge table ###############
    if (nrow(contact_summary) > 0) {
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

        # Incorporate microbiological/pathological data as node attributes ###############
        # This algorithm assumes a single sample per pathway.
        samples <- .read_samples(samples)  # Returns NULL if parameter samples = NULL
        if (! is.null(samples)) {
            samples_unique <- samples |> distinct(Subject, Pathway, .keep_all = TRUE)  # Only select the first row for duplicated combinations of Subject and Pathway
            nodes <- nodes |>
                left_join(samples_unique, by = c("Subject", "Pathway")) |>
                arrange(Pathway, Subject)  # Incorporate sample information as node attributes, while edges are kept unchanged
        }
    } else {
        edges <- tibble()
        nodes <- tibble()
    }
    return(new("Network", V = nodes, E = edges))
}

.read_samples <- function(sample_data = NULL) {
    if (! is.null(sample_data)) {
        if (file_exists(sample_data)) {
            ESSENTIAL_COLUMNS <- c("Sample", "Subject", "Pathway")
            sample_table <- read_tsv(file = sample_data, show_col_types = FALSE, progress = FALSE)
            if (nrow(sample_table) > 0 & ncol(sample_table) >3) {
                if (all(ESSENTIAL_COLUMNS %in% names(sample_table))) {
                    samples <- sample_table[, c(ESSENTIAL_COLUMNS, setdiff(names(sample_table), ESSENTIAL_COLUMNS))]
                } else {
                    samples <- NULL  # Not all essential columns are present in the input spreadsheet
                }
            } else {
                samples <- NULL
            }
        } else {
            samples <- NULL  # The input file does not exist.
        }
    } else {
        samples <- NULL
    }
    return(samples)
}
