#' @title Create a undirected Cytoscape-compatible network files from a contact summary
#'
#' @description This function creates a node table and edge table that can be visualised in Cytoscape or other compatible
#' software for network analysis.
#'
#' @param contact_summary The tibble or data frame of contact summaries in pathopath's output (slot "summary").
#'
#' @return A network object with two tibble slots: V for nodes and E for edges.
#'
#' @author Yu Wan, \email{yu.wan@liverpool.ac.uk}
#' @author Mohammad Saiful Islam Sajib, \email{saiful.sajib@chrfbd.org}
#'
#' @importFrom tibble tibble
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr select mutate rename bind_rows distinct
#' @export create_network
#
#  Copyright (C) 2025 Yu Wan <yu.wan@liverpool.ac.uk>, Mohammad Saiful Islam Sajib <saiful.sajib@chrfbd.org>
#  Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
#  Creation: 16 November 2025; the latest update: 16 November 2025

create_network <- function(contact_summary) {
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
    } else {
        edges <- tibble()
        nodes <- tibble()
    }
    return(new("Network", V = nodes, E = edges))
}
