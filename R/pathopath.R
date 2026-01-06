#' @title Construct a contact network from epidemiological data
#'
#' @description This is the main function of the Pathopath package
#'
#' @param movement_table Path to a tab-delimited spreadsheet of subjects' movement records arranged in pathways. The Location
#' column stores location accessions, which are unique identifiers of locations at a user-specified level. The pathway accessions
#' must be unique across the input data. Note that by definition, subject and pathway accessions are in
#' one-to-one mapping.
#' @param genotype_table Optional path to a tab-delimited spreadsheet of isolates' genotypical data, with three mandatory column
#' names Subject, Pathway, and Sample. The genotypical data will be added to node attributes for network visualisation and analysis.
#' @param dt Delta t, ±dt days (inclusive) to determine an indirect contact. Set it to zero to turn off the detection of
#' indirect contacts. Default: 3.
#'
#' @return A Pathopath object of six slots: (1) pathways, a list of pathway tibbles named by subject names; (2) migrations,
#' a tibble summarising the number of migrations per pathway; (3) contacts, a tibble of contact status (Direct, Indirect, and None)
#' between pathways of different subjects across shared locations; (4) summary, a tibble summarising direct and indirect
#' contacts for each unique combination of Subject_1, Subject_2, Pathway_1, and Pathway_2; (5) network, a Network object
#' with two slots - V for nodes and E for edges - for exportation as node and edge tables compatible with Cytoscape; (6) parameters,
#' a named list storing arguments of the pathopath function for reproducibility and recalculation for contacts.
#'
#' @author Yu Wan, \email{yu.wan@liverpool.ac.uk}
#' @author Mohammad Saiful Islam Sajib, \email{saiful.sajib@chrfbd.org}
#'
#' @importFrom tibble tibble
#' @importFrom dplyr group_by n_distinct summarise
#' @export
#
#  Copyright (C) 2025 Yu Wan <yu.wan@liverpool.ac.uk>, Mohammad Saiful Islam Sajib <saiful.sajib@chrfbd.org>
#  Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
#  Creation: 12 November 2025; the latest update: 6 January 2026

pathopath <- function(movement_table = NULL, genotype_table = NULL, dt = 3) {
    # Parse the input spreadsheet into a named list of pathways and count the number of migrations per pathway
    pathways <- read_pathways(movement_table)

    # Identify and quantify direct and indirect contacts between pathways. The results include absence of contacts
    # between pathways because indirect contacts depend on the dt parameter.
    contacts <- compute_contacts(pathways = pathways@pathways, dt = dt)

    # Summarise contacts for unique combinations of Subject_1, Subject_2, Pathway_1, and Pathway_2.
    # Function .summarise_contacts returns an empty tibble when no contact is present in the "contacts" tibble.
    contact_summary <- contacts |>
        dplyr::filter(Contact != "None") |>
        .summarise_contacts()

    # Create the output object
    return(new("Pathopath",
               pathways = pathways@pathways,
               migrations = pathways@migrations,
               contacts = contacts,
               summary = contact_summary,
               network = create_network(contact_summary, genotype_table),
               parameters = list(movement_table = movement_table,
                                 genotype_table = genotype_table,
                                 dt = dt)))
}

# Functional modules ##########
setClass(Class = "Network",
         slots = list(V = "data.frame",
                      E = "data.frame"))

setClass(
    # The output class of function pathpath
    Class = "Pathopath",
    slots = list(
        pathways = "list",
        migrations = "data.frame",
        contacts = "data.frame",
        summary = "data.frame",
        network = "Network",
        parameters = "list"))

.summarise_contacts <- function(contacts) {
    # This function summarises contacts for each combination of Subject_1, Subject_2, Pathway_1, and Pathway_2.
    # For each combination, it calculates the total lengths of direct and indirect contacts, respectively, according to
    # the Length column in tibble contacts, and counts the number of unique locations where direct and indirect contacts
    # occurred. It also counts the number of direct and indirect contacts, respectively, for each combination.
    # Columns in the output tibble:
    #  *_num: number of direct/indirect contacts
    #  *_len: total length of direct/indirect contacts
    #  *_loc: number of locations where direct/indirect contacts took place
    if (nrow(contacts) > 0) {
        contact_summary <- contacts |>
            group_by(Subject_1, Subject_2, Pathway_1, Pathway_2) |>
            summarise(Direct_contact_num = sum(Contact == "Direct"),
                      Direct_contact_len = sum(Length[Contact == "Direct"], na.rm = TRUE),
                      Direct_contact_loc = n_distinct(Location[Contact == "Direct"]),
                      Indirect_contact_num = sum(Contact == "Indirect"),
                      Indirect_contact_len = sum(Length[Contact == "Indirect"], na.rm = TRUE),
                      Indirect_contact_loc = n_distinct(Location[Contact == "Indirect"]),
                      .groups = "drop")
    } else {
        contact_summary <- tibble()
    }
    return(contact_summary)
}
