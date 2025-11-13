#' @title Construct a contact network from epidemiological data
#'
#' @description This is the main function of the Pathopath package
#'
#' @param movement_data Path to a spreadsheet of subjects' movement records in the tab-delimited format. The Location
#' column stores location accessions, which are unique identifiers of locations at a user-specified level.
#'
#' @author Yu Wan, \email{yu.wan@liverpool.ac.uk}
#' @author Mohammad Saiful Islam Sajib, \email{saiful.sajib@chrfbd.org}
#'
#' @return An object of Class ContactNetwork
#'
#' @export pathopath
#
#  Copyright (C) 2025 Yu Wan <yu.wan@liverpool.ac.uk>, Mohammad Saiful Islam Sajib <saiful.sajib@chrfbd.org>
#  Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
#  Creation: 12 November 2025; the latest update: 13 November 2025

pathopath <- function(movement_data = NULL) {
    movements <- import_movements(movement_data)
    pathways <- build_pathways(movements)
    migrations <- summarise_migrations(pathways)  # Summary of migrations (namely, number of location transitions - 1) per subject
    contacts <- compute_contacts(pathways)
    contact_network <- new("Pathopath",
                           movements = movements,
                           pathways = pathways,
                           migrations = migrations)
    return(contact_network)
}

# Define the output class ###############
setClass(
    "Pathopath",
    slots = list(
        movements = "data.frame",
        pathways = "list",
        migrations = "data.frame"
    )
)
