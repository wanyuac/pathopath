#' @title Construct a contact network from epidemiological data
#'
#' @description This is the main function of the Pathopath package
#'
#' @param pathway_data Path to a spreadsheet of subjects' movement records arranged in pathways in the tab-delimited
#' format. The Location column stores location accessions, which are unique identifiers of locations at a user-specified level.
#'
#' @author Yu Wan, \email{yu.wan@liverpool.ac.uk}
#' @author Mohammad Saiful Islam Sajib, \email{saiful.sajib@chrfbd.org}
#'
#' @return A named list of two elements "pathways" and "contacts".
#' @export
#
#  Copyright (C) 2025 Yu Wan <yu.wan@liverpool.ac.uk>, Mohammad Saiful Islam Sajib <saiful.sajib@chrfbd.org>
#  Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
#  Creation: 12 November 2025; the latest update: 14 November 2025

pathopath <- function(pathway_data = NULL) {
    pathways <- read_pathways(pathway_data)
    contacts <- compute_contacts(pathways)
    return(new("Pathopath",
               pathways = pathways@pathways,
               migrations = pathways@migrations,
               contacts = contacts))
}


# The output class of function pathpath ##########
setClass(
    "Pathopath",
    slots = list(
        pathways = "list",
        migrations = "data.frame",
        contacts = "data.frame"
    )
)
