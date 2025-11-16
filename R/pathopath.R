#' @title Construct a contact network from epidemiological data
#'
#' @description This is the main function of the Pathopath package
#'
#' @param pathway_data Path to a spreadsheet of subjects' movement records arranged in pathways in the tab-delimited
#' format. The Location column stores location accessions, which are unique identifiers of locations at a user-specified level.
#' @param dt Delta t, ±dt days (inclusive) to determine an indirect contact. Set it to zero to turn off the detection of
#' indirect contacts. Default: 3.
#'
#' @author Yu Wan, \email{yu.wan@liverpool.ac.uk}
#' @author Mohammad Saiful Islam Sajib, \email{saiful.sajib@chrfbd.org}
#'
#' @return A Pathopath object of three slots: (1) pathways, a list of pathway tibbles named by subject names; (2) migrations,
#' a tibble summarising the number of migrations per pathway; (3) contacts, a tibble of contact status (Direct, Indirect, and None)
#' between pathways of different subjects across shared locations.
#' @export
#
#  Copyright (C) 2025 Yu Wan <yu.wan@liverpool.ac.uk>, Mohammad Saiful Islam Sajib <saiful.sajib@chrfbd.org>
#  Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
#  Creation: 12 November 2025; the latest update: 16 November 2025

pathopath <- function(pathway_data = NULL, dt = 3) {
    pathways <- read_pathways(pathway_data)  # Parse the input spreadsheet into a named list of pathways and count the number of migrations per pathway
    contacts <- compute_contacts(pathways = pathways@pathways, dt = dt)  # Identify and quantify direct and indirect contacts between pathways
    return(new("Pathopath",
               pathways = pathways@pathways,
               migrations = pathways@migrations,
               contacts = contacts))
}


# The output class of function pathpath ##########
setClass(
    Class = "Pathopath",
    slots = list(
        pathways = "list",
        migrations = "data.frame",
        contacts = "data.frame"
    )
)
