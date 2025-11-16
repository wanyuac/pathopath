#' @title Construct a contact network from epidemiological data
#'
#' @description This is the main function of the Pathopath package
#'
#' @param pathway_data Path to a spreadsheet of subjects' movement records arranged in pathways in the tab-delimited
#' format. The Location column stores location accessions, which are unique identifiers of locations at a user-specified level.
#' @param dt Delta t, ± dt days to determine an indirect contact. If both are set, `indirect_cutoff` wins.
#' @param d0 Duration for a point overlap in time. For example, if two subject were only at the same location
#' for one day or even one hour. Default: 1 (day).
#'
#' @author Yu Wan, \email{yu.wan@liverpool.ac.uk}
#' @author Mohammad Saiful Islam Sajib, \email{saiful.sajib@chrfbd.org}
#'
#' @return A named list of two elements "pathways" and "contacts".
#' @export
#
#  Copyright (C) 2025 Yu Wan <yu.wan@liverpool.ac.uk>, Mohammad Saiful Islam Sajib <saiful.sajib@chrfbd.org>
#  Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
#  Creation: 12 November 2025; the latest update: 15 November 2025

pathopath <- function(pathway_data = NULL, d0 = 1, dt = 0) {
    pathways <- read_pathways(pathway_data)  # Parse the input spreadsheet into a named list of pathways and count the number of migrations per pathway
    contacts <- compute_contacts(pathways@pathways, d0, dt)  # Identify and quantify direct and indirect contacts between pathways
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
