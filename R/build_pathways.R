#' @title Build patient pathways from movements
#'
#' @description
#' This function takes as input the tibble "movements" in the output Movements object of import_movements.
#' The tibble has four columns: Admission, Location, Time_start, and Time_end.
#'
#' @param movements Output tibble of function import_movements
#'
#' @return A named list with each element named by a unique subject ID and the element's
#' value is a tibble of this subject's Admission, movements, sorted by Admission and Time_start
#' in an ascending order.
#' @export build_pathways
#
#  Copyright (C) 2025 Yu Wan <yu.wan@liverpool.ac.uk>, Mohammad Saiful Islam Sajib <saiful.sajib@chrfbd.org>
#  Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
#  Creation: 12 November 2025; the latest update: 14 November 2025

build_pathways <- function(movements) {
    # Named list of tibbles: one tibble per Subject
    # Each pathway comprises all movements during an admission.
    # This function assumes rows in movements are sorted by Subject, Admission, and Time_start in
    # an ascending order, namely, the output from function import_movements.
    pathways <- movements %>%
        rename(Pathway = Admission) %>%
        group_split(Subject)

    names(pathways) <- unique(movements$Subject)

    # Drop the Subject column in each tibble
    pathways <- lapply(pathways, function(df) select(df, -Subject))

    return(pathways)  # To-do: quality check of pathways (e.g., not gaps or conflicts in time while accounting for Admission)
}
