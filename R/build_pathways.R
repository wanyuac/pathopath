#' @title Build patient pathways from movements
#'
#' @description
#' This function takes as input the output tibble of import_movements. The tibble
#' has four columns: Subject, Location, Time_start, and Time_end.
#'
#' @param movements Output tibble of function import_movements
#'
#' @return A named list with each element named by a unique subject ID and the element's
#' value is a tibble of this subject's movements, sorted by Time_start from the earliest
#' to the latest.
#' @export build_pathways
#
#  Copyright (C) 2025 Yu Wan <yu.wan@liverpool.ac.uk>, Mohammad Saiful Islam Sajib <saiful.sajib@chrfbd.org>
#  Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
#  Creation: 12 November 2025; the latest update: 13 November 2025

build_pathways <- function(movements) {
    if (! "tbl" %in% class(movements)) {
        movements <- as_tibble(movements) %>%
            arrange(Subject, Time_start)  # Sort rows by Subject and Time_start in an ascending order
    }

    # Named list of tibbles: one tibble per Subject
    pathways <- movements %>%
        group_split(Subject)

    names(pathways) <- movements %>%
        distinct(Subject) %>%
        pull(Subject)  # Equivalent to movements$Subject after the distinct function.

    # Drop the Subject column in each tibble
    pathways <- lapply(pathways, function(df) select(df, -Subject))

    return(pathways)  # To-do: quality check of pathways (e.g., not gaps or conflicts in time while accounting for Admission)
}
