#' @title Read the pathways spreadsheet
#'
#' @description
#' Read the TSV file recording subjects' pathway identifiers, locations, and the period at each location. Pathways are
#' defined by users. For instance, a pathway can be an admission of a patient. Required columns: (1) Subject: unique
#' subject identifier; (2) Pathway: unique pathway identifier; (3) Location: unique location identifier; (4) Time_start:
#' when the subject appeared at a specific location; (5) Time_end: when the subject left this location.
#'
#' @param movement_table Path to a tab-delimited spreadsheet of user-defined movement pathways.
#' @param count_migrations A logical switch turn on/off counting the number of migrations per pathway. Default: TRUE.
#'
#' @return An object of Class pathways, which consists of a list "pathways" with subject names as indices and a data
#' frame "migrations". The migration tibble is NULL when argument count_migrations = FALSE.
#' @importFrom tibble tibble
#' @importFrom readr read_tsv
#' @importFrom fs file_exists
#' @importFrom tidyselect all_of
#' @importFrom dplyr mutate arrange select group_split group_by summarise
#'
#' @author Yu Wan, \email{yu.wan@liverpool.ac.uk}
#' @author Mohammad Saiful Islam Sajib, \email{saiful.sajib@chrfbd.org}
#'
#' @export read_pathways
#
#  Copyright (C) 2025 Yu Wan <yu.wan@liverpool.ac.uk>, Mohammad Saiful Islam Sajib <saiful.sajib@chrfbd.org>
#  Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
#  Creation: 13 November 2025; the latest update: 6 January 2026

read_pathways <- function(movement_table, count_migrations = TRUE) {
    if (! is.null(movement_table)) {
        if (file_exists(movement_table)) {  # fs::file_exists
            ESSENTIAL_COLUMNS <- c("Subject", "Pathway", "Location", "Time_start", "Time_end")
            DATE_FORMAT <- "%Y-%m-%d"
            movements <- read_tsv(file = movement_table, show_col_types = FALSE, progress = FALSE)
            if (nrow(movements) > 0 & ncol(movements) > 4) {  # The location spreadsheet must not be empty and contain at least five columns.
                if (all(ESSENTIAL_COLUMNS %in% names(movements))) {  # Check if all essential columns are present
                    movements <- movements |>
                        select(all_of(ESSENTIAL_COLUMNS)) |>  # Drop unnecessary columns and fix the order of columns
                        mutate(
                            Time_start = as.Date(Time_start, format = DATE_FORMAT),
                            Time_end = as.Date(Time_end, format = DATE_FORMAT)
                        ) |>
                        arrange(Subject, Pathway, Time_start)  # Sort rows by these three columns in an ascending order
                    pathways <- .build_pathways(movements)
                    if (count_migrations) {
                        migrations <- .summarise_migrations(movements)
                    } else {
                        migrations <- tibble()  # Return an empty tibble to comply with the definition of Class Pathways
                    }
                } else {
                    stop("Error: essential columns Subject, Pathway, Location, Time_start, or Time_end were not found in movement records.")
                }
            } else {
                stop("Error: the TSV file of pathways is empty or has less than five columns.")
            }
        } else {
            stop("Error: the pathway TSV file was not found.")
        }
    } else {
        stop("Error: Parameter 'pathways' is not specified.")
    }
    return(new("Pathways",
               pathways = pathways,
               migrations = migrations))
}

setClass(
    Class = "Pathways",
    slots = list(
        pathways = "list",
        migrations = "data.frame"
    )
)

.build_pathways <- function(movements) {
    # This function assumes rows in movements are sorted by Subject, Pathway ID, and Time_start in
    # an ascending order.
    pathways <- movements |> group_split(Subject)
    names(pathways) <- unique(movements$Subject)  # Named list of tibbles: one tibble per Subject
    pathways <- lapply(pathways, function(df) select(df, -Subject))  # Drop the redundant Subject column from each tibble
    return(pathways)  # To-do: quality check of pathways (e.g., not gaps or conflicts in time while accounting for Admission)
}

.summarise_migrations <- function(movements) {
    # Calculate the number of location transitions - 1 per pathway of each subject
    migration_summary <- movements |>
        group_by(Subject, Pathway) |>
        summarise(
            Migrations = dplyr::n() - 1L,
            Time_start = min(Time_start),
            Time_end = max(Time_end),
            .groups = "drop"
        )
    return(migration_summary)
}

