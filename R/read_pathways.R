#' @title Read the pathways spreadsheet
#'
#' @description
#' Read the TSV file recording subjects' pathway identifiers, locations, and the period at each location. Pathways are
#' defined by users. For instance, a pathway can be an admission of a patient. Required columns: (1) Subject: unique
#' subject identifier; (2) Pathway: unique pathway identifier; (3) Location: unique location identifier; (4) Time_start:
#' when the subject appeared at a specific location; (5) Time_end: when the subject left this location.
#'
#' @param pathway_data Path to a spreadsheet of user-defined pathways in the tab-delimited format (TSV).
#' @param count_migrations A logical switch turn on/off counting the number of migrations per pathway. Default: TRUE.
#'
#' @author Yu Wan, \email{yu.wan@liverpool.ac.uk}
#' @author Mohammad Saiful Islam Sajib, \email{saiful.sajib@chrfbd.org}
#'
#' @return An object of Class pathways, which consists of a list "pathways" with subject names as indices and a data
#' frame "migrations". The migration tibble is NULL when argument count_migrations = FALSE.
#' @importFrom dplyr mutate arrange select group_split group_by summarise
#' @export read_pathways
#
#  Copyright (C) 2025 Yu Wan <yu.wan@liverpool.ac.uk>, Mohammad Saiful Islam Sajib <saiful.sajib@chrfbd.org>
#  Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
#  Creation: 13 November 2025; the latest update: 15 November 2025

read_pathways <- function(pathway_data, count_migrations = TRUE) {
    if (! is.null(pathway_data)) {
        ESSENTIAL_COLUMNS <- c("Subject", "Pathway", "Location", "Time_start", "Time_end")
        DATE_FORMAT <- "%Y-%m-%d"
        if (fs::file_exists(pathway_data)) {  # fs::file_exists
            movements <- readr::read_tsv(file = pathway_data, show_col_types = FALSE, progress = FALSE)  # dplyr::read_tsv
            if (nrow(movements) > 0 & ncol(movements) > 4) {  # The location spreadsheet must not be empty and contain at least five columns.
                if (setequal(x = names(movements), y = ESSENTIAL_COLUMNS)) {  # Check if all essential columns are present
                    movements <- movements |>
                        select(all_of(ESSENTIAL_COLUMNS)) |>  # Drop unnecessary columns and fix the order of columns
                        mutate(
                            Time_start = as.Date(movements$Time_start, format = DATE_FORMAT),
                            Time_end = as.Date(movements$Time_end, format = DATE_FORMAT)
                        ) |>
                        arrange(Subject, Pathway, Time_start)  # Sort rows by these three columns in an ascending order
                    pathways <- .build_pathways(movements)
                    if (count_migrations) {
                        migrations <- .summarise_migrations(movements)
                    } else {
                        migrations <- NULL
                    }
                } else {
                    stop("Error: essential columns Subject, Pathway, Location, Time_start, or Time_end were not found in movement records.")
                }
            } else {
                stop("Error: the TSV file of locations is empty or has less than five columns.")
            }
        } else {
            stop("Error: the pathway TSV file was not found.")
        }
    } else {
        stop("Error: Parameter 'pathway_data' is not specified.")
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
    # This function assumes rows in movements are sorted by Subject, Pathway, and Time_start in
    # an ascending order, namely, the output from function import_movements.
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
            Time_start_earliest = min(Time_start),
            Time_end_latest = max(Time_end),
            .groups = "drop"
        )
    return(migration_summary)
}

