#' @title Import the movements spreadsheet
#'
#' @description
#' Read the TSV file recording subjects' locations and the period at each location.
#' Required columns: (1) Subject: unique subject identifier; (2) Admission: unique admission identifier;
#' (3) Location: unique location identifier; (4) Time_start: when the subject appeared at a specific location;
#' (5) Time_end: when the subject left this location.
#'
#' Each admission defines a movement pathway.
#'
#' @param movement_tsv Path to a spreadsheet of subjects' location records in the tab-delimited format (TSV).
#' @param count_migrations A logical switch turn on/off counting the number of migrations per admission
#' (pathway). Default: TRUE.
#'
#' @author Yu Wan, \email{yu.wan@liverpool.ac.uk}
#' @author Mohammad Saiful Islam Sajib, \email{saiful.sajib@chrfbd.org}
#'
#' @return An object of Class Movements, which consists of two tibble slots: movements and migrations. The
#' migration tibble is NULL when argument count_migrations = FALSE.
#' @export import_movements
#
#  Copyright (C) 2025 Yu Wan <yu.wan@liverpool.ac.uk>, Mohammad Saiful Islam Sajib <saiful.sajib@chrfbd.org>
#  Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
#  Creation: 13 November 2025; the latest update: 14 November 2025

import_movements <- function(movement_tsv, count_migrations = TRUE) {
    if (! is.null(movement_tsv)) {
        ESSENTIAL_COLUMNS <- c("Subject", "Admission", "Location", "Time_start", "Time_end")  # The Admission column is not essential for contact detection but it could be used for quality assessment.
        DATE_FORMAT <- "%Y-%m-%d"
        if (file_exists(movement_tsv)) {  # fs::file_exists
            movements <- read_tsv(file = movement_tsv, show_col_types = FALSE, progress = FALSE)  # dplyr::read_tsv
            if (nrow(movements) > 0 & ncol(movements) > 4) {  # The location spreadsheet must not be empty and contain at least five columns.
                if (setequal(x = names(movements), y = ESSENTIAL_COLUMNS)) {  # Check if all essential columns are present
                    movements <- movements %>%
                        dplyr::select(all_of(ESSENTIAL_COLUMNS)) %>%  # Drop unnecessary columns and fix the order of columns
                        mutate(
                            Time_start = as.Date(movements$Time_start, format = DATE_FORMAT),
                            Time_end = as.Date(movements$Time_end, format = DATE_FORMAT)
                        ) %>%
                        arrange(Subject, Admission, Time_start)  # Sort rows by these three columns in an ascending order
                    if (count_migrations) {
                        migrations <- .summarise_migrations(movements)
                    } else {
                        migrations <- NULL
                    }
                } else {
                    stop("Error: essential columns Subject, Admission, Location, Time_start, or Time_end were not found in movement records.")
                }
            } else {
                stop("Error: the TSV file of locations is empty or has less than five columns.")
            }
        } else {
            stop("Error: the movement TSV file was not found.")
        }
    } else {
        stop("Error: Parameter 'movement_tsv' is not specified.")
    }

    return(new("Movements",
               movements = movements,
               migrations = migrations))
}

# Functional modules ###############
.summarise_migrations <- function(movements) {
    # Calculate the number of location transitions - 1 per pathway (admission) of each subject
    migration_summary <- movements %>%
        group_by(Subject, Admission) %>%
        summarise(
            Migrations = n() - 1L,
            Time_start_earliest = min(Time_start),
            Time_end_latest = max(Time_end),
            .groups = "drop"
        ) %>%
        rename(Pathway = Admission)  # Rename the column for convenience of data linkage between results of the pathopath pipeline

    return(migration_summary)
}

setClass(
    "Movements",
    slots = list(
        movements = "data.frame",
        migrations = "data.frame"
    )
)
