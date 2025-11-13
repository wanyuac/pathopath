#' @title Import the movements spreadsheet
#'
#' @description
#' Read the TSV file recording subjects' locations and the period at each location.
#' Required columns: Subject, Location, Time_start, Time_end, Admission.
#'
#' @param movement_tsv Path to a spreadsheet of subjects' location records in the tab-delimited format (TSV).

#' @author Yu Wan, \email{yu.wan@liverpool.ac.uk}
#' @author Mohammad Saiful Islam Sajib, \email{saiful.sajib@chrfbd.org}
#'
#' @return Tibble movements
#' @export
#
#  Copyright (C) 2025 Yu Wan <yu.wan@liverpool.ac.uk>, Mohammad Saiful Islam Sajib <saiful.sajib@chrfbd.org>
#  Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
#  Creation: 13 November 2025; the latest update: 13 November 2025

import_movements <- function(movement_tsv) {
    if (! is.null(movement_tsv)) {
        ESSENTIAL_COLUMNS <- c("Subject", "Location", "Time_start", "Time_end", "Admission")  # The Admission column is not essential for contact detection but it could be used for quality assessment.
        DATE_FORMAT <- "%Y-%m-%d"
        if (file_exists(movement_tsv)) {  # fs::file_exists
            movements <- read_tsv(file = movement_tsv, show_col_types = FALSE, progress = FALSE)  # dplyr::read_tsv
            if (nrow(movements) > 0 & ncol(movements) > 4) {  # The location spreadsheet must not be empty and contain at least five columns.
                if (setequal(x = names(movements), y = ESSENTIAL_COLUMNS)) {  # Check if all essential columns are present
                    movements <- movements[, ESSENTIAL_COLUMNS]  # Drop unnecessary columns and fix the order of columns
                    movements$Time_start <- as.Date(movements$Time_start, format = DATE_FORMAT)
                    movements$Time_end <- as.Date(movements$Time_end, format = DATE_FORMAT)
                } else {
                    stop("Error: essential columns Subject, Location, Time_start, Time_end, or Admission were not found in movement records.")
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

    return(movements)
}

