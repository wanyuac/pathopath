#' @title Summarie migration history
#'
#' @description
#' Summary tibble of migrations: counts and time span per Subject
#'
#' @param pathways Output list from function build_pathways
#'
#' @return A summary tibble of four columns: Subject, Migration_count, Time_in_earlies, Time_out_latest.
#' @export
#
#  Copyright (C) 2025 Yu Wan <yu.wan@liverpool.ac.uk>, Mohammad Saiful Islam Sajib <saiful.sajib@chrfbd.org>
#  Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
#  Creation: 12 November 2025; the latest update: 13 November 2025

summarise_migrations <- function(pathways) {
    movements <- bind_rows(pathways, .id = "Subject")  # Combine list of tibbles; .id uses the list names as Subject

    migration_summary <- movements %>%
        group_by(Subject) %>%
        summarise(
            Migrations = n() - 1L,
            Time_start_earliest = min(Time_start),
            Time_end_latest = max(Time_end),
            .groups = "drop"
        )

    return(migration_summary)
}
