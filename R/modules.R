# Modules of the main function pathopath.
# Copyright (C) 2025 Yu Wan <wanyuac@gmail.com> Mohammad Saiful Islam Sajib <saiful.sajib@chrfbd.org>
# Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
# Creation: 12 November 2025; the latest update: 12 November 2025

.import_epi_data <- function(epi_tsv) {
    require(readr)
    epi_data <- read_tsv(file = epi_tsv, show_col_types = FALSE, progress = FALSE)
    if (nrow(epi_data) == 0) {
        epi_data <- NULL
    }
    return(epi_data)
}
