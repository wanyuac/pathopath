#' @title Construct a contact network from epidemiological data
#'
#' @description This is the main function of the Pathopath package
#'
#' @param epi_tsv Path to a spreadsheet of epidemiological data in the tab-delimited format (TSV)
#'
#' @author Yu Wan, \email{yu.wan@liverpool.ac.uk}
#' @author Mohammad Saiful Islam Sajib, \email{saiful.sajib@chrfbd.org}
#'
#' @return An object of Class ContactNetwork
#'
#' @export

pathopath <- function(epi_tsv = NULL) {
    if (! is.null(epi_tsv)) {
        epi_data <- .import_epi_data(epi_tsv)
    }
    return(epi_data)
}
