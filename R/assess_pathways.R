#' @title Quality assessment of input movement data.
#' @description This function evaluates the quality of movement data under four
#' criteria: 1. Time_start <= Time_end at each location; 2. location uniqueness;
#' 3. pathway integrity; 4. pathway uniqueness. See Subsection "Requirements for
#'  the quality of input movement data" on github.com/wanyuac/pathopath for details.
#' @param pathways The named list "pathways" in the output of function read_pathways.
#' @return A logical value indicating whether the data pass (TRUE) or fail (FALSE)
#' this assessment.
#' @author Yu Wan, \email{yu.wan@liverpool.ac.uk}
#' @export
#
#  Copyright (C) 2025-2026 Yu Wan <yu.wan@liverpool.ac.uk>, Mohammad Saiful Islam Sajib <saiful.sajib@chrfbd.org>
#  Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
#  Creation: 9 January 2026; the latest update: 9 January 2026

assess_pathways <- function(pathways) {
    pass <- logical()  # Its length will equal the number of subjects in the pathways parameter.
    subjects <- names(pathways)
    for (s in subjects) {
        pathways_s <- pathways[[s]]  # Pathway table of Subject s
        # Check timestamps of each spatiotemporal tuple
        if (nrow(pathways_s) > 1) {  # Multiple movements of the same pathway or multiple pathways
            qualities_s <- as.logical(mapply(.assess_timestamps,
                                             pathways_s$Time_start, pathways_s$Time_end, pathways_s$Location, pathways_s$Pathway,
                                             MoreArgs = list(subject = s), USE.NAMES = FALSE))
            # Check pathway integrity
            # Note that the location uniqueness is assured if both the .assess_timestamps and .assess_pathway_integrity functions return TRUE.
            pathway_IDs <- sort(unique(pathways_s$Pathway))  # Unique pathway identifiers of the current subject
            if (length(pathway_IDs) > 1) {  # Multiple pathways
                for (p in pathway_IDs) {
                    qualities_s <- append(qualities_s, .assess_pathway_integrity(pathway = subset(pathways_s, Pathway == p), pathway_id = p, subject = s))
                }
            } else {  # Single pathway, multiple locations
                qualities_s <- append(qualities_s, .assess_pathway_integrity(pathway = pathways_s, pathway_id = p, subject = s))
            }
            pass <- append(pass, all(qualities_s))  # Summarises all quality results of subject s and returns one logical value
        } else {  # nrow(pathways_s) == 1; the value of nrow(pathways_s) cannot be <1 as defined by the read_pathways function.
            pass <- append(pass, .assess_timestamps(time_start = pathways_s$Time_start, time_end = pathways_s$Time_end,
                                                    pathway = pathways_s$Pathway, subject = s))  # No other assessment is required under this condition.
        }
    }
    return(all(pass))
}

.assess_timestamps <- function(time_start, time_end, location, pathway, subject) {
    # Check timestamps of a location-time tuple (single row) in a pathway to ensure
    # that time_start is not greater than time_end at each location
    quality <- time_start <= time_end
    if (!quality) {
        warning(paste0("Error: Time_start at location ", location, " in pathway ", pathway, " of subject ", subject, " exceeds Time_end."), call. = FALSE)
    }
    return(quality)
}

.assess_pathway_integrity <- function(pathway, pathway_id, subject) {
    # Check the integrity of a pathway, given that the read_pathways function has sorted
    # locations by Time_start for each pathway (see command "arrange(Subject, Pathway,
    # Time_start))" in function read_pathways.
    num_locations <- nrow(pathway)  # Number of locations in this pathway
    if (num_locations > 1) {
        unbroken_pathway <- TRUE  # Initial (default) value
        loc1_end_time <- pathway$Time_end[1]
        for (i in 2 : num_locations) {
            loc2_start_time <- pathway$Time_start[i]  # Variable type: Date
            time_diff <- loc2_start_time - loc1_end_time  # Pathopath currently only works on dates, so we only consider full days in this function.
            if (time_diff < 0) {  # The subject cannot move to the next place without existing the first place, so this condition indicates a mistake.
                warning(paste0("Error: incorrect timestamps are identified in pathway ", pathway_id, " of subject ", subject, "."), call. = FALSE)
                unbroken_pathway <- FALSE
                break  # There is no need to check other locations in this pathway given the warning message above.
            } else if (time_diff > 1) {  # A time gap exists between two consecutive locations.
                warning(paste0("Error: at least one time gap between consecutive locations are identified in pathway ", pathway_id, " of subject ", subject, "."), call. = FALSE)
                unbroken_pathway <- FALSE
                break  # There is no need to check other locations in this pathway given the warning message above.
            } else if (i < num_locations) {
                loc1_end_time <- pathway$Time_end[i]
            }
        }
    } else {
        unbroken_pathway <- TRUE  # Pathways consisting of single locations do not have any gaps by definition.
    }
    return(unbroken_pathway)
}
