#' @title Determine direct and indirect contacts from pathways
#'
#' @description
#' A direct contact is defined as an overlap between pathways of distinct subjects at the same location and on the same day.
#' An indirect contact is defined as an overlap between pathways of distinct subjects at the same location and not overlap
#' in time unless allowing ±dt days.
#'
#' @param pathways Named list "pathways" in the output of function read_pathways.
#' @param dt Delta t, ± dt days to determine an indirect contact. Set it to zero to turn off the detection of indirect contacts.
#' Default: 3.
#' @return Edges tibble with from, to, contact_type, direct_days, indirect_days, duration, delta_t, d0. Attribute 'resolution' is set.
#' @importFrom tibble tibble
#' @importFrom purrr map list_rbind
#' @importFrom dplyr filter bind_rows across mutate
#' @export compute_contacts
#
#  Copyright (C) 2025 Yu Wan <yu.wan@liverpool.ac.uk>, Mohammad Saiful Islam Sajib <saiful.sajib@chrfbd.org>
#  Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
#  Creation: 13 November 2025; the latest update: 16 November 2025

compute_contacts <- function(pathways, dt = 3) {
    # This function creates pairwise combinations of subject names as list of two-element character vectors using the
    # combn function, and then applies the .pairwise_contacts function to detect contacts between each pair of subjects.
    if (dt < 0) {
        stop("Error: dt cannot be negative.")
    }
    if (length(pathways) >= 2L) {
        subject_pairs <- combn(x = names(pathways), m = 2, simplify = FALSE)
        contacts <- list_rbind(map(subject_pairs, ~ .pairwise_contacts(.x, pathways = pathways, dt = dt)))
    } else {
        stop("Input error: the list of pathways has less than two elements.")
    }
    return(contacts)
}

.pairwise_contacts <- function(subject_pair, pathways, dt) {
    # Find and quantify contacts between all pairs of pathways of two specific subjects s1 and s2
    s1 <- subject_pair[[1]]  # Subject 1
    s2 <- subject_pair[[2]]  # Subject 2
    pt1 <- pathways[[s1]]  # Pathway tibble of subject 1, which may have >1 pathways.
    pt2 <- pathways[[s2]]  # Pathway tibble of subject 2, which may have >1 pathways.
    pathway_accessions_s1 <- unique(pt1$Pathway)  # All pathway accessions (number >= 1) of subject 1
    pathway_accessions_s2 <- unique(pt2$Pathway)  # All pathway accessions (number >= 1) of subject 2
    contacts_s1_s2 <- NULL  # Initiate the output tibble for subjects 1 and 2
    for (p_acc_s1 in pathway_accessions_s1) {
        for (p_acc_s2 in pathway_accessions_s2) {
            new_contacts <- .find_contacts(subject_1 = s1,
                                           subject_2 = s2,
                                           pathway_s1 = filter(pt1, Pathway == p_acc_s1),
                                           pathway_s2 = filter(pt2, Pathway == p_acc_s2),
                                           dt = dt)  # This function returns a NULL value when two pathways do not cross over in space.
            if (!is.null(new_contacts)) {
                contacts_s1_s2 <- bind_rows(contacts_s1_s2, new_contacts)
            }
        }
    }
    return(contacts_s1_s2)
}

setClass(
    Class = "Periods",
    slots = list(
        start = "Date",
        end = "Date",
        length = "integer"
    )
)

.find_contacts <- function(subject_1, subject_2, pathway_s1, pathway_s2, dt) {
    # Find and quantify contacts across shared locations between two specific pathways
    contacts <- NULL  # Initiate the output tibble for pathways s1 and s2
    loc_p1 <- unique(pathway_s1$Location)  # Locations in the current pathway of subject 1
    loc_p2 <- unique(pathway_s2$Location)  # Locations in the current pathway of subject 2
    loc_overlap <- intersect(x = loc_p1, y = loc_p2)
    if (length(loc_overlap) > 0) {  # No contact occurred when subjects had never been to the same places, returning a NULL value from this function.
        for (l in loc_overlap) {  # Iterate through shared locations and report contact status at each location
            periods_p1_s1 <- filter(pathway_s1, Location == l)  # All Time_start and Time_end tuples of subject 1 at this location in the current pathway
            periods_p2_s2 <- filter(pathway_s2, Location == l)  # All Time_start and Time_end tuples of subject 2 at this location in the current pathway
            period_indices_s1 <- seq_len(nrow(periods_p1_s1))  # The length of this vector equals the number of periods of subject 1 at this location.
            period_indices_s2 <- seq_len(nrow(periods_p2_s2))  # The length of this vector equals the number of periods of subject 2 at this location.
            for (i in period_indices_s1) {
                for (j in period_indices_s2) {
                    times_i <- periods_p1_s1[i, c("Time_start", "Time_end")]
                    times_j <- periods_p2_s2[j, c("Time_start", "Time_end")]
                    i_start <- min(times_i$Time_start[1], times_i$Time_end[1])  # Use the mininum and maximum of time stamps to avoid human errors where Time_start > Time_end.
                    i_end <- max(times_i$Time_start[1], times_i$Time_end[1])
                    j_start <- min(times_j$Time_start[1], times_j$Time_end[1])  # Same attempt to fix human errors
                    j_end <- max(times_j$Time_start[1], times_j$Time_end[1])
                    period_1 <- new("Periods",
                                    start = i_start,
                                    end = i_end,
                                    length = as.integer(i_end - i_start) + 1L)  # Length of period 1, including the start day; note that "+1" only applies when time is measured by days.
                    period_2 <- new("Periods",
                                    start = j_start,
                                    end = j_end,
                                    length = as.integer(j_end - j_start) + 1L)  # Length of period 2, including the start day; note that "+1" only applies when time is measured by days.
                    new_contact <- .identify_contact(period_1, period_2, dt)  # Results: direct, indirect, or no contact
                    new_row <- tibble(Subject_1 = subject_1,
                                      Subject_2 = subject_2,
                                      Location = l,
                                      Contact = new_contact$status,
                                      Length = new_contact$length,
                                      Pathway_1 = periods_p1_s1$Pathway[1],
                                      Pathway_2 = periods_p2_s2$Pathway[1],
                                      Time_start_1 = i_start,
                                      Time_end_1 = i_end,
                                      Time_start_2 = j_start,
                                      Time_end_2 = j_end,
                                      Period_1 = period_1@length,
                                      Period_2 = period_2@length,
                                      Interval = new_contact$interval,
                                      dt = dt)
                    contacts <- bind_rows(contacts, new_row)  # Add a new row regardless of whether a direct or indirect contact is detected.
                }
            }
        }
    }
    return(contacts)  # Result "contacts" is either NULL or an non-empty tibble.
}

.identify_contact <- function(period_1, period_2, dt) {
    L_max <- as.integer(max(period_1@end, period_2@end) - min(period_1@start, period_2@start)) + 1  # Maximum end-to-end span of these two periods, including the gap between periods when they do not overlap.
    L12 <- period_1@length + period_2@length  # L1 + L2
    interval <- L_max - L12
    if (interval < 0) {
        # L_max < L12. Detected a direct contact since two periods overlap, including a point contact (e.g., only for one day)
        # Two patients were at the same location for just one day when interval = -1.
        contact_status <- "Direct"  # Note that direct contacts do not depend on the dt parameter.
        contact_length <- -interval  # The length of this direct contact is measured by days.
        interval <- 0  # Reset the interval to zero for direct contacts to make it more intuitive than negative values.
    } else if (interval < dt) {
        # An indirect contact is identified. Since interval >=0 here, no indirect contact will be detected when dt = 0.
        # "interval < dt" is equivalent to "interval < dt - 1" when time is measured by days. For instance, let dt = 3
        # and an interval of zero, which means patient 2 entered a specific location on the next day of patient 1's
        # depature from this location. Similarly, patient 2 entered the location on the third day after patient 1's
        # departure when interval = 2, so an indirect contact should be reported between these two patients when dt = 3.
        contact_status <- "Indirect"
        # Lengths of indirect contacts depend on the dt parameter.
        # An example of the calculation of indirect-contact lengths: patient 1 left location A on 2025-11-01 and patient 2
        # arrived at this location on 2025-11-03, then the interval or gap in time is one day. Given dt = 3, which virtually
        # extends patient 1's presence to 2025-11-04, then the length of this indirect contact is two days since
        # patient 2's arrival.
        contact_length <- dt - interval
    } else {
        # No attempt to detect indirect contacts when dt = 0. So users need to consider their dt values when interpreting "None".
        contact_status <- "None"
        contact_length <- 0
    }
    return(list(status = contact_status,
                length = contact_length,
                interval = interval))
}
