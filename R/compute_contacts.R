#' @title Determine direct and indirect contacts from pathways
#'
#' @description
#' A direct contact is defined as an overlap between pathways of distinct subjects at the same location and on the same day.
#' An indirect contact is defined as an overlap between pathways of distinct subjects at the same location and not overlap
#' in time unless allowing ±dt days.
#'
#' @param pathways Named list "pathways" in the output of function read_pathways.
#' @param dt Delta t, ± dt days to determine an indirect contact. Default: 0 (do no detect indirect contacts)
#' @return Edges tibble with from, to, contact_type, direct_days, indirect_days, duration, delta_t, d0. Attribute 'resolution' is set.
#' @importFrom purrr map list_rbind
#' @importFrom dplyr filter bind_rows
#' @importFrom tibble tibble
#' @export compute_contacts
#
#  Copyright (C) 2025 Yu Wan <yu.wan@liverpool.ac.uk>, Mohammad Saiful Islam Sajib <saiful.sajib@chrfbd.org>
#  Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
#  Creation: 13 November 2025; the latest update: 15 November 2025

compute_contacts <- function(pathways, d0 = 1, dt = 0) {
    subject_pairs <- combn(x = names(pathways), m = 2, simplify = FALSE)  # Creates pairwise combinations of subject names as list of two-element character vectors
    contacts <- list_rbind(map(subject_pairs,
                               ~ .pairwise_contacts(.x, pathways = pathways, d0 = d0, dt = dt)))  # Find and quantify contacts between all pathways of any pair of subjects
    return(contacts)
}

.pairwise_contacts <- function(subject_pair, pathways, d0 = 1, dt = 0) {
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
                                           d0 = d0,
                                           dt = dt)
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

.find_contacts <- function(subject_1, subject_2, pathway_s1, pathway_s2, dt = 0) {
    # Find and quantify contacts between two specific pathways
    contacts <- NULL  # Initiate the output tibble for pathways s1 and s2
    loc_p1 <- unique(pathway_s1$Location)  # Locations in the current pathway of subject 1
    loc_p2 <- unique(pathway_s2$Location)  # Locations in the current pathway of subject 2
    loc_overlap <- intersect(x = loc_p1, y = loc_p2)
    if (length(loc_overlap) > 0) {  # No contact occurred when subjects had never been to the same places.
        for (l in loc_overlap) {  # Iterate through shared locations
            periods_p1_s1 <- filter(pathway_s1, Location == l)  # All Time_start and Time_end tuples of subject 1 at this location in the current pathway
            periods_p2_s2 <- filter(pathway_s2, Location == l)  # All Time_start and Time_end tuples of subject 2 at this location in the current pathway
            period_indices_s1 <- 1 : length(periods_p1_s1)  # The length of this vector equals the number of periods of subject 1 at this location.
            period_indices_s2 <- 1 : length(periods_p2_s2)  # The length of this vector equals the number of periods of subject 2 at this location.
            for (i in period_indices_s1) {
                for (j in period_indices_s2) {
                    times_i <- periods_p1_s1[i, c("Time_start", "Time_end")]
                    times_j <- periods_p2_s2[j, c("Time_start", "Time_end")]
                    i_start <- min(times_i)  # Use the mininum and maximum of time stamps to avoid human errors where Time_start > Time_end.
                    i_end <- max(times_i)
                    j_start <- min(times_j)
                    j_end <- max(times_j)
                    period_1 <- new("Periods",
                                    start = i_start,
                                    end = i_end,
                                    length = as.integer(i_end - i_start) + 1)  # Length of period 1, including the start day
                    period_2 <- new("Periods",
                                    start = j_start,
                                    end = j_end,
                                    length = as.integer(j_end - j_start) + 1)  # Length of period 2, including the start day
                    direct_contact <- .find_direct_contact(period_1, period_2)
                    if (direct_contact > 0) {  # A direct contact is found.
                        new_contact <- tibble(Subject_1 = subject_1,
                                              Subject_2 = subject_2,
                                              Location = l,
                                              Contact_length = direct_contact,
                                              Direct_contact = TRUE,
                                              Pathway_1 = periods_p1_s1$Pathway[1],
                                              Pathway_2 = periods_p2_s2$Pathway[1],
                                              Time_start_1 = i_start,
                                              Time_end_1 = i_end,
                                              Period_1 = period_1@length,
                                              Time_start_2 = j_start,
                                              Time_end_2 = j_end,
                                              Period_2 = period_2@length,
                                              dt = NA)  # Direct contacts do not depend on dt.
                        contacts <- bind_rows(contacts, new_contact)
                    } else if (dt > 0) {  # When the user tolerates a maximum of dt days of indirect contact
                        indirect_contact <- .find_indirect_contact(period_1, period_2, dt)
                        if (indirect_contact > 0) {  # An indirect contact is found
                            new_contact <- tibble(Subject_1 = subject_1,
                                                  Subject_2 = subject_2,
                                                  Location = l,
                                                  Contact_length = indirect_contact,
                                                  Direct_contact = FALSE,
                                                  Pathway_1 = periods_p1_s1$Pathway[1],
                                                  Pathway_2 = periods_p2_s2$Pathway[1],
                                                  Time_start_1 = i_start,
                                                  Time_end_1 = i_end,
                                                  Period_1 = period_1@length,
                                                  Time_start_2 = j_start,
                                                  Time_end_2 = j_end,
                                                  Period_2 = period_2@length,
                                                  dt = dt)
                            contacts <- bind_rows(contacts, new_contact)
                        }
                    }
                }
            }
        }
    }
    return(contacts)  # Result "contacts" is either NULL or an non-empty tibble.
}

.find_direct_contact(period_1, period_2) {
    len_max <- as.integer(max(period_1@end, period_2@end) - min(period_1@start, period_2@start)) + 1  # Maximum span of these two periods, including the gap between periods when they do not overlap.
    len_sum <- period_1@length + period_2@length
    if (len_max < len_sum) {  # Detected a direct contact since two periods overlap, including a point contact (e.g., only for one day)
        contact_length <- len_sum - len_max
    } else {
        contact_length <- 0
    }
    return(contact_length)
}

.find_indirect_contact(period_1, period_2, dt = 3) {
    contact_length <- 0
    return(contact_length)
}
