#' @title Incorporation of additional patient movements into existing results
#'
#' @description
#' This function adds patient movements into existing results without requiring reprocessing the whole dataset.
#'
#' @param movements Path to a tab-delimited spreadsheet of new movements to be incorporated into the existing results.
#' @param samples Path to a tab-delimited spreadsheet of microbiological/pathological data. Given the complexity of such
#' data, currently users need to provide curated data for all samples, including previously analysed samples.
#' @param previous_results A previous output object of pathopath.
#'
#' @author Yu Wan, \email{yu.wan@liverpool.ac.uk}
#'
#' @importFrom tibble tibble
#' @importFrom dplyr group_by n_distinct summarise intersect bind_rows arrange
#' @importFrom purrr transpose map
#' @export
#
#  Copyright (C) 2025-2026 Yu Wan <yu.wan@liverpool.ac.uk>, Mohammad Saiful Islam Sajib <saiful.sajib@chrfbd.org>
#  Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
#  Creation: 5 January 2026; the latest update: 7 January 2026

add_movements <- function(movements = NULL, samples = NULL, previous_results = NULL) {
    new_movements <- read_pathways(movements)
    # To-do: QC of new data with regards to existing results

    # Identify and quantify contacts within the new movements
    if (length(new_movements@pathways) > 1L) {
        new_contacts_intra <- compute_contacts(pathways = new_movements@pathways,
                                               dt = previous_results@parameters[["dt"]])
    } else {
        new_contacts_intra <- NULL
    }

    # Identify and quantify contacts between new and previously recorded movements
    additional_subjects <- names(new_movements@pathways)  # Note that some subjects may exist in previous_subjects
    previous_subjects <- names(previous_results@pathways)
    subject_pairs <- expand.grid(sa = additional_subjects, sp = previous_subjects) |>
        mutate(across(c(sa, sp), as.character)) |>
        filter(sa != sp)  # Generate pairwise combinations of additional and previous subjects
    contacts_list <- mapply(.pairwise_contacts_across_datasets,
                            subject_pairs$sa, subject_pairs$sp,
                            MoreArgs = list(pathways_add = new_movements@pathways,
                                            pathways_prev = previous_results@pathways,
                                            dt = previous_results@parameters[["dt"]]),
                            SIMPLIFY = FALSE)

    # Combine all new contacts and previously identified ones
    contacts <- bind_rows(new_contacts_intra, bind_rows(contacts_list), previous_results@contacts) |>
        arrange(Subject_1, Subject_2, Location)

    # Re-summarise all contacts, which could be more efficient than meticulously updating individual rows
    # according to updated contacts.
    contact_summary <- contacts |>
        dplyr::filter(Contact != "None") |>
        .summarise_contacts()

    # Update previous results with new information. Here, the whole network is regenerated for simplicity
    # rather than modifying some edges and nodes according to updated contacts, because the network
    # generation is computationally light.
    return(new("Pathopath",
               pathways = .combine_pathways(pathways_add = new_movements@pathways,
                                            pathways_prev = previous_results@pathways),
               migrations = .update_migration_summary(migrations_add = new_movements@migrations,
                                                      migrations_prev = previous_results@migrations),
               contacts = contacts,
               summary = contact_summary,
               network = create_network(contact_summary, samples),
               parameters = list(movements = append(previous_results@parameters[["movements"]], movements),
                                 samples = samples,
                                 dt = previous_results@parameters[["dt"]])))
}

.pairwise_contacts_across_datasets <- function(sa, sp, pathways_add, pathways_prev, dt = 3) {
    # This function supports "add_movements" and is adapted from function ".pairwise_contacts".
    pt_a <- pathways_add[[sa]]  # Pathway tibble of subject sa, which may have >1 pathways.
    pt_p <- pathways_prev[[sp]]
    pathway_accessions_sa <- unique(pt_a$Pathway)  # All pathway accessions (number >= 1) of subject sa
    pathway_accessions_sp <- unique(pt_p$Pathway)  # All pathway accessions (number >= 1) of subject sp
    contacts <- NULL
    for (p_acc_sa in pathway_accessions_sa) {
        for (p_acc_sp in pathway_accessions_sp) {
            new_contacts <- .find_contacts(subject_1 = sa,
                                           subject_2 = sp,
                                           pathway_s1 = filter(pt_a, Pathway == p_acc_sa),
                                           pathway_s2 = filter(pt_p, Pathway == p_acc_sp),
                                           dt = dt)  # Function .find_contacts is defined in script compute_contacts.R.
            if (!is.null(new_contacts)) {
                contacts <- bind_rows(contacts, new_contacts)
            }
        }
    }
    return(contacts)  # NULL is returned if no contact is identified.
}

.combine_pathways <- function(pathways_add, pathways_prev) {
    subjects_prev <- names(pathways_prev)
    for (sa in names(pathways_add)) {
        if (sa %in% subjects_prev) {  # The current subject is already in the previous dataset.
            pathways_prev[[sa]] <- bind_rows(pathways_prev[[sa]], pathways_add[[sa]]) |>
                arrange(Pathway, Time_start)
        } else {  # The current subject is new to the previous dataset.
            pathways_prev[[sa]] <- pathways_add[[sa]]
        }
    }
    return(pathways_prev)  # Return the updated pathway records
}

.update_migration_summary <- function(migrations_add, migrations_prev) {
    subjects_add <- unique(migrations_add$Subject)
    subjects_prev <- unique(migrations_prev$Subject)
    subjects_new <- setdiff(x = subjects_add, y = subjects_prev)  # Unique subjects in the new migration table

    if (length(subjects_new) == length(subjects_add)) {  # The simpliest scenario: All subjects in migrations_add are not present in migrations_prev
        migrations_updated <- bind_rows(migrations_prev, migrations_add)
    } else {  # By definition, length(subjects_new) <= length(subjects_add), so "<" if not "=" here.
        subjects_shared <- setdiff(x = subjects_add, y = subjects_new)  # Subjects shared between the new and previous migration tables. Equivalent to intersect(subjects_add, subjects_prev) but faster to computer.

        # Split previous migration table by shared and unique subjects
        rows_selected <- migrations_prev$Subject %in% subjects_shared
        migrations_to_update <- subset(migrations_prev, rows_selected)  # For shared subjects
        migrations_to_keep <- subset(migrations_prev, !rows_selected)  # Migration summaries of subjects that are only present in migrations_prev

        # Split additional migration table by shared and unique subjects
        rows_selected <- migrations_add$Subject %in% subjects_shared
        migrations_to_merge <- subset(migrations_add, rows_selected)  # To merge with migrations_to_update
        migrations_to_combine <- subset(migrations_add, !rows_selected)  # Migration summaries of subjects that are only present in migrations_add

        # Append unique migrations records (according to subjects) from migrations_add to migrations_to_keep
        migrations_updated <- bind_rows(migrations_to_keep, migrations_to_combine)

        # Update migration summaries for shared subjects
        # This step needs to consider both the subject ID and pathway ID.
        for (s in subjects_shared) {
            migs_sa <- subset(migrations_to_merge, Subject == s)  # By definition, migs_sa and migs_sp each may comprise multiple rows when subjects have multiple pathways.
            migs_sp <- subset(migrations_to_update, Subject == s)  # From migrations_prev
            pathways_a <- migs_sa$Pathway  # Pathway IDs; unique by definition (see function read_pathways), so we do not need to use the unique function
            pathways_p <- migs_sp$Pathway
            pathways_shared <- intersect(x = pathways_a, y = pathways_p)  # Pathway IDs shared between migs_sa and migs_sp
            n <- length(pathways_shared)
            if (n > 0) {  # When migs_sa and migs_sp share some pathways
                if (n < length(pathways_p)) {  # When migs_sp has pathways IDs not present in migs_sa
                    rows_selected <- migs_sp$Pathway %in% pathways_shared
                    migs_sp_shared <- subset(migs_sp, rows_selected)
                    migs_sp_unique <- subset(migs_sp, !rows_selected)
                    migrations_updated <- bind_rows(migrations_updated, migs_sp_unique)
                } else {  # n == length(pathways_p)
                    migs_sp_shared <- migs_sp
                }
                for (p in pathways_a) {
                    add_row <- subset(migs_sa, Pathway == p)  # A single row
                    if (p %in% pathways_shared) {  # Additional movements (n >= 1) of a previously analysed pathway of Subject s
                        prev_row <- subset(migs_sp_shared, Pathway == p)  # This tibble must comprises a single row.
                        new_row <- tibble(Subject = s,
                                          Pathway = p,
                                          Migrations = prev_row$Migrations[1] + add_row$Migrations[1],
                                          Time_start = min(prev_row$Time_start[1], add_row$Time_start[1]),
                                          Time_end = max(prev_row$Time_end[1], add_row$Time_end[1]))
                    } else {  # A new pathway from the additional dataset
                        new_row <- add_row
                    }
                    migrations_updated <- bind_rows(migrations_updated, new_row)
                }
            } else {  # When migs_sa and migs_sp (of Subject s) do not share any pathway
                migrations_updated <- bind_rows(migrations_updated, migs_sp, migs_sa)
            }
        }
        migrations_updated <- arrange(migrations_updated, Subject, Pathway)
    }
    return(migrations_updated)
}
