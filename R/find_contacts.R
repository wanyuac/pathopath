#' @title Determine direct and indirect contacts from pathways
#'
#' @description
#' A direct contact is defined as an overlap between pathways of distinct subjects at the same location and on the same day.
#' An indirect contact is defined as an overlap between pathways of distinct subjects at the same location and not overlap
#' in time unless allowing ±dt days.
#'
#' @param pathways Output of function read_pathways.
#' @param dt Delta t, ± dt days to determine an indirect contact. If both are set, `indirect_cutoff` wins.
#' @param d0 Duration for a point overlap in time. For example, if two subject were only at the same location
#' for one day or even one hour. Default: 1 (day).
#' @return Edges tibble with from, to, contact_type, direct_days, indirect_days, duration, delta_t, d0. Attribute 'resolution' is set.
#' @export
#
#  Copyright (C) 2025 Yu Wan <yu.wan@liverpool.ac.uk>, Mohammad Saiful Islam Sajib <saiful.sajib@chrfbd.org>
#  Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
#  Creation: 13 November 2025; the latest update: 13 November 2025

find_contacts <- function(pathways, d0 = 1, dt = 0) {

}

.find_contacts <- function(pathways, indirect_cutoff = NULL, dt = NULL, d0 = 1) {
  if (!is.null(dt) && is.null(indirect_cutoff)) indirect_cutoff <- dt
  if (is.null(indirect_cutoff)) indirect_cutoff <- 0
  stopifnot(indirect_cutoff >= 0, d0 >= 0)

  td <- pathways |>
    dplyr::rowwise() |>
    dplyr::mutate(date = list(seq.Date(t_start, t_end, by = 'day'))) |>
    tidyr::unnest(date) |>
    dplyr::ungroup()

  direct <- td |>
    dplyr::inner_join(td, by = c('loc_id','date'), relationship = 'many-to-many') |>
    dplyr::filter(patient_id.x < patient_id.y) |>
    dplyr::count(patient_id.x, patient_id.y, name = 'direct_days')

  p2 <- pathways |>
    dplyr::rename(start2 = t_start, end2 = t_end, pid2 = patient_id)
  p1 <- pathways |>
    dplyr::rename(start1 = t_start, end1 = t_end, pid1 = patient_id)

  ind <- p1 |>
    dplyr::inner_join(p2, by = 'loc_id', relationship = 'many-to-many') |>
    dplyr::filter(pid1 < pid2) |>
    dplyr::filter(
      (end1 < start2 & as.integer(start2 - end1) <= indirect_cutoff) |
      (end2 < start1 & as.integer(start1 - end2) <= indirect_cutoff)
    ) |>
    dplyr::mutate(gap = dplyr::if_else(end1 < start2, as.integer(start2 - end1), as.integer(start1 - end2))) |>
    dplyr::mutate(indirect_days = pmax(indirect_cutoff - gap, 0L)) |>
    dplyr::filter(indirect_days > 0L) |>
    dplyr::group_by(pid1, pid2) |>
    dplyr::summarise(indirect_days = sum(indirect_days), .groups = 'drop')

  edges <- direct |>
    dplyr::full_join(ind, by = c('patient_id.x' = 'pid1', 'patient_id.y' = 'pid2')) |>
    dplyr::mutate(
      direct_days   = dplyr::coalesce(direct_days, 0L),
      indirect_days = dplyr::coalesce(indirect_days, 0L),
      contact_type  = dplyr::case_when(
        direct_days > 0 ~ 'direct',
        indirect_days > 0 ~ 'indirect',
        TRUE ~ NA_character_
      ),
      duration = dplyr::case_when(
        direct_days > 0 ~ direct_days + d0,
        indirect_days > 0 ~ indirect_days + d0,
        TRUE ~ 0
      ),
      delta_t = indirect_cutoff,
      d0 = d0
    ) |>
    dplyr::filter(duration > 0) |>
    dplyr::rename(from = patient_id.x, to = patient_id.y)

  attr(edges, "resolution") <- attr(pathways, "resolution", exact = TRUE)
  attr(edges, "delta_t") <- indirect_cutoff
  attr(edges, "d0") <- d0
  edges
}
