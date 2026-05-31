#' @keywords internal
#' @noRd
#' @author Yu Wan, \email{yu.wan@liverpool.ac.uk}
#
#  Create a separate classes file to avoid ROxygen2's scoping problem, which occurs when setClass commands are placed before
#  the function definition, creating a "[class name].rd" file instead of the function documentation.
#
#  Copyright (C) 2025-2026 Yu Wan <yu.wan@liverpool.ac.uk>, Mohammad Saiful Islam Sajib <saiful.sajib@chrfbd.org>
#  Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
#  Creation: 10 May 2026; the latest update: 31 May 2026

# For pathopath.R ###############
setClass(
    Class = "Network",
    slots = list(
        V = "data.frame",
        E = "data.frame"
        )
    )

setClass(
    # The output class of function pathpath
    Class = "Pathopath",
    slots = list(
        pathways = "list",
        migrations = "data.frame",
        contacts = "data.frame",
        summary = "data.frame",
        network = "Network",
        parameters = "list"
        )
    )

# For read_pathways.R ###############
setClass(
    Class = "Pathways",
    slots = list(
        pathways = "list",
        migrations = "data.frame"
        )
    )

# For compute_contacts.R ###############
setClass(
    Class = "Periods",
    slots = list(
        start = "Date",
        end = "Date",
        length = "integer"
        )
    )

# For hd_clustering.R ###############
setClass(
    Class = "Clusters",
    slots = list(
        clusters = "data.frame",
        cluster_counts = "data.frame",
        tree = "list",
        distances = "matrix"
        )
    )

# For find_components.R ###############
setClass(
    Class = "Components",
    slots = list(
        V = "data.frame",
        E = "data.frame",
        membership = "data.frame",
        component_size = "data.frame"
    )
)
