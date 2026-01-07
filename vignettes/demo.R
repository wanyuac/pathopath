# Demonstration of pathopath's utility
#  Copyright (C) 2025-2026 Yu Wan <yu.wan@liverpool.ac.uk>, Mohammad Saiful Islam Sajib <saiful.sajib@chrfbd.org>
#  Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
#  Creation: 16 November 2025; the latest update: 7 January 2026

library(dplyr)
library(readr)
library(tidyr)
library(pathopath)

# Section 1: build an initial contact network ###############
pp <- pathopath(movements = "vignettes/input_movements.tsv", dt = 3)

# Explore individual slots ===============
print(slotNames(pp))  # Five slot names

pathways <- pp@pathways
migrations <- pp@migrations
contacts <- pp@contacts |> filter(Contact != "None")
contact_summary <- pp@summary
V <- pp@network@V
E <- pp@network@E

# Export results as individual data files ===============
write_tsv(migrations, file = "vignettes/output_migrations.tsv")
write_tsv(contacts, file = "vignettes/output_contacts.tsv")
write_tsv(contact_summary, file = "vignettes/output_contact_summary.tsv")
write_tsv(V, file = "vignettes/output_network_nodes.tsv")
write_tsv(E, file = "vignettes/output_network_edges.tsv")
saveRDS(pp, file = "vignettes/output_pp.rds")
saveRDS(pathways, file = "vignettes/output_pathways.rds")

# Section 2: add movements to the network ###############
updated_results <- add_movements(movements = "vignettes/input_movements_additional.tsv", previous_results = pp)
print(slotNames(updated_results))
saveRDS(updated_results, file = "vignettes/updated_results.rds")
