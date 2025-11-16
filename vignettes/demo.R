# Demonstration of pathopath's utility
#  Copyright (C) 2025 Yu Wan <yu.wan@liverpool.ac.uk>, Mohammad Saiful Islam Sajib <saiful.sajib@chrfbd.org>
#  Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
#  Creation: 16 November 2025; the latest update: 16 November 2025

library(dplyr)
library(readr)
library(tidyr)
library(pathopath)

pp <- pathopath(pathway_data = "vignettes/movements.tsv", dt = 3)

# Explore individual slots
print(slotNames(pp))  # Five slot names

pathways <- pp@pathways
migrations <- pp@migrations
contacts <- pp@contacts |> filter(Contact != "None")
contact_summary <- pp@summary
V <- pp@network@V
E <- pp@network@E

# Export results
write_tsv(contacts, file = "vignettes/contacts.tsv")
write_tsv(contact_summary, file = "vignettes/contact_summary.tsv")
write_tsv(V, file = "vignettes/network_nodes.tsv")
write_tsv(E, file = "vignettes/network_edges.tsv")

saveRDS(pp, file = "vignettes/pp.rds")
