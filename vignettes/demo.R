# Demonstration of pathopath's utility
#  Copyright (C) 2025 Yu Wan <yu.wan@liverpool.ac.uk>, Mohammad Saiful Islam Sajib <saiful.sajib@chrfbd.org>
#  Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
#  Creation: 16 November 2025; the latest update: 16 November 2025

library(pathopath)
library(dplyr)

pp <- pathopath(pathway_data = "vignettes/movements.tsv", dt = 3)

# Explore individual slots
pathways <- pp@pathways
migrations <- pp@migrations
Contacts <- pp@contacts |> filter(Contact != "None")
