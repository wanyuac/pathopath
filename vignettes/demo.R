# Demonstration of pathopath's utility
#  Copyright (C) 2025-2026 Yu Wan <yu.wan@liverpool.ac.uk>, Mohammad Saiful Islam Sajib <saiful.sajib@chrfbd.org>
#  Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
#  Creation: 16 November 2025; the latest update: 31 May 2026

library(dplyr)
library(readr)
library(tidyr)
library(pathopath)

# Section 1: build an initial contact network ###############
pp <- pathopath(movements = "vignettes/input_movements.tsv",
                samples = "vignettes/input_samples.tsv",
                dt = 3)

# Explore individual slots ===============
print(slotNames(pp))  # Six slot names

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
pp_updated <- add_movements(movements = "vignettes/input_movements_additional.tsv",
                            samples = "vignettes/input_samples_updated.tsv",
                            previous_results = pp)
print(slotNames(pp_updated))  # Six slot names
saveRDS(pp_updated, file = "vignettes/output_pp_updated.rds")
View(pp_updated@network@V)

# Section 3: demonstrate error messages from the quality assessment of input movement data ###############
incorrect_pathways <- read_pathways(movement_table = "vignettes/input_movements_with_mistakes.tsv")

# Section 4: clustering of nodes based on contact lengths ###############
# Scenario 1: no edge pruning
contact_clusters <- contact_clustering(E = pp@network@E, V = pp@network@V)
g <- igraph::graph_from_data_frame(d = contact_clusters@E, directed = FALSE, vertices = contact_clusters@V)
View(contact_clusters@membership)
plot(g)

# Scenario 2: pruning edges for a maximum of contact length of 10 days
contact_clusters_10 <- contact_clustering(E = pp@network@E, V = pp@network@V, l_max = 10)
g_10 <- igraph::graph_from_data_frame(d = contact_clusters_10@E, directed = FALSE, vertices = contact_clusters_10@V)
View(contact_clusters_10@membership)
plot(g_10)

# Section 5: clustering of nodes based on Hamming distances ###############

# Method 1: complete-linkage hierarchical clustering ===============
library(ape)

hc <- h_clustering(mat = "vignettes/input_distance_matrix.tsv", thresholds = 2)
View(hc@clusters)  # Show membership of subjects under distance thresholds 0 and 2
View(hc@cluster_counts)  # Show the number of clusters under distance thresholds 0 and 2
plot(hc@tree[["tree"]])  # Draw the dendrogram
write_tsv(hc@clusters, file = "vignettes/output_h_clustering.tsv")
write.tree(hc@tree[["tree"]], file = "vignettes/output_h_clustering_dend.newick")

# Method 2: component discovery in a distance network ===============
comp <- dn_clustering(mat = "vignettes/input_distance_matrix.tsv", thresholds = 2)
View(comp@clusters)  # Membership of subjects
View(comp@cluster_counts)  # Number of components under distance thresholds 0 and 2
write_tsv(comp@clusters, file = "vignettes/output_dn_clustering.tsv")
