## BacterialCore.py outputs a PCG table already. However, the initial communities we
## are using here to start our simulations don't have the same PCGs as the endpoint
## communities (already adapted to glucose) and that we are interested in.
##
## Therefore, this script takes those OTUs in our original samples and classifies them
## as Node35562, Node27828 or "others" using ape to check if they are within those
## nodes or not

# tr0
table_file <- "simulation_parameters/original_table.from_biom_0.99.txt"
outfile    <- "simulation_parameters/pcg_leaves.txt"
# tr1
table_file <- "simulation_parameters/tr1_X2_X6_table.from_biom_0.99.txt"
outfile    <- "simulation_parameters/pcg_leaves_tr1.txt"

library(ape)

# Load the tree
tree_file <- "simulation_parameters/99_otus_nodes.tree"
phy_tree <- read.tree(tree_file)


leaves_node27828 <- extract.clade(phy_tree, "Node27828")$tip.label
leaves_node35562 <- extract.clade(phy_tree, "Node35562")$tip.label

# Load the abund table
table_data <- read.table(table_file, skip = 1, sep = "\t", row.names = 1)

# Save intersections
leaves_string_node27828 <- paste(intersect(leaves_node27828, rownames(table_data)), collapse = ";")
leaves_string_node35562 <- paste(intersect(leaves_node35562, rownames(table_data)), collapse = ";")

# Save the rest as "others"
leaves_string_others <- paste(rownames(table_data)[!(rownames(table_data) %in% c(leaves_string_node27828, leaves_string_node35562))], collapse = ";")

# Create the pcgtable
out_data <- data.frame(Core = c("Node27828", "Node35562", "others"),
                       Leaves = c(leaves_string_node27828, leaves_string_node35562, leaves_string_others),
                       Percentage = c(0.714995632150111, 0.210402382979819, 0.07460198487007))

write.table(out_data, file = outfile, sep = "\t", quote = FALSE, row.names = F)
