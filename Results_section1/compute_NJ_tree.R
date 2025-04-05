"""
==============================================================================
Dependencies:
    - ape           : For working with phylogenetic trees.
    - spaa          : For similarity table to similarity matrix conversion.

Version Information:
    - R:            4.2.1
    - ape:          5.8
    - spaa:         0.2.2
        
Notes:
    - The script assumes ADMIXTOOLS output format with the value of f3 in the 4th column, 
    the outgroup in the 3rd column, and population pairs in the 1st and 2nd columns.

==============================================================================
"""

library(ape)
library(spaa)

df = read.table("outgroupf3.results", header=FALSE)

# Keep only pairs and f3 value
df <- df[, -c(3, 5:7)]

# Fix symmetries 
for (i in 1:length(df[, 1])) {
  df <- rbind(df, data.frame("V1"=df[i, 2], "V2"=df[i, 1], "V4"=df[i, 3]))
}

# Transform from similarity table to similarity matrix
M_dist <- as.matrix(list2dist(df))

# Invert values to obtain a distance matrix 
M_dist <- 1/M_dist

# Fix diagonal 
diag(M_dist) <- 0

# Exclusions (comment out if not needed)
exclusions <- grepl("Sumidouro|Anzick|12000BP", row.names(M_dist)) 
# Filter M_dist
M_dist <- M_dist[!exclusions, !exclusions]

# Convert to distance matrix
M_dist_r <- as.dist(M_dist)

# Compute NJ tree
tree <- nj(M_dist_r)

# Root tree
rooted_tree <- root(tree, outgroup="USA_Ancient_Beringian.SG", resolve.root=TRUE)

# Write tree
write.tree(rooted_tree, file = "outf3_rooted_noSumidouro_noAnzick_noLosRieles12000BP.tree", append = FALSE,
           digits = 10, tree.names = FALSE)
