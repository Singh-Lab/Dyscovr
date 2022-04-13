# Script to process String database interactions for a given gene of interest and extract
# network interactions for gene set enrichment analysis
# Written by: Sara Geraghty, April 2022

# Download String networks from website: https://string-db.org/cgi/network?taskId=bSXrcGoHJTto&sessionId=beGJRBX9wxZ4

# Import the String network file
string_nw <- read.csv(paste0(getwd(), "/Desktop/string_interactions.tsv"),
                      header = TRUE, check.names = FALSE, sep = "\t")
colnames(string_nw)[1] <- "node1"

# Limit to given GOI
goi <- "TP53"
# goi <- "PIK3CA"
string_nw_goi <- string_nw[(string_nw$node1 == goi),]  # will already be subsetted by our max # of interactions

# Visualize the scores for this GOI
string_nw_goi <- string_nw_goi[order(string_nw_goi$combined_score),]
hist(string_nw_goi$combined_score)

# Get the list of genes that made the cutoff for top interactions with our GOI
interactors <- string_nw_goi$node2

# Write these interactors to a file
write.table(interactors, paste0(getwd(), paste0("/Desktop/", paste0(goi, "_interactors_string_top100.txt"))),
            quote = FALSE, row.names = FALSE)

