# https://www.r-bloggers.com/working-with-venn-diagrams/

library(gtools)
install.packages("VennDiagram")


# Here we upload files with gene names (column 2) and comparison name (column 1). One file for each comparison
genes_list <- lapply(list.files(pattern = "*placenta_.txt"), FUN = read.table, header = T) # liver placenta brain


# Here we produce list of character vectors to venn (from column 2 of each file) and names (from column 1 of each file). We need named list for intersection comparison
vector_list <- list()
name_vector <- NA
for (n in seq(from = 1, to = length(genes_list))){
  vector_list[[n]] <- as.character(genes_list[[n]][[2]])
  name_vector[n] <- as.character(genes_list[[n]][[1,1]])
}
names(vector_list) <- name_vector

# Here is the crucial part: We produce list of genes belonging to each intersection of comparisons!
inters <- attr(gplots::venn(vector_list, show.plot=FALSE),"intersections")
rlist::list.save(x = inters, file = "placenta_intersections_names_list.json")
# liver_ placenta_ brain_
# convert json to excel here https://json-csv.com/
