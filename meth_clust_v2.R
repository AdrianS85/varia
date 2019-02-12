library(data.table)
library(foreach)



#Load .cov files
cov_tibble_list3 <- lapply(list.files(pattern = "^P_(.*)bismark.cov$"), FUN = data.table::fread, sep = "\t")



# Here we copy the list, to have the original one for the reference
c_cov_tibble_list3 <- data.table::copy(cov_tibble_list3)



# Here we prepare ID, cov and meth columns
lapply(X = c_cov_tibble_list3, FUN = function(X) {
  X[, ':=' (org_coord = paste0(V1, "_", V2), next_coord = paste0(V1, "_", V2 + 1), cov = V5 + V6, meth = V4)] #rename and prepare new columns
  X[, (grep(pattern = "V.", x = colnames(X), value = TRUE)) := NULL] # rm old columns
} )



# Here we prepare the "next coordinate" table. In this one dont need "next coord" column
n_cov_tibble_list3 <- data.table::copy(c_cov_tibble_list3)
lapply(X = n_cov_tibble_list3, FUN = function(X) {
  X[, next_coord := NULL]
} )



# Here we merge orginal and +1 coordinate for each .cov. This could use parallel, but it would be easy and its fast either way
merged_cov_tibble_list3  <- mapply(
  function(x,y)
    { merge(x, y, by.x = "next_coord", by.y = "org_coord", all.x = T)}, 
  x = c_cov_tibble_list3, 
  y = n_cov_tibble_list3,
  SIMPLIFY = F)



# Here we establish filenames
cov_names <- stringr::str_remove(list.files(pattern = "^P_(.*)bismark.cov$"), "_S(.*)")



# Here we tidy up and calculate weighted mean of cytosines in each CpG
lapply(X = merged_cov_tibble_list3, FUN = function(x) {
  as.data.table(x)
  setcolorder(x, neworder= c("org_coord", "cov.x", "meth.x", "next_coord", "cov.y", "meth.y"))
  x[!is.na(meth.y), meth_weight := (((cov.x * meth.x) + (cov.y * meth.y)) / (cov.x + cov.y)) ] # Here we calculate the weighted mean methylation level
  x[is.na(meth.y), meth_weight := meth.x ]
} )



# Here we meticulously rename columns with sample IDs
foreach(n = seq_along(merged_cov_tibble_list3)) %do% 
  setnames(
    x = merged_cov_tibble_list3[[n]], 
    old = colnames(merged_cov_tibble_list3[[n]]), 
    new = c("ID", 
            paste0("cov_org", "_", cov_names[n]),
            paste0("meth_org", "_", cov_names[n]),
            "next_coord",
            paste0("cov_next", "_", cov_names[n]),
            paste0("meth_next", "_", cov_names[n]),
            paste0("meth_weight", "_", cov_names[n])
  ) )




# Here we write CpG - combined files
foreach(n = seq_along(merged_cov_tibble_list3)) %do% 
  data.table::fwrite(x = merged_cov_tibble_list3[[n]], file = paste0("all_info_", cov_names[n], ".tsv"), sep = "\t")



# Here we produce only-CpGID/weighted mean of cytosines in each CpG files
lapply(X = merged_cov_tibble_list3, FUN = function(X) {
  X[, (grep(pattern = "next", x = colnames(X), value = TRUE)) := NULL]
  X[, (grep(pattern = "_org", x = colnames(X), value = TRUE)) := NULL]
} )


# Merge and collapse list into single table
comb_cov_tibble_list <- Reduce(function(x, y) merge(x, y, by = "ID", all = T), merged_cov_tibble_list3)


data.table::fwrite(c_cov_tibble_list, "brain_all_cov.tsv") # brain liver placenta



# Prepare list of all significant CpG sites from all comparisons? Perhaps from one comparison at a time?
rnbead_sites <- lapply(list.files(pattern = "^br_*"), FUN = readr::read_csv, col_names = T, col_types = "-cc---n-n") ## Change ^br to "_site_"
names_sites <- str_remove(list.files(pattern = "^br_"), "_diff(.*)")  ## Change ^br to "_site_"



# Filter the list to get only pval < 0.05 and diff higher than 5%. This wierd diffmeth.p.val > 1 is cause RnBeads throws some numbers as "e-9" or smth. And they are not read properly by R.
filtered_rnbead_sites <- rnbead_sites %>% 
  map(filter, (diffmeth.p.val < 0.05 | diffmeth.p.val > 1) & abs(mean.diff) > 0.05) %>%
  map(mutate, ID = paste(Chromosome, Start, sep = "_"))
# Add comparison naming column to site file                          
for (n in seq(from = 1, to = length(names_sites))){
  filtered_rnbead_sites[[n]]$comparison <- names_sites[n]
}



# Here we get singular table with all CpGs along with their diff and pval
merged_filtered_rnbead_sites <- Reduce(function(x, y) merge(x, y, by = "ID", all=TRUE), filtered_rnbead_sites)
write_tsv(merged_filtered_rnbead_sites, "brain_diff_pval_sites.tsv") # brain liver placenta



#Now we merge diffsites and raw methylation values
# Create pvaldiff with only IDs
clus_merged_filtered_rnbead_sites <- data.frame(merged_filtered_rnbead_sites$ID, stringsAsFactors = F)
colnames(clus_merged_filtered_rnbead_sites) <- "ID"

for_clustering <- merge(clus_merged_filtered_rnbead_sites, c_cov_tibble_list, by = "ID", all.x = T)
data.table::fwrite(for_clustering, "brain_for_clustering_cov.tsv") # brain liver placenta                               

#http://bonsai.hgc.jp/~mdehoon/software/cluster/cluster3.pdf                                      
# INTO BASH: ../cluster -f brain_for_clustering_cov.tsv -m a -g 7 

### INTO R

### PREPARE DIFFERENTIALLY METHYLATED CPG FILES
### Move diffmeth cpg files into the working folder. Rename them so that they are recognizable
rename 's/^diffMethTable/b1_diffMethTable/' diffMethTable*
