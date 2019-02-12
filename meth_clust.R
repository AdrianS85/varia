############ CLUSTERING ###########
### Get all bedgraph files into single matrix
#In folder with all .bedgraph.gz files
### Get bedtools
wget https://github.com/arq5x/bedtools2/releases/download/v2.27.1/bedtools-2.27.1.tar.gz
tar jxvf bedtools-2.27.1.tar.gz
cp bedtools2/bin/* /usr/local/bin

### Prepare sorted bedgraph files
gzip -dk *gz
ls *bedGraph > AllLongNames; sed 's/_R1.*//' AllLongNames > AllShortNames; paste AllShortNames AllLongNames > AllCompareNames
parallel --colsep '\t' "sort -k 1,1 -k2,2n {2} > {1}.clust.sort.bedGraph" :::: AllCompareNames ##sort a BED file by chromosome then by start position https://bedtools.readthedocs.io/en/latest/content/tools/sort.html


### PREPARE list of all captured CpG sites in all bedgraph files in given tissue
ls B*.clust.sort.bedGraph | sudo tee BrainMultiinterOrder
sudo bedtools multiinter -header -i B*.clust.sort.bedGraph | sudo tee BrainAllCpgs ##Requires that each interval file is sorted by chrom/start
ls L*.clust.sort.bedGraph | sudo tee LiverMultiinterOrder
bedtools multiinter -header -i L*.clust.sort.bedGraph > LiverAllCpgs
ls L*.clust.sort.bedGraph | sudo tee PlacentaMultiinterOrder
bedtools multiinter -header -i P*.clust.sort.bedGraph > PlacentaAllCpgs

### INTO R:
x <- readr::read_tsv(file = "BrainAllCpgs", col_types = "cnnnc") ### Read only relevant columns
x <- x %>% dplyr::()
x2 <- x[, 4:6]

### PREPARE .bedgraph files
## https://stackoverflow.com/questions/14096814/merging-a-lot-of-data-frames




### NEW ###

#Load .cov files
cov_tibble_list3 <- lapply(list.files(pattern = "^P_(.*)bismark.cov$"), FUN = data.table::fread, sep = "\t")



# Here we copy the list, to have the original one for the reference
library(data.table)
c_cov_tibble_list3 <- data.table::copy(cov_tibble_list3)




cl = parallel::makeCluster(2, type = "PSOCK")
parallel::parLapply(cl, X = c_cov_tibble_list3, FUN = function(X) {
  X[, ':=' (org_coord = paste0(V1, "_", V2), next_coord = paste0(V1, "_", V2 + 1), cov = V5 + V6, meth = V4)] #rename and prepare new columns
  X[, (grep(pattern = "V.", x = colnames(X), value = TRUE)) := NULL] # rm old columns
} )
parallel::stopCluster(cl)




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

bc_merged_cov_tibble_list3 <- copy(merged_cov_tibble_list3)
merged_cov_tibble_list3 <- copy(bc_merged_cov_tibble_list3)

#Here we tidy up tables a bit and add weighted mean methylation!
lapply(X = merged_cov_tibble_list3, FUN = function(x) {
  as.data.table(x)
  setcolorder(x, neworder= c("org_coord", "cov.x", "meth.x", "next_coord", "cov.y", "meth.y"))
  setnames(x = x, old = c("org_coord", "cov.x", "meth.x", "next_coord", "cov.y", "meth.y"), new = c("org_coord", "cov_org", "meth_org", "next_coord", "cov_next", "meth_next"))
  x[!is.na(meth_next), meth_weight := (((cov_org * meth_org) + (cov_next * meth_next)) / (cov_org + cov_next)) ] # Here we calculate the weighted mean methylation level
  x[is.na(meth_next), meth_weight := meth_org ]
} )

### NEW ###









cov_tibble_list <- lapply(list.files(pattern = "^B_(.*)bismark.cov$"), FUN = readr::read_tsv, col_names = F, col_types = "cnnn") ## B_ L_ P_
a_cov_tibble_list <- cov_tibble_list %>% map(mutate, ID = paste(X1, X2, sep = "_"))
b_cov_tibble_list <- a_cov_tibble_list %>% map(select, ID, X4)
rm(a_cov_tibble_list)
## Name bedgraph columns with filenames
cov_names <- str_remove(list.files(pattern = "^B_(.*)bismark.cov$"), "_S(.*)") ## B_ L_ P_
for (n in seq(from = 1, to = length(cov_names))) { colnames(b_cov_tibble_list[[n]]) <- c("ID", cov_names[n]) }

c_cov_tibble_list <- Reduce(function(x, y) merge(x, y, by = "ID", all = T), b_cov_tibble_list)
write_tsv(c_cov_tibble_list, "brain_all_cov.tsv") # brain liver placenta

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
write_tsv(for_clustering, "brain_for_clustering_cov.tsv") # brain liver placenta                               

 #http://bonsai.hgc.jp/~mdehoon/software/cluster/cluster3.pdf                                      
 # INTO BASH: ../cluster -f brain_for_clustering_cov.tsv -m a -g 7 
                                       
### INTO R

### PREPARE DIFFERENTIALLY METHYLATED CPG FILES
### Move diffmeth cpg files into the working folder. Rename them so that they are recognizable
rename 's/^diffMethTable/b1_diffMethTable/' diffMethTable*


############ CLUSTERING ###########
