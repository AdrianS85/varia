library(data.table)
library(foreach)

#### Prepare list of all .cov data ####

#Load .cov files
cov_tibble_list3 <- lapply(list.files(pattern = "^P_(.*)bismark.cov$"), FUN = data.table::fread, sep = "\t")


# Here we copy the list, to have the original one for the reference
c_cov_tibble_list3 <- data.table::copy(cov_tibble_list3)


# Here we prepare ID, cov and meth columns
lapply(X = c_cov_tibble_list3, 
       FUN = function(X) { 
         X[, ':=' (org_coord = paste0(V1, "_", V2), next_coord = paste0(V1, "_", V2 + 1), cov = V5 + V6, meth = V4)] #rename and prepare new columns
         X[, (grep(pattern = "V.", x = colnames(X), value = TRUE)) := NULL] # rm old columns
}  )


# Here we prepare the "next coordinate" table.
n_cov_tibble_list3 <- data.table::copy(c_cov_tibble_list3)
lapply(X = n_cov_tibble_list3, FUN = function(X) {
  X[, next_coord := NULL]
} )


# Here we merge orginal and +1 coordinate for each .cov. This could use parallel, but its fast as it is, fortunatly
merged_cov_tibble_list3  <- mapply(
  function(x,y)
    { merge(x, y, by.x = "next_coord", by.y = "org_coord", all.x = T)}, 
  x = c_cov_tibble_list3, 
  y = n_cov_tibble_list3,
  SIMPLIFY = F)
rm(c_cov_tibble_list3, n_cov_tibble_list3)
gc()


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


# Here we write CpG - combined files (could try %dopar%, but something didnt work I think)
foreach(n = seq_along(merged_cov_tibble_list3)) %do% 
  data.table::fwrite(x = merged_cov_tibble_list3[[n]], file = paste0("all_info_", cov_names[n], ".tsv"), sep = "\t")


# Here we produce only-CpGID/weighted mean of cytosines in each CpG files
lapply(X = merged_cov_tibble_list3, FUN = function(X) {
  X[, (grep(pattern = "next", x = colnames(X), value = TRUE)) := NULL]
  X[, (grep(pattern = "_org", x = colnames(X), value = TRUE)) := NULL]
} )


# Merge and collapse list into single table
comb_cov_tibble_list <- Reduce(function(x, y) merge(x, y, by = "ID", all = T), merged_cov_tibble_list3)


data.table::fwrite(comb_cov_tibble_list, "placenta_all_cov.tsv", sep = "\t") # brain liver placenta
           
#### Prepare list of all .cov data ####

                               
                               
                               
                               
                               
#### Prepare list of all .cov data PLUS ONE, because RnBeads reports only first cytosine from pair and in .cov they are mixed ####
                               
#Here we need to add the same meth levels with site +1, because sometimes in .cov we have only data on second cytosine, and RnBeads still annotates such data as the first cytosine                                       
additional_sites <- as.data.table(comb_cov_tibble_list)

                               
# Here we add a column with location minus (plus?) one, so that we can capture sites which had only second cytosine sequenced
additional_sites[, plus_one := (paste0(stringr::str_split_fixed(string = additional_sites$ID, pattern = "_", n = 2)[,1], "_", as.numeric(stringr::str_split_fixed(string = additional_sites$ID, pattern = "_", n = 2)[,2]) - 1) )  ]
additional_sites$ID <- NULL
additional_sites[, ":=" (ID = plus_one, plus_one = NULL)]
ordered_additional_sites <- additional_sites %>% dplyr::select(ID, dplyr::everything())
rm(additional_sites)
     
                               
# Add additional logic to add to ***comb_cov_tibble_list*** only those rows from ***placenta_more_sites_comb_cov_tibble_list.tsv*** which are absent in ***comb_cov_tibble_list***
ordered_additional_sites[ID %in% comb_cov_tibble_list$ID, ID := "NA"]
                               
# Merge orginal and PLUS ONE sites and write them                             
more_sites_comb_cov_tibble_list <- rbind(ordered_additional_sites, comb_cov_tibble_list)                                       
data.table::fwrite(more_sites_comb_cov_tibble_list, "placenta_more_sites_comb_cov_tibble_list.tsv", sep = "\t")

#### Prepare list of all .cov data PLUS ONE, because RnBeads reports only first cytosine from pair and in .cov they are mixed ####
                               
                               

                               
                               
                               
#### Prepare list of all significant CpG sites from all comparisons? Perhaps from one comparison at a time? ####
                               
rnbead_sites <- lapply(list.files(pattern = "^pl_*"), FUN = readr::read_csv, col_names = T, col_types = "-cc---n-n--n--n---nn") ## Change ^br to "_site_"

                               
## Here we change tables do data.tables for quicker manipulations later on
rnbead_sites <- parallel::mclapply(X = rnbead_sites, mc.cores = 12, FUN = function(X) { data.table::as.data.table(X) } )
names_sites <- stringr::str_remove(list.files(pattern = "^pl_"), "_diff(.*)")  ## Change ^br to "_site_"


# Filter the list to get only pval < 0.05 and diff higher than 5% and more. This wierd diffmeth.p.val > 1 is cause RnBeads throws some numbers as "e-9" or smth. And they are not read properly by R.
filtered_rnbead_sites <- parallel::mclapply(X = rnbead_sites, mc.cores = 12,
FUN = function(X) {
subset(X, (diffmeth.p.val < 0.001 | diffmeth.p.val > 1) & abs(mean.diff) > 0.1 & X[[7]] < 3 & X[[8]]  < 3 ) 
})

                               
#Add ID column 
filtered_rnbead_sites <- parallel::mclapply(X = filtered_rnbead_sites , mc.cores = 12,
FUN = function(X) {
X[, ID := paste(Chromosome, Start, sep = "_")]
X[, ':=' (Chromosome = NULL, Start = NULL)]
} )
   
                               
# Add comparison naming column to site file and save the object                         
for (n in seq(from = 1, to = length(names_sites))){
  filtered_rnbead_sites[[n]]$comparison <- names_sites[n]
}
rlist::list.save(x = filtered_rnbead_sites, file = "placenta_filtered_rnbead_sites_list.json") # brain liver placenta


# Extract ID columns into separate tables                                
ID_filtered_rnbead_sites <- list()
for (n in seq(from = 1, to = length(names_sites))){
  ID_filtered_rnbead_sites[[n]] <- as.data.table(as.character(filtered_rnbead_sites[[n]]$ID) ) }

                               
# Here we get singular table with all CpGs
merged_ID_filtered_rnbead_sites <- data.frame(unique(rlist::list.rbind(ID_filtered_rnbead_sites)), stringsAsFactors = F)
colnames(merged_ID_filtered_rnbead_sites) <- "ID"

#### Prepare list of all significant CpG sites from all comparisons? Perhaps from one comparison at a time? ####                             

                               
                               

                               
                               
# Merge .cov-derived data with RnBeads-derived IDs and write them
for_clustering <- merge(merged_ID_filtered_rnbead_sites, more_sites_comb_cov_tibble_list, by = "ID", all.x = T)
data.table::fwrite(for_clustering, "placenta_for_clustering_cov.tsv", sep = "\t") # brain liver placenta

# Try clustering                               
h_test <- hclust(d = dist(for_clustering), method = "average") ### Cant do it cause missing data
                               
#https://www.r-bloggers.com/drawing-heatmaps-in-r/
#http://genomicsclass.github.io/book/pages/clustering_and_heatmaps.html
#https://www.rdocumentation.org/packages/RFLPtools/versions/1.6/topics/write.hclust
#http://bonsai.hgc.jp/~mdehoon/software/cluster/cluster3.pdf                                      
# INTO BASH: ../cluster -f brain_for_clustering_cov.tsv -m a -g 7 

### INTO R

### PREPARE DIFFERENTIALLY METHYLATED CPG FILES
### Move diffmeth cpg files into the working folder. Rename them so that they are recognizable
rename 's/^diffMethTable/b1_diffMethTable/' diffMethTable*
