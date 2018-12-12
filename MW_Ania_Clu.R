#https://bioconductor.org/packages/release/bioc/html/multiClust.html   
#clustComp - clustComp is an open source Bioconductor package that implements different techniques for the comparison of two gene expression clustering results. These include flat versus flat and hierarchical versus flat comparisons.


# ROWNAMES_FINAL_COMP_check_TEST_CLU <- as.data.frame(FINAL_COMP_check_TEST_CLU)
# rownames(ROWNAMES_FINAL_COMP_check_TEST_CLU) <- ROWNAMES_FINAL_COMP_check_TEST_CLU$ID_Symbol
# ROWNAMES_FINAL_COMP_check_TEST_CLU$ID_Symbol <- NULL
# 
# matrixed_annot_the_list$ID_Symbol <- NULL
#          
# test <- annot_the_list[[4]]
# class(annot_the_list[[2]])
# annot_the_list <- merge(T_CLU, COMP, by = "ID", all.x = T) %>%
#   mutate()
# 
# ?select_helpers
# 
# ## Tutaj mamy plik z sondami z nowymi anotacjami
# NEW_ANNOT_PROBE <- readr::read_tsv("From OnlyVagotomy v ania k.txt", col_types = "cc")
# 
# 
# ## Niestety trzeba jeszcze zanotowa? COMPARISONS with newer probe names
# TEST_CLU <- readr::read_tsv("CIV.txt", col_types = "cc")
# # Tutaj dodamy sondy do gen?w
# T_CLU <- merge(TEST_CLU, NEW_ANNOT_PROBE, by = "Gene_Symbol", all.x = T)
# 
# 
# COMP_check_TEST_CLU <- merge(T_CLU, COMP, by = "ID", all.x = T) %>%
#   mutate(ID_Symbol = paste(Gene_Symbol.x, ID, sep = "_"))
# 
# FINAL_COMP_check_TEST_CLU <- COMP_check_TEST_CLU %>%
#   mutate(
#     ID = NULL,
#     Gene_Symbol.x = NULL,
#     Gene_Symbol.y = NULL,
#     Cluster = NULL,
#     Genes = NULL
#   ) %>%
#   select(ID_Symbol, everything())
# 
# ROWNAMES_FINAL_COMP_check_TEST_CLU <- as.data.frame(FINAL_COMP_check_TEST_CLU)
# rownames(ROWNAMES_FINAL_COMP_check_TEST_CLU) <- ROWNAMES_FINAL_COMP_check_TEST_CLU$ID_Symbol
# ROWNAMES_FINAL_COMP_check_TEST_CLU$ID_Symbol <- NULL
# 
# matrixROWNAMES_FINAL_COMP_check_TEST_CLU <- as.matrix(ROWNAMES_FINAL_COMP_check_TEST_CLU)



#https://ocw.mit.edu/courses/electrical-engineering-and-computer-science/6-0002-introduction-to-computational-thinking-and-data-science-fall-2016/lecture-videos/index.htm





library(tidyverse)

#setwd("E:/Projekty/W - MWieczorek/Results for MW/For Anna new clustering vagotomy/Nowe od Anny 1118/Clustering")
setwd("/media/adrians/USB DISK1/Projekty/W - MWieczorek/Results for MW/For Anna new clustering vagotomy/Nowe od Anny 1118/Clustering")













#### Here we prepare file with gene/probenames and expression values ####

COMP <- readr::read_tsv("OnlyVagxClusteringF.txt", col_types = "cnnnnnnnnnnnnnnnnnnnnn") %>%
  mutate(Gene_Symbol = str_remove(Genes, ".* ")) %>%
  mutate(Transcript_ID = str_remove(Genes, " .*")) #%>%
  #filter(Gene_Symbol != "NA")

#### Here we prepare file with gene/probenames and expression values ####



#### Here we prepare input files annotated with expression values ####

# Here we have dataframe list with probe-gene names as first column

## This is what we want to cluster. First column with Cluster number, second column with Gene_Symbol, third with Transcript_ID
the_list <- lapply(list.files(pattern = "^C.*txt"), FUN = read_tsv, col_names = T, col_types = "ccc")
the_names <- str_remove(list.files(pattern = "^C.*txt"), ".txt")

## Here we add expression values to input files
annot_the_list <- the_list %>%
  map(.f = merge, y = COMP, by = "Transcript_ID", all.x = T) %>%
  map(mutate, ID_Symbol = paste(Gene_Symbol.x, Transcript_ID, sep = "_")) %>%
  map(select, ID_Symbol, matches(match = "CNT"), matches(match = "SHAM"), matches(match = "VAGX"))

## Here I make matrix-like dataframe list with probe-gene names in rownames
matrixed_annot_the_list <- annot_the_list 
for (n in seq(1,9)){
  rownames(matrixed_annot_the_list[[n]]) <- matrixed_annot_the_list[[n]]$ID_Symbol
  matrixed_annot_the_list[[n]]$ID_Symbol <- NULL
  matrixed_annot_the_list[[n]] <- as.matrix(matrixed_annot_the_list[[n]])
}

# Here lets output the tables with raw values
dir.create("ClusterExpressionValues")
for (n in seq(1,9)) {
  write.table(x = annot_the_list[[n]], file = paste0("ClusterExpressionValues/", the_names[n], "_Expression_Values.txt"), sep = "\t", dec = ",")
}

#### Here we prepare input files annotated with expression values ####



###### FUNCTION FOR GETTING ENRICHMENTS ######
### tibble_list - tibble list with tibbles for all 5 comparisons, 
### Currently I work with tibbles with 3 columns, with the first one being gene names and the third one being comparison_id. We should be able to use different ones
function_enrich_list <- function(tibble_list = enrich_overview_list, gene_names_column = 1){
  ### Create list of interesting databases
  dbs_Onto_Path <- c("GO_Molecular_Function_2018", "GO_Cellular_Component_2018", "GO_Biological_Process_2018", "MGI_Mammalian_Phenotype_2017", "Human_Phenotype_Ontology", "KEGG_2016", "WikiPathways_2016", "Panther_2016", "Reactome_2016", "BioCarta_2016", "NCI-Nature_2016", "ARCHS4_Kinases_Coexp", "HumanCyc_2016", "BioPlex_2017", "SILAC_Phosphoproteomics")
  dbs_Regul <- c("Genome_Browser_PWMs", "TRANSFAC_and_JASPAR_PWMs", "Transcription_Factor_PPIs", "ChEA_2016",  "TF-LOF_Expression_from_GEO", "PPI_Hub_Proteins", "ENCODE_TF_ChIP-seq_2015", "ENCODE_Histone_Modifications_2015", "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X", "CORUM", "Pfam_InterPro_Domains", "Phosphatase_Substrates_from_DEPOD", "TF_Perturbations_Followed_by_Expression", "ARCHS4_TFs_Coexp", "miRTarBase_2017", "TargetScan_microRNA_2017", "Enrichr_Submissions_TF-Gene_Coocurrence", "Epigenomics_Roadmap_HM_ChIP-seq")
  dbs_Drug_Tissue_Other <- c("ARCHS4_IDG_Coexp", "DrugMatrix", "RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO", "OMIM_Disease", "Jensen_DISEASES", "DSigDB",  "Jensen_COMPARTMENTS", "ARCHS4_Tissues", "Tissue_Protein_Expression_from_Human_Proteome_Map", "Tissue_Protein_Expression_from_ProteomicsDB", "Mouse_Gene_Atlas", "ESCAPE", "Chromosome_Location", "MSigDB_Computational")
  dbs <- c(dbs_Onto_Path, dbs_Regul, dbs_Drug_Tissue_Other)
  
  ##  dbs_Removed <- c("Jensen_TISSUES", "dbGaP", "Genes_Associated_with_NIH_Grants", "GeneSigDB")
  ##  Error in scan(file = file, what = what, sep = sep, quote = quote, dec = dec,  : 
  ##  line 1 did not have 2 elements OR line 109 did not have 9 elements
  ##  https://www.biostars.org/p/349240/
  
  ### Initialize list, so I can work on it
  enriched <- list()
  
  for(n in seq(from = 1, to = length(comparisons_ids))){
    enriched[[n]] <- enrichR::enrichr(tibble_list[[n]][[gene_names_column]], dbs) %>%
      map(filter, P.value <= 0.05) %>%
      map(as.tibble)
  }
  
  ### Initialize list, so I can work on it
  enriched_names <- list()
  ### Initialize list counter for the loop
  merged_number_2 <- 1
  
  ### Add information on databaseses and comparisons used
  for(n in seq(from = 1, to = length(comparisons_ids))){
    for(d in seq(from = 1, to = length(dbs))){
      enriched_names[[merged_number_2]] <- enriched[[n]][[d]]%>%
        mutate(Database = dbs[d], Comparison = comparisons_ids[n])
      merged_number_2 <- merged_number_2 + 1  
    }
  }
  
  merged_enriched_names <- rlist::list.rbind(enriched_names)
  
  ### Create list of tibbles based on comparison names
  final_list_2 <- list()
  for (c in seq(from = 1, to = length(comparisons_ids))){
    final_list_2[[c]] <- merged_enriched_names %>%
      filter(Comparison == comparisons_ids[c])
  }
  
  return(final_list_2)
}
###### FUNCTION FOR GETTING ENRICHMENTS ######



###### FUNCTION FOR WRITING TABLES FROM LISTS OF REGIONS ######
### tibble_list = list of tibbles from RnBeads differential module results, analyzed as neccecary
### analysis_name = name of analysis, tissue, region, string
function_write_list <- function(annotated_tibble_list, analysis_name){
  
  for(n in seq(from = 1, to = length(comparisons_ids))){
    write.table(annotated_tibble_list[n], paste0(annotated_tibble_list[[n]]$Comparison[1], "_", analysis_name, ".txt"), sep="\t")
  } 
  
  return()  
}
###### FUNCTION FOR WRITING TABLES FROM LISTS OF REGIONS ######



#### Here we will enrich original gene name lists ####

library(enrichR)

dbs_Onto_Path <- c("GO_Molecular_Function_2018", "GO_Cellular_Component_2018", "GO_Biological_Process_2018", "MGI_Mammalian_Phenotype_2017", "Human_Phenotype_Ontology", "KEGG_2016", "WikiPathways_2016", "Panther_2016", "Reactome_2016", "BioCarta_2016", "NCI-Nature_2016", "ARCHS4_Kinases_Coexp", "HumanCyc_2016", "BioPlex_2017", "SILAC_Phosphoproteomics")
dbs_Regul <- c("Genome_Browser_PWMs", "TRANSFAC_and_JASPAR_PWMs", "Transcription_Factor_PPIs", "ChEA_2016",  "TF-LOF_Expression_from_GEO", "PPI_Hub_Proteins", "ENCODE_TF_ChIP-seq_2015", "ENCODE_Histone_Modifications_2015", "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X", "CORUM", "Pfam_InterPro_Domains", "Phosphatase_Substrates_from_DEPOD", "TF_Perturbations_Followed_by_Expression", "ARCHS4_TFs_Coexp", "miRTarBase_2017", "TargetScan_microRNA_2017", "Enrichr_Submissions_TF-Gene_Coocurrence", "Epigenomics_Roadmap_HM_ChIP-seq")
dbs_Drug_Tissue_Other <- c("Jensen_TISSUES", "ARCHS4_IDG_Coexp", "DrugMatrix", "RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO", "OMIM_Disease", "Jensen_DISEASES", "DSigDB",  "Jensen_COMPARTMENTS", "ARCHS4_Tissues", "Tissue_Protein_Expression_from_Human_Proteome_Map", "Tissue_Protein_Expression_from_ProteomicsDB", "Mouse_Gene_Atlas", "ESCAPE", "Chromosome_Location", "MSigDB_Computational", "dbGaP", "Genes_Associated_with_NIH_Grants", "GeneSigDB")
dbs <- c(dbs_Onto_Path, dbs_Regul, dbs_Drug_Tissue_Other)

comparisons_ids = seq(1,9)

enrichr_list1 <- function_enrich_list(tibble_list = the_list, gene_names_column = 2)

dir.create("SimpleEnrichment")
for(n in seq(1,9)){
  write.table(enrichr_list1[[n]], paste0("SimpleEnrichment/", the_names[n], "_enrichment.txt"), sep="\t")
} 

#### Here we will enrich original gene name lists ####



#### Basic hierachical clustering just to visualize expression patterns #### 

dir.create("BasicClustering")
for (n in seq(1,9)) {
  png(filename = paste0("BasicClustering/", the_names[n], "_Cluster_Summary.png"), width = 1080, height = 5000, units = "px")
  par(cex.axis = 1, oma=c(4,0,0,0))
  gplots::heatmap.2(matrixed_annot_the_list[[n]], Colv = F, trace = "none")
  dev.off()
}

#### Basic hierachical clustering just to visualize expression patterns #### 



############ 
install.packages("viridis")
library(viridis)
palette("viridis")
library(RColorBrewer)
col_for_heatmaps <- colorRampPalette(brewer.pal(3, "RdYlBu"))

par(cex.axis = 0.5, mar = c(5, 5, 5, 5) + 5)  ## Set graphical options for the plot

############ runibic/biclust ############ 
#runibic - UniBic biclustering algorithm for gene expression data scientif reports 2016 uses biclust package. looks good.
#Biclustering is considered NP-hard as it investigates relations between multiple rows that occur in different subsets of columns. The running time of the algorithms is usually highly dependent on the size of the input data.
#One of the recent breakthroughs in gene expression analysis was development of UniBic (Wang et al., 2016). The algorithm originally implemented in C managed to capture biologically meaningful trend-preserving patterns and proved to outperform multiple other methods. 
#The algorithm, originally released as sequential, has proven to outperform multiple popular biclustering state-of-the-art biclustering methods (Wang et al., 2016).
#Results returned from runibic are wrapped into a Biclust object, which can be used for further examination, including visualization and analysis provided by biclust package.
#Similarly, the biclust method could be applied to any matrix extracted from ExpressionSet using exprs() function.
#Apart from allowing analysis of genomic data from historical ExpressionSet, runibic package is compatible with SummarizedExperiment class (Morgan et al., 2017).
# https://uhdspace.uhasselt.be/dspace/bitstream/1942/22018/1/GuideBiclustGUI.pdf
#BiocManager::install(c("runibic", "biclust", "gplots"))
library(runibic)



###### Clustering 
#set_runibic_params(t = 0.85, q = 0, f = 1, nbic = 100L, div = 0L, useLegacy = FALSE)

runibicRES2 <- matrixed_annot_the_list %>%
  map(biclust::biclust, method = BCUnibic())

###### Clustering 



###### Visualizations 

dir.create("ClusterSummaryRunibic")
for (n in seq(1,9)) {
  png(filename = paste0("ClusterSummaryRunibic/", the_names[n], "_Cluster_Summary.png"), width = 1920, height = 1080, units = "px")
  par(cex.axis = 1, oma=c(0,4,0,8))
  biclustmember(runibicRES2[[n]], matrixed_annot_the_list[[n]], main = the_names[n])
  dev.off()
}

#plotclust(runibicRES2, matrixROWNAMES_FINAL_COMP_check_TEST_CLU, noC = 6)
#heatmapBC(matrixROWNAMES_FINAL_COMP_check_TEST_CLU, runibicRES2)

# Simultaneous representation of biclusters. multidimensional scaling based on the gene and condition profiles
#biclust::bubbleplot(matrixed_annot_the_list[[1]], runibicRES2[[1]], showLabels = T)

# Returns actual heatmap
#biclust::drawHeatmap(matrixed_annot_the_list[[1]], runibicRES2[[1]], number = 1, beamercolor = F, paleta = col_for_heatmaps(1000))

# Które próbki poszły do którego z klastrów. 
# Basically if in a column of of these stacked rectangles, a rectangle is coloured, it means that this condition is part of that particular bicluster. Now, if the Mid box is not checked, a coloured rectangle consists out of two parts, left and right.  The left colour represents the mean of this condition for all the genes within the biclusters.  However, the right colour contains the global mean value for this condition.  If the Mid option is checked though, the rectangle exist out of three colours with the global mean in the middle and the bicluster mean on the left and right.

###### Visualizations 



###### Quality measures  

# Differential co-expression framework to quantify goodness of biclusters and compare biclustering algorithms
BC_QM_CK <- biclust::ChiaKaruturi(matrixed_annot_the_list[[1]], runibicRES2[[1]], number = 3)

BC_QM_CV <- biclust::constantVariance(matrixed_annot_the_list, runibicRES2, number = 3)
BC_QM_AV <- biclust::additiveVariance(matrixed_annot_the_list, runibicRES2, number = 3)
BC_QM_MV <- biclust::multiplicativeVariance(matrixed_annot_the_list, runibicRES2, number = 3)
BC_QM_SV <- biclust::signVariance(matrixed_annot_the_list, runibicRES2, number = 3)

BC_QM_COF <- biclust::computeObservedFstat(matrixed_annot_the_list[[3]], runibicRES2[[1]], number = 1)

BC_QM_DCR <- biclust::diagnoseColRow(matrixed_annot_the_list[[3]], runibicRES2[[1]], number = 1, nResamplings = 200)
biclust::ddiagnosticPlot(BC_QM_DCR)


BC_QM_DT <- biclust::diagnosticTest(runibicRES2, matrixed_annot_the_list, number = 1:runibicRES2@Number, statistics = "Tukey", save_F = T)
biclust::diagnosticPlot2(BC_QM_DT, number = 1, StatVal = T, binwidth = NULL)

###### Quality measures 



###### Get results 

## Here we write result tables with all clusters including genes and comparisons inside them
dir.create("ClustersRunibic")
for (n in seq(1,9)) {
  writeBiclusterResults(fileName = paste0("ClustersRunibic/", the_names[n], "_ClustersRunibic.txt"), 
                        bicResult = runibicRES2[[n]], 
                        bicName = paste0(the_names[n]),
                        geneNames = rownames(matrixed_annot_the_list[[n]]), 
                        arrayNames = colnames(matrixed_annot_the_list[[n]]),
                        append =FALSE, 
                        delimiter="\t")}

## Here we can get given cluster with gene names and expression values
test1 <- bicluster(matrixed_annot_the_list[[1]], runibicRES2[[1]], number= 1)

#writeclust(runibicRES2[[1]], row = TRUE, noC = 3)

###### Get results 

############ runibic/biclust ############ 







############ cogena ############ 
#cogena - bardzo dobre do klasycznego hierarchicznego klastrowania, daje zrozumiały, ładny output
#An example dataset of Psoriasis. This dataset is used for illustration of the usage of cogena package. It has been normalised the expression profling using rma method, filtered some non-informative genes using MetaDE package and analysed the differentially expressed genes using limma package with the cut-off adjuested p-value 0.05 and abs(logFC) >=1
# http://homes.di.unimi.it/valentini/SlideCorsi/MB0910/HierarchicalClustering.pdf
# https://www.biostars.org/p/74223/

BiocManager::install("cogena")
library(cogena)

# Getting (obsolete, 3 year-old) annotation files
annoGMT <- "c5.bp.v5.0.symbols.gmt.xz"
annofile <- system.file("extdata", annoGMT, package="cogena")

# Creating our own sampleLabel. It is a named factor vector
ourSampleLabel <- factor(str_remove(colnames(matrixROWNAMES_FINAL_COMP_check_TEST_CLU), "[_].*"))
names(ourSampleLabel) <- colnames(matrixROWNAMES_FINAL_COMP_check_TEST_CLU)

# Creating data matrix with unique, avaraged gene names
FINAL_COMP_check_TEST_CLU$ID_Symbol <- str_remove(FINAL_COMP_check_TEST_CLU$ID_Symbol, "[_].*") # Get gene names without Probe IDs

# Get medians of all non-uniwue gene names
AVG.NAMES_FINAL_COMP_check_TEST_CLU <- FINAL_COMP_check_TEST_CLU %>%
  group_by(ID_Symbol) %>%
  summarise_all(median) 

# Put gene names into rownames and make it a matrix
AVG.NAMES_FINAL_COMP_check_TEST_CLU <- as.data.frame(AVG.NAMES_FINAL_COMP_check_TEST_CLU)
rownames(AVG.NAMES_FINAL_COMP_check_TEST_CLU) <- AVG.NAMES_FINAL_COMP_check_TEST_CLU$ID_Symbol
AVG.NAMES_FINAL_COMP_check_TEST_CLU$ID_Symbol <- NULL
AVG.NAMES_FINAL_COMP_check_TEST_CLU <- as.matrix(AVG.NAMES_FINAL_COMP_check_TEST_CLU)

# Set gene names to all upper, because cogena databases use this type of symbols (human?)
rownames(AVG.NAMES_FINAL_COMP_check_TEST_CLU) <- toupper(rownames(AVG.NAMES_FINAL_COMP_check_TEST_CLU))

### Clustering and enrichment. nClust oznacza chyba ile razy mamy sklastrować dane, czy coś takiego. Koniecznie przeczytać o co chodzi.
genecl_result <- cogena::coExp(obj = AVG.NAMES_FINAL_COMP_check_TEST_CLU, 
                      nClust = (2:6), 
                      clMethods = "hierarchical", 
                      metric = "abscorrelation", 
                      method = "average", 
                      ncore = 2)


clen_res <- clEnrich(genecl_result, annofile=annofile, sampleLabel=ourSampleLabel)

#test1 <- sota(AVG.NAMES_FINAL_COMP_check_TEST_CLU, maxCycles = 4, distance = "correlation")

### Clustering and enrichment.



### Reports

# summary(clen_res)
# enrichment.table <- enrichment(clen_res, "hierarchical", "5")
# geneInCluster(clen_res, "hierarchical", "5", "1")
# gene2set(annofile, rownames(AVG.NAMES_FINAL_COMP_check_TEST_CLU))
# geneclusters(genecl_result, "hierarchical", 4)

geneExp <- geneExpInCluster(clen_res, "hierarchical", "5") # This is what we want i think. Gene names + clusters + expression
geneExpDT <- geneExp[[1]]

### Reports



### Figures

dev.off() # naprawia rysowanie obrazków
heatmapCluster(clen_res, "hierarchical", "6", maintitle="2") 

heatmapCmap(clen_res, "hierarchical", "3", orderMethod="All")

heatmapPEI(clen_res, "hierarchical", "2") # 

corInCluster(clen_res, "hierarchical", "5", "1") # Figures with correlation between gene expression patterns in cluster

### Figures

############ cogena ############ 







############ MCbiclust ############ 

# MCbiclust Nucleic Acids Res. 2017 "genes from large datasets with maximal correlated gene expression, allowing regulation of complex networks to be examined" for large scale gene expression sets. easy input
BiocManager::install("rpart")
library(MCbiclust)
#For reproducibility set.seed has been used to set R’s pseudo-random number generator. It should also be noted that the for gem the data matrix can not contain all the genes, since FindSeed involves the calculation of correlation matrices which are not computationally efficient to compute if they involve greater than ~1000 genes.

# I guess this does the clustering...?
set.seed(1)
seed <- FindSeed(matrixed_annot_the_list[[1]], seed.size = 10, iterations = 5000, initial.seed = 100, messages = 1000)

# I guess this calculates correlation matrix for samples for heatmap.2 function...? Dunno why?
CCLE.random.cor <- cor(t(matrixed_annot_the_list[,seed])) 
gplots::heatmap.2(CCLE.random.cor, trace = "none")

# automatically selects those genes that are most strongly associated with the pattern. It works by using hierarchical clustering to select the genes into n different groups, and then discarding any of these groups that fail to have a correlation score greater than the correlation score from all the genes together.
CCLE.hicor.genes <- as.numeric(HclustGenesHiCor(matrixROWNAMES_FINAL_COMP_check_TEST_CLU, seed, cuts = 8))
CCLE.mito.cor2 <- cor(t(matrixROWNAMES_FINAL_COMP_check_TEST_CLU[CCLE.hicor.genes, seed]))
CCLE.heat <- gplots::heatmap.2(CCLE.mito.cor2,trace = "none")

#Getting gene names from heatmap...
CCLE.groups <- list(labels(CCLE.heat$rowDendrogram[[1]]),
                    labels(CCLE.heat$rowDendrogram[[2]]))

#In this example a distinct correlation pattern was found. However this was only examined for genes involved in mitochondrial function. Non-mitochondrial genes are likely also involved in this pattern and it is important to identify them.
CCLE.cor.vec <- CVEval(gem.part = matrixROWNAMES_FINAL_COMP_check_TEST_CLU,
                       gem.all = matrixROWNAMES_FINAL_COMP_check_TEST_CLU,
                       seed = seed, splits = 10)

GOEnrichmentAnalysis

#Already all the genes in the data set have had the correlation calculated to the pattern found. One more task that can be readily done is to order the samples according to the strength of correlation. 
CCLE.samp.sort <- SampleSort(matrixROWNAMES_FINAL_COMP_check_TEST_CLU[as.numeric(CCLE.hicor.genes),], seed = seed)

#Once the samples have been sorted it is possible to summarise the correlation pattern found using principal component analysis (PCA).PCA is a method of dimensional reduction, and converts a data set to a new set of variables known as the principal components. These are designed to be completely uncorrelated or orthogonal to each other. In this way the principal components are new variables that capture the correlations between the old variables, and are in fact a linear combination of the old variables. PC1 captures the highest variance within the data, so if PCA is run on the found bicluster with very strong correlations between the genes, PC1 will be a variable that summarises this correlation.
pc1.vec <- PC1VecFun(top.gem = matrixROWNAMES_FINAL_COMP_check_TEST_CLU[as.numeric(CCLE.hicor.genes),], seed.sort = CCLE.samp.sort, n = 10)

#So far MCbiclust outputs a ranked list of genes and samples. In many cases it is however necessary to determine which genes and samples are within the bicluster and which are not. 
CCLE.bic <- ThresholdBic(cor.vec = CCLE.cor.vec,
                         sort.order = CCLE.samp.sort,
                         pc1 = pc1.vec, samp.sig = 0.05)

#Once this thresholded bicluster has been found it is important to properly align the PC1 vector and the correlation vector such that samples with a high PC1 values are those samples with up-regulated genes that have positive CV values. This is not strictly necessary to do, but makes the interpretation of MCbiclust simpler.
pc1.vec <- PC1Align(gem = matrixROWNAMES_FINAL_COMP_check_TEST_CLU, pc1 = pc1.vec,
                    sort.order = CCLE.samp.sort,
                    cor.vec = CCLE.cor.vec, bic = CCLE.bic)


############ MCbiclust ############ 







############ pam ############ 
cluster::pam()


############ pam ############ 

