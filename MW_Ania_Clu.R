https://bioconductor.org/packages/release/bioc/html/multiClust.html   
#clustComp - clustComp is an open source Bioconductor package that implements different techniques for the comparison of two gene expression clustering results. These include flat versus flat and hierarchical versus flat comparisons.



https://ocw.mit.edu/courses/electrical-engineering-and-computer-science/6-0002-introduction-to-computational-thinking-and-data-science-fall-2016/lecture-videos/index.htm
BiocManager::install("lattice")
BiocManager::install("multiClust", version = "3.8")    
BiocManager::install("runibic", version = "3.8")    
library(multiClust)



library(tidyverse)

#setwd("E:/Projekty/W - MWieczorek/Results for MW/For Anna new clustering vagotomy/Nowe od Anny 1118/Clustering")
setwd("/media/adrians/USB DISK1/Projekty/W - MWieczorek/Results for MW/For Anna new clustering vagotomy/Nowe od Anny 1118/Clustering")


COMPARISONS <- readr::read_tsv("OnlyVagxClusteringF.txt", col_types = "cnnnnnnnnnnnnnnnnnnnnn")
## Tutaj nadal s? zduplikowane geny
COMP <- COMPARISONS %>%
  mutate(Gene_Symbol = str_remove(Genes, ".* ")) %>%
  mutate(ID = str_remove(Genes, " .*")) #%>%
  #filter(Gene_Symbol != "NA")


## Tutaj mamy plik z sondami z nowymi anotacjami
NEW_ANNOT_PROBE <- readr::read_tsv("From OnlyVagotomy v ania k.txt", col_types = "cc")


## Niestety trzeba jeszcze zanotowa? COMPARISONS with newer probe names
TEST_CLU <- readr::read_tsv("CIV.txt", col_types = "cc")
# Tutaj dodamy sondy do gen?w
T_CLU <- merge(TEST_CLU, NEW_ANNOT_PROBE, by = "Gene_Symbol", all.x = T)


COMP_check_TEST_CLU <- merge(T_CLU, COMP, by = "ID", all.x = T) %>%
  mutate(ID_Symbol = paste(Gene_Symbol.x, ID, sep = "_"))

FINAL_COMP_check_TEST_CLU <- COMP_check_TEST_CLU %>%
  mutate(
    ID = NULL,
    Gene_Symbol.x = NULL,
    Gene_Symbol.y = NULL,
    Cluster = NULL,
    Genes = NULL
  ) %>%
  select(ID_Symbol, everything())

ROWNAMES_FINAL_COMP_check_TEST_CLU <- as.data.frame(FINAL_COMP_check_TEST_CLU)
rownames(ROWNAMES_FINAL_COMP_check_TEST_CLU) <- ROWNAMES_FINAL_COMP_check_TEST_CLU$ID_Symbol
ROWNAMES_FINAL_COMP_check_TEST_CLU$ID_Symbol <- NULL

matrixROWNAMES_FINAL_COMP_check_TEST_CLU <- as.matrix(ROWNAMES_FINAL_COMP_check_TEST_CLU)


gplots::heatmap.2(matrixROWNAMES_FINAL_COMP_check_TEST_CLU, Colv = F, margins = c(8, 8), trace = "none")


#





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
BiocManager::install(c("runibic", "biclust"))
library(runibic)



###### Clustering 
set_runibic_params(t = 0.85, q = 0, f = 1, nbic = 100L, div = 0L, useLegacy = FALSE)

runibicRES2 <- biclust::biclust(matrixROWNAMES_FINAL_COMP_check_TEST_CLU, method = BCUnibic())
###### Clustering 



###### Visualizations  
# Simultaneous representation of biclusters. multidimensional scaling based on the gene and condition profiles
biclust::bubbleplot(matrixROWNAMES_FINAL_COMP_check_TEST_CLU, runibicRES2, showLabels = T)

# Returns actual heatmap
biclust::drawHeatmap(matrixROWNAMES_FINAL_COMP_check_TEST_CLU, bicResult = runibicRES2, number = 1, beamercolor = F, paleta = col_for_heatmaps(1000))

# Które próbki poszły do którego z klastrów. 
# Basically if in a column of of these stacked rectangles, a rectangle is coloured, it means that this condition is part of that particular bicluster. Now, if the Mid box is not checked, a coloured rectangle consists out of two parts, left and right.  The left colour represents the mean of this condition for all the genes within the biclusters.  However, the right colour contains the global mean value for this condition.  If the Mid option is checked though, the rectangle exist out of three colours with the global mean in the middle and the bicluster mean on the left and right.
biclustmember(runibicRES2, matrixROWNAMES_FINAL_COMP_check_TEST_CLU)


#plotclust(runibicRES2, matrixROWNAMES_FINAL_COMP_check_TEST_CLU, noC = 6)
#heatmapBC(matrixROWNAMES_FINAL_COMP_check_TEST_CLU, runibicRES2)
###### Visualizations 



###### Quality measures  

# Differential co-expression framework to quantify goodness of biclusters and compare biclustering algorithms
BC_QM_CK <- biclust::ChiaKaruturi(matrixROWNAMES_FINAL_COMP_check_TEST_CLU, runibicRES2, number = 3)

BC_QM_CV <- biclust::constantVariance(matrixROWNAMES_FINAL_COMP_check_TEST_CLU, runibicRES2, number = 3)
BC_QM_AV <- biclust::additiveVariance(matrixROWNAMES_FINAL_COMP_check_TEST_CLU, runibicRES2, number = 3)
BC_QM_MV <- biclust::multiplicativeVariance(matrixROWNAMES_FINAL_COMP_check_TEST_CLU, runibicRES2, number = 3)
BC_QM_SV <- biclust::signVariance(matrixROWNAMES_FINAL_COMP_check_TEST_CLU, runibicRES2, number = 3)

BC_QM_COF <- biclust::computeObservedFstat(matrixROWNAMES_FINAL_COMP_check_TEST_CLU, runibicRES2, number = 1)

BC_QM_DCR <- biclust::diagnoseColRow(matrixROWNAMES_FINAL_COMP_check_TEST_CLU, runibicRES2, number = 1, nResamplings = 200)
biclust::ddiagnosticPlot(BC_QM_DCR)


BC_QM_DT <- biclust::diagnosticTest(runibicRES2, matrixROWNAMES_FINAL_COMP_check_TEST_CLU, number = 1:runibicRES2@Number, statistics = "Tukey", save_F = T)
biclust::diagnosticPlot2(BC_QM_DT, number = 1, StatVal = T, binwidth = NULL)

###### Quality measures 



###### Get results 

# Here we get actual result tables
summary(runibicRES2)

writeBiclusterResults(fileName = "ass", bicResult = runibicRES2, 
                      bicName = "asshole",
                      geneNames = rownames(matrixROWNAMES_FINAL_COMP_check_TEST_CLU), 
                      arrayNames = colnames(matrixROWNAMES_FINAL_COMP_check_TEST_CLU),
                      append =FALSE, 
                      delimiter="\t")

writeclust(runibicRES2, row = TRUE, noC = 3)

test1 <- bicluster(matrixROWNAMES_FINAL_COMP_check_TEST_CLU, runibicRES2, number= 1:runibicRES2@Number)

###### Get results 
############ runibic/biclust ############ 



############ MCbiclust ############ 

# MCbiclust Nucleic Acids Res. 2017 "genes from large datasets with maximal correlated gene expression, allowing regulation of complex networks to be examined" for large scale gene expression sets. easy input
BiocManager::install("MCbiclust")
library(MCbiclust)
data(CCLE_small)
xxx <- head(CCLE_small)
#For reproducibility set.seed has been used to set R’s pseudo-random number generator. It should also be noted that the for gem the data matrix can not contain all the genes, since FindSeed involves the calculation of correlation matrices which are not computationally efficient to compute if they involve greater than ~1000 genes.

# I guess this does the clustering...?
set.seed(1)
seed <- FindSeed(matrixROWNAMES_FINAL_COMP_check_TEST_CLU, seed.size = 10, iterations = 5000, initial.seed = 100, messages = 1000)

# I guess this calculates correlation matrix for samples for heatmap.2 function...? Dunno why?
CCLE.random.cor <- cor(t(matrixROWNAMES_FINAL_COMP_check_TEST_CLU[,seed])) 
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



############ cogena ############ 
#cogena - bardzo dobre do klasycznego hierarchicznego klastrowania, daje zrozumiały, ładny output
# http://homes.di.unimi.it/valentini/SlideCorsi/MB0910/HierarchicalClustering.pdf
# https://www.biostars.org/p/74223/

BiocManager::install("cogena")
library(cogena)
data(Psoriasis, package = "cogena")

enecl_result <- cogena::coExp(obj = matrixROWNAMES_FINAL_COMP_check_TEST_CLU, 
                      nClust = (2:6), 
                      clMethods = "hierarchical", 
                      metric = "abscorrelation", 
                      method = "average", 
                      ncore = 2)


############ cogena ############ 
