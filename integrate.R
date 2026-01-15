# "Normal" means the pediatric kidney adjacent to the tumor and "Human" means
# the developing kidney

require(Seurat)
rna <- readRDS('data/pat4580d5/pat4580d5_snRNAseq_normxpatientxhfk_integrated_V1.rds')

require(SingleCellExperiment)
require(scuttle)
rna <- as.SingleCellExperiment(rna)
rownames(rna) <- make.names(rownames(rna)) # replace "-" with "."
assay(rna,'logcounts') <- log1p(t(t(counts(rna)) / colSums(counts(rna)))*1000)
require(scry)
rna <- devianceFeatureSelection(rna, batch = factor(rna$orig.ident))


# just patient 4580d5 for now
s1 <- readRDS("data/pat4580d5/pat4580d5_seqFISH_ROI1.rds")
s1 <- as.SingleCellExperiment(s1)
s1 <- assay(s1,'logcounts')
colnames(s1) <- paste0('s1_', colnames(s1))
s2 <- readRDS("data/pat4580d5/pat4580d5_seqFISH_ROI2.rds")
s2 <- as.SingleCellExperiment(s2)
s2 <- assay(s2,'logcounts')
colnames(s2) <- paste0('s2_', colnames(s2))
s3 <- readRDS("data/pat4580d5/pat4580d5_seqFISH_ROI3.rds")
s3 <- as.SingleCellExperiment(s3)
s3 <- assay(s3,'logcounts')
colnames(s3) <- paste0('s3_', colnames(s3))
s4 <- readRDS("data/pat4580d5/pat4580d5_seqFISH_slice2.rds")
s4 <- as.SingleCellExperiment(s4)
s4 <- assay(s4,'logcounts')
colnames(s4) <- paste0('s4_', colnames(s4))
# dimension reduction "SG" is the original coordinates

# checks
all(rownames(s1) == rownames(s2))
all(rownames(s1) == rownames(s3))
all(rownames(s1) == rownames(s4))

rownames(rna)[rownames(rna)=='MARCH11'] <- 'MARCHF11'

all(rownames(s1) %in% rownames(rna))

# Genes to keep:
# - 2000 most variable, 
# - ~250 from spatial, 
# - TFs from Ng et al.
# - L-R pair genes from CellChat
TFs <- readRDS('data/TFs.rds')
LRs <- readRDS('data/LRs.rds')

keep <- which(rowData(rna)$binomial_deviance >= sort(rowData(rna)$binomial_deviance, decreasing = TRUE)[2000] |
                  rownames(rna) %in% rownames(s1) |
                  rownames(rna) %in% TFs |
                  rownames(rna) %in% LRs)

assay_list <- list(rna = assay(rna,'logcounts')[keep,],
                   s1 = s1,
                   s2 = s2,
                   s3 = s3,
                   s4 = s4)

require(StabMap)

mosaicDataUpSet(assay_list, plot = FALSE)

# cleanup
rm(rna,s1,s2,s3,s4,keep,LRs,TFs)

# get joint embedding
stab <- stabMap(assay_list, reference_list = c("rna"),
                suppressMessages = FALSE, maxFeatures = nrow(assay_list$rna),
                plot = TRUE)




