library("GEOquery")
library("plyr")
library("stringr")

geo = getGEO("GSE242252", destdir = "./GSE242252")    # download the main GEO data
geo <- geo[[1]]

# extract the relevant pData
coldata <- pData(geo)[, c("title", "geo_accession", "description")]

# extract groups
coldata$group <- as.factor(coldata$description)
coldata$group <- revalue(coldata$group, c("Healthy control" = "HC", "Herpes Zoster" = "HZ", "Herpes Zoster resolved (1Y)" = "HZres"))

# extract patient identifier of HZ samples
coldata$id <- str_extract(coldata$title, "HZ_\\d+")
coldata[1:22, ]$id <- str_extract(coldata[1:22, ]$title, "CO_\\d+")
coldata$id <- as.factor(coldata$id)

getGEOSuppFiles("GSE242252", baseDir = ".", makeDirectory = TRUE)    # download supplementary files
untar("GSE242252/GSE242252_RAW.tar", exdir = "GSE242252")
countsFiles = list.files("GSE242252", pattern = "txt.gz$", full = TRUE)
library("R.utils")
sapply(countsFiles, gunzip, overwrite = TRUE)    # unzip counts files

countsFiles <- list.files("GSE242252", pattern = "txt$", full = TRUE)
accessionNos <- str_extract(countsFiles, regex("GSM\\d{7}"))

# read in counts
countdata <- read.delim(countsFiles[1], header = FALSE)
colnames(countdata)[2] <- accessionNos[1]

for (i in 2:length(countsFiles)) {
  temp <- read.delim(countsFiles[i], header = FALSE)
  if (sum(countdata[, 1] != temp[, 1]) > 0) {
    print("Gene IDs are not identical!")
    break
  }
  countdata <- cbind(countdata, temp[, 2])
  colnames(countdata)[i+1] <- accessionNos[i]
}

# Remove version from ENSEMBL ID
countdata[, 1] <- str_replace(countdata[, 1], pattern = "\\.\\d+$", replacement = "")

rownames(countdata) <- countdata[, 1]
countdata <- countdata[, -1]

# Remove data from samples not investigated here:
# Remove 9 HC samples used as internal controls in the original publication
countdata <- countdata[, -which(str_detect(coldata$title, regex("^ERC_")))]
coldata <- coldata[-which(str_detect(coldata$title, regex("^ERC_"))), ]

countdata <- countdata[, -which(str_detect(coldata$title, regex("^CO_")) & !(str_detect(coldata$title, regex("-1Y_"))))]
coldata <- coldata[-which(str_detect(coldata$title, regex("^CO_")) & !(str_detect(coldata$title, regex("-1Y_")))), ]

# Remove HC samples and unpaired HZ samples
countdata <- countdata[, -which(coldata$group == "HC",)]
coldata <- coldata[-which(coldata$group == "HC"), ]
# Count the occurrences of each unique id
value_counts <- table(coldata$id)
# Find columns with values that occur only once
unique_values <- names(value_counts[value_counts == 1])
countdata <- countdata[, -which(coldata$id %in% unique_values,)]
coldata <- coldata[-which(coldata$id %in% unique_values), ]
coldata$group <- droplevels(coldata$group)
coldata$id <- droplevels(coldata$id)

table(coldata$group)
# -> Paired data of 23 patients with an acute HZ and one year later following resolution of HZ remaining.

# Create DESeqDataSet
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~ id + group)

# Prefiltering the dataset: At least 2/3 of the samples need to have a read count >5
smallestGroupSize <- 2/3*nrow(coldata)
keep <- rowSums(counts(dds) >= 5) >= smallestGroupSize
dds <- dds[keep,]
nrow(dds)
# -> 16998 genes to analyze

dds <- DESeq(dds)

# Annotate genes
library("org.Hs.eg.db")
library("biomaRt")
library("tidyr")
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

annotate <- function(DESeqResult, returnLists = FALSE) {
  # Remove genes with no log2FC
  DESeqResult <- DESeqResult[!(is.na(DESeqResult$log2FoldChange)), ]

  DESeqResult$symbol <- unname(mapIds(org.Hs.eg.db,
                                      keys=rownames(DESeqResult),
                                      column="SYMBOL",
                                      keytype="ENSEMBL",
                                      multiVals="first"))

  # Try to replace missing symbols using Biomart annotation
  missing <- is.na(DESeqResult$symbol)
  missing <- rownames(DESeqResult[missing, ])
  
  info <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                filters = "ensembl_gene_id",
                values = missing,
                mart = ensembl)
  info <- info[info$hgnc_symbol != "", ] # -> 3810 additional genes
  
  ensembls <- as.data.frame(cbind(rownames(DESeqResult), DESeqResult$symbol))
  colnames(ensembls) <- c("ensembl_gene_id", "hgnc_symbol.1")
  
  ensembls <- join(ensembls, info, by = "ensembl_gene_id", type = "left", match = "all")
  ensembls <- unite(ensembls, "hgnc_symbol", 2:3, sep="", na.rm=TRUE)

  DESeqResult$symbol <- ensembls[, "hgnc_symbol"]

  DESeqResult$genetype <- unname(mapIds(org.Hs.eg.db,
                                 keys=rownames(DESeqResult),
                                 column="GENETYPE",
                                 keytype="ENSEMBL",
                                 multiVals="first"))
  DESeqResult$genename <- unname(mapIds(org.Hs.eg.db,
                                 keys=rownames(DESeqResult),
                                 column="GENENAME",
                                 keytype="ENSEMBL",
                                 multiVals="first"))
  DESeqResult$gene_id <- unname(mapIds(org.Hs.eg.db,
                                keys=rownames(DESeqResult),
                                column="ENTREZID",
                                keytype="ENSEMBL",
                                multiVals="first"))
  DESeqResult
}


# Results HZ vs. HZres
resHZ <- results(dds, contrast = c("group", "HZ", "HZres"), independentFiltering = FALSE, alpha = 0.05)
resHZ <- annotate(resHZ)
summary(resHZ, alpha = 0.05)
# 484 genes up, 680 genes down
write.csv(as.data.frame(resHZ[order(resHZ$padj, decreasing = FALSE),]), file = "HZ_vs_HZres_paired.csv")

# Genes significantly altered in HZ (The numbers between this analysis and the original publication differ
# as we performed a paired analysis of the 23 patients that reappeared one year after the HZ, whereas in the
# original publication an unpaired analysis of 26 initial HZ patients and the 23 patients that reappeared
# was performed. The top DEGs are the same, however.

# Include 'plotVolcano' function from 'functions.R'
# Include function 'removeNoise' from file 'functions.R'
plotVolcano(removeNoise(resHZ), title = "HZ active vs HZ after 1Y", pCutoff = 0.05, FCcutoff = 0.5)


# Compare differential expression patterns between VZV stimulated rPBMC in vitro and HZ reactivation in vivo
# Include 'geneList' function from 'functions.R'
V_vs_O <- geneList("V_vs_O.csv") # The file 'V_vs_0.csv' has been created in 'Analysis of in vitro DEGs.R'
V_vs_O$symbol <- rownames(V_vs_O) # needed for merging
colnames(V_vs_O) <- c("lg2FC.x", "padj.x", "symbol") # in vitro data for the x-axis of the plot

HZres <- geneList("HZ_vs_HZres_paired.csv")
HZres$symbol <- rownames(HZres)
colnames(HZres) <- c("lg2FC.y", "padj.y", "symbol") # in vivo data for the y-axis of the plot

# Get genes present in both expression data sets
merged <- merge(V_vs_O, HZres, by = "symbol")
# -> 11520 genes overlap between V_vs_O and HZ_vs_HZres_paired

# Only plot genes significantly altered by HZ
plotData <- merged[merged$padj.y < 0.05, ]

# Proportion of genes significantly altered in vitro of the genes significantly altered in HZ
sum((plotData$padj.x < 0.01) & (abs(plotData$lg2FC.x) > 1))/nrow(plotData)*100
# -> 11.6 %

# Color dots according to gene expression in the in vitro setting
# No significant change in vitro -> color 0
# Genes significantly different between V vs. O -> color +1
dot.color <- as.integer((plotData$padj.x < 0.01) & (abs(plotData$lg2FC.x) > 1))
# Genes significantly different in vitro and in vivo but with differing signs in the lg2FC -> color = 2
dot.color2 <- as.integer(dot.color == 1 & sign(plotData$lg2FC.x) != sign(plotData$lg2FC.y))
dot.color <- dot.color + dot.color2
dot.color <- as.factor(dot.color)

print(paste("Genes altered in vivo but not in vitro:", sum(dot.color == 0)))
# -> 896
print(paste("Genes altered in vivo and in vitro in same direction:", sum(dot.color == 1)))
# -> 59
print(paste("Genes altered in vivo and in vitro in opposite direction:", sum(dot.color == 2)))
# -> 59
# Names for legend
levels(dot.color) <- c("No effect", "Change in the same direction", "Change in the opposite direction")

# Make dots with lower V vs O lg2FC more transparent
alpha <- abs(plotData[, "lg2FC.x"])
alpha[alpha > 2] = 2
alpha <- (alpha/max(alpha))

# set symbol of genes not significantly different in both investigations to NA, so only the relevant genes will be labelled
plotData[dot.color == 'No effect', ]$symbol <- NA

library("ggrepel")
paste(min(plotData$lg2FC.x), max(plotData$lg2FC.x), min(plotData$lg2FC.y), max(plotData$lg2FC.y))

ggplot(plotData, 
       aes(x = lg2FC.x, y = lg2FC.y)
  ) +
  scale_x_continuous(expression(paste("Effects of VZV antigen stimulation on PBMC ", italic("in vitro"), " [lg2FC]")), 
                     limits = c(-5, 5),
                     breaks = c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5),
                     expand = c(0, 0.2)
  ) +
  scale_y_continuous(expression(paste("Effects of herpes zoster on PBC ", italic("in vivo"), " [lg2FC]")), 
                     limits = c(-0.8, 1.25), 
                     expand = c(0, 0.1)
  ) +
  theme(text = element_text(size=15)
  ) +
  geom_point(
    aes(color = dot.color),
    size = 1, 
    alpha = alpha
  ) +
  scale_color_manual(
    values = c("orange2", "green3", "red2")
  ) +
  geom_text_repel(
    aes(label = symbol),
    size = 2.5,
    min.segment.length = 0.5, 
    segment.color = "grey",
    max.overlaps = 10, 
    direction = "both"
  ) + 
  labs (
    color = expression(paste(italic("In vitro"), " effect: "))
  ) +
  theme(legend.position = "bottom"
  )

# Save genes significantly (p < 0.05) altered in HZ and differentially expressed in the same direction
# in VZV stimulated rPBMC to file 'HZ & VO gene set.csv' (407 genes).
# Add the lg2FC and p-values of the V vs. O comparison for gene enrichment analysis
geneSet <- merged[merged$padj.y < 0.05 & sign(merged$lg2FC.x) == sign(merged$lg2FC.y), 1:3]
geneSet$genetype <- unname(mapIds(org.Hs.eg.db,
                                  keys=geneSet$symbol,
                                  column="GENETYPE",
                                  keytype="SYMBOL",
                                  multiVals="first"))
geneSet$genename <- unname(mapIds(org.Hs.eg.db,
                                  keys=geneSet$symbol,
                                  column="GENENAME",
                                  keytype="SYMBOL",
                                  multiVals="first"))
geneSet$gene_id <- unname(mapIds(org.Hs.eg.db,
                                  keys=geneSet$symbol,
                                  column="ENTREZID",
                                  keytype="SYMBOL",
                                  multiVals="first"))
geneSet <- geneSet[order(geneSet$lg2FC, decreasing = TRUE), ]
write.csv(geneSet, file = "HZ & VO gene set.csv", row.names = FALSE)


# Are the in vitro DEGs due to VZV antigen representative for HZ in vivo?
# GSEA with genes significantly altered in HZ vs. HZres as the in vivo HZ gene set.

# Are the in vitro DEGs due to VZV antigen representative for VZV infection/reactivation as known from the literature?
# Genes upregulated in chickenpox/HZ or being associated with an increased risk for severe chickenpox/HZ
# as known from the literature as the literature gene set.

# In vivo HZ gene set -> Genes significantly altered in HZ with padj < 0.05 (n = 1164).
HZGeneSet <- read.csv(file = "HZ_vs_HZres_paired.csv", header = TRUE)
HZGeneSet <- HZGeneSet[HZGeneSet$padj < 0.05, ]
# Limit HZ gene set to genes that are represented in the in vitro gene expression data (n = 1014)
# Include 'geneList' function from 'functions.R'
V_vs_O <- geneList("V_vs_O.csv") # The file 'V_vs_0.csv' has been created in 'Analysis of in vitro DEGs.R'
HZGeneSet <- HZGeneSet[HZGeneSet$symbol %in% rownames(V_vs_O), ]

# Define gene sets with genes upregulated or downregulated in HZ
HZGenesUp <- HZGeneSet[HZGeneSet$log2FoldChange > 0, ] # 448 genes upregulated
HZGenesDown <- HZGeneSet[HZGeneSet$log2FoldChange < 0, ] # 566 genes downregulated

# Literature gene set (46 genes, 2 are not expressed in the in vitro setting)
library("mclm")
LitGeneSet <- read_txt("VZV_literature_genes")
# Remove empty lines and comments
LitGeneSet <- LitGeneSet[nzchar(LitGeneSet)]
LitGeneSet <- LitGeneSet[!grepl("#", LitGeneSet)]

# For GSEA we need the ENTREZ ID instead of the symbol
LitGeneSet <- mapIds(org.Hs.eg.db,
                     keys = LitGeneSet,
                     column = "ENTREZID",
                     keytype = "SYMBOL",
                     multiVals = "first")

T2G <- rbind(data.frame(term = rep("HZup", nrow(HZGenesUp)),
                        gene = HZGenesUp$gene_id),
             data.frame(term = rep("HZdown", nrow(HZGenesDown)),
                        gene = HZGenesDown$gene_id),
             data.frame(term = rep("Literature", length(LitGeneSet)),
                        gene = LitGeneSet))

# Include 'geneListGSEA' function from 'functions.R'
gene_list <- geneListGSEA("V_vs_O.csv")
library("clusterProfiler")
gseaVO <- GSEA(gene_list, 
               maxGSSize = 1000,
               pvalueCutoff = 1,
               TERM2GENE = T2G)
gseaVO
# Enrichment of genes in the gene set upregulated in HZ with a p-value of 0.0001.
# No enrichment of genes downregulated in V vs. O in the set of genes downregulated in HZ.
# There just is significant enrichment of genes upregulated in vitro in the set of genes downregulated in vivo.
# Highly significant enrichment of the VZV infection/reactivation associated genes from the literature.

# GSEA plot
library("enrichplot")
gseaplot2(gseaVO, geneSetID = 3, title = "Genes upregulated in HZ in vivo", 
          pvalue_table = TRUE,
          pvalue_table_rownames = NULL,
          pvalue_table_columns = c("NES", "pvalue"))
gseaplot2(gseaVO, geneSetID = 2, title = "Genes downregulated in HZ in vivo", 
          pvalue_table = TRUE,
          pvalue_table_rownames = NULL,
          pvalue_table_columns = c("NES", "pvalue"))
gseaplot2(gseaVO, geneSetID = 1, title = "Genes associated with VZV infection/reactivation
according to literature", 
          pvalue_table = TRUE,
          pvalue_table_rownames = NULL,
          pvalue_table_columns = c("NES", "pvalue"))


# Get the core HZ upregulated gene set
HZup_core_genes <- geneInCategory(gseaVO)["HZup"][[1]]
# -> 68 genes out of the 448 genes in the HZup gene set

idx <- which(HZGenesUp$gene_id %in% HZup_core_genes)
HZup_core_genes <- HZGenesUp[idx, ]

# Format dataframe to save it to file
HZup_core_genes$X <- NULL
HZup_core_genes$baseMean <- NULL
HZup_core_genes$log2FoldChange <- NULL
HZup_core_genes$lfcSE <- NULL
HZup_core_genes$stat <- NULL
HZup_core_genes$pvalue <- NULL
HZup_core_genes$padj <- NULL

# Add effect of filgotinib on these 68 genes
gene_list <- geneList("F_vs_V.csv")
gene_list$symbol <- rownames(gene_list)
HZup_core_genes <- merge(HZup_core_genes, gene_list, by = "symbol")

write.csv(HZup_core_genes, file = "Core_gene_set_HZup.csv")
# -> Supplementary table S1


# Plot ridgeplots of the HZ up core gene set and the literature gene set

# Get p-values for the gsea in the F_vs_V and F_vs_O comparisons
gene_list <- geneListGSEA("F_vs_V.csv")
gseaFV <- GSEA(gene_list, 
               maxGSSize = 1000,
               pvalueCutoff = 1,
               TERM2GENE = T2G)
gene_list <- geneListGSEA("F_vs_O.csv")
gseaFO <- GSEA(gene_list, 
               maxGSSize = 1000,
               pvalueCutoff = 1,
               TERM2GENE = T2G)

# Get all lg2FCs of the in vitro assay
# Include 'multipleComparisons' function from 'functions.R'
allLg2FCs <- multipleComparisons()

# Get lg2FCs of the core gene set for different conditions
idx <- which(rownames(allLg2FCs) %in% HZup_core_genes$gene_id)
lg2FCs <- allLg2FCs[idx, ]

# Uncorrected p-values gsea within the 3 comparisons
pvals <- c(gseaVO["HZup"]$pvalue,
           gseaFV["HZup"]$pvalue,
           gseaFO["HZup"]$pvalue)

# Perform BH-correction.
pvals <- p.adjust(pvals, method = "BH")
# p values: VO: 0.0002383726 FV: 0.0001855521 FO: 0.9631284916

# Median lg2FC of the filgotinib effect
median(lg2FCs[, 2])
# -> -0.884501
# Median decrease in gene expression in %:
(1-2^median(lg2FCs[, 2]))*100
# -> 45.8 %

# Include function 'myridgeplot' in 'functions.R'
myridgeplot(coreLg2FCs = lg2FCs, 
            pvals = pvals, 
            title = expression(paste(italic("In vitro"), " differential expression of genes significantly upregulated in HZ")), 
            subtitle = expression(paste("Core set of genes differentially upregulated upon VZV antigen stimulation ", italic("in vitro"))),
            vline = median(lg2FCs[, 2]))

# Get lg2FCs of the literature gene set for different conditions
idx <- which(rownames(allLg2FCs) %in% LitGeneSet)
lg2FCs <- allLg2FCs[idx, ]

# Uncorrected p-values gsea within the 3 comparisons
pvals <- c(gseaVO["Literature"]$pvalue,
           gseaFV["Literature"]$pvalue,
           gseaFO["Literature"]$pvalue)

# Perform BH-correction.
pvals <- p.adjust(pvals, method = "BH")
# p values: VO: 1.500000e-10 FV: 1.875913e-06 FO: 1.500000e-10

# Median lg2FC of the VZV effect
median(lg2FCs[, 1])
# -> 1.608022
# Median lg2FC of the filgotinib effect
median(lg2FCs[, 2])
# -> -0.7491651
# Median decrease in gene expression in %:
(1-2^median(lg2FCs[, 2]))*100
# -> 40.5 %

myridgeplot(coreLg2FCs = lg2FCs, 
            pvals = pvals, 
            title = expression(paste(italic("In vitro"), " differential expression of genes relevant in VZV infection / HZ")), 
            subtitle = "According to the literature",
            vline = median(lg2FCs[, 2]))
