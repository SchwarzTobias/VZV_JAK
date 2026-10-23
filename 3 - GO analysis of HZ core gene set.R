# Which GO terms do the genes significantly upregulated in HZ and in VZV stimulated rPBMC represent?

# Perform GO overrepresentation analysis with HZ core up gene set
topGenes <- read.csv("Core_gene_set_HZup.csv", header=TRUE)
# Include 'geneListGSEA' function from file 'Functions.R'.
# We need the list of gene IDs from the in vitro assays which we load from file using the geneListGSEA function.
geneNames <- names(geneListGSEA("F_vs_V.csv"))
geneList <- factor(as.integer(geneNames %in% topGenes$gene_id))
names(geneList) <- geneNames

library("topGO")

GOdata <- new("topGOdata", 
              description = "HZ core up geneset", 
              ontology = "BP",
              allGenes = geneList,
              nodeSize = 10,
              annot = annFUN.org, 
              mapping = "org.Hs.eg.db",
              ID = "entrez")
GOdata

numSigGenes(GOdata)
# -> GO BP: 63 of the 68 genes in the HZup core gene set are feasible genes -> relevant for the GeneRatio

# Perform overrepresentation analysis using Fisher’s exact test and the elim algorithm.
resultFisher.elim <- runTest(GOdata, 
                             algorithm = "elim", 
                             statistic = "fisher")

sum(score(resultFisher.elim) < 0.05)
# -> 118 GO BP terms with a p-value < 0.05

showSigOfNodes(GOdata, score(resultFisher.elim), 
               firstSigNodes = 15, useInfo = 'def')

fisher.elimRes <- GenTable(GOdata,
                           Fisher.elim = resultFisher.elim,
                           orderBy = "Fisher.elim", 
                           topNodes = 118,
                           numChar=1000)
fisher.elimRes$GeneRatio <- fisher.elimRes$Significant/numSigGenes(GOdata)
# In the fisher.elimRes object, the p-values are stored in the variable 'Fisher.elim' as char and not num
fisher.elimRes$pvalue <- score(resultFisher.elim, whichGO = fisher.elimRes$GO.ID)
fisher.elimRes$Fisher.elim <- NULL

# Remove terms with only 1 or 2 significant genes
fisher.elimRes <- fisher.elimRes[fisher.elimRes$Significant > 2, ]
# -> 70 terms remaining

# Save median lg2FCs and pval of filgo effect as well as gene set size of pathways to file

# Calculate adjusted p-values for filgotinib effect on pathways
# Get unadjusted p-values from gsea F vs. V.
gene_list <- geneListGSEA("F_vs_V.csv")
gseaFV <- gseGO(geneList      = gene_list,
                ont           = "BP",
                OrgDb         = "org.Hs.eg.db",
                keyType       = "ENTREZID",
                minGSSize     = 10,
                maxGSSize     = 20000,
                pvalueCutoff  = 1,
                pAdjustMethod = "BH",
                verbose       = TRUE,
                by            = "fgsea")

pvals <- as.data.frame(gseaFV)[fisher.elimRes$GO.ID, "pvalue"]
# Perform BH correction for the pathways (n=104)
pvals <- p.adjust(pvals, method = "BH")

# Include 'multipleComparisons' function from 'functions.R'
allLg2FCs <- multipleComparisons()
ann.genes <- genesInTerm(GOdata, fisher.elimRes$GO.ID)
medians <- sapply(ann.genes, function(X) { median(allLg2FCs[X, 'F_vs_V']) })

fisher.elimRes$medianFilgoLg2FC <- medians
fisher.elimRes$pvalFilgoEffect <- pvals
write.csv(fisher.elimRes, "Filgotinib effect on GO BP pathways relevant in HZ.csv", row.names = FALSE)

# Number of pathways on which filgotinib has a significant effect
sum(!is.na(fisher.elimRes$pvalFilgoEffect) & fisher.elimRes$pvalFilgoEffect < 0.05)
# -> 28

# Dotplot
library("ggplot2")
library("dplyr")
library("DOSE")
library("stringr")
library("ggpubr")

# Plot the GO terms with the top 15 p-values in the overrepresentation analysis
dataToPlot <- fisher.elimRes[1:15, ]
# Order according to Gene Ratio first, followed by p-value
dataToPlot <- dataToPlot[order(dataToPlot$GeneRatio, -dataToPlot$pvalue, decreasing = FALSE), ]

# Add column containing ordered GOs as a factor
dataToPlot <- dataToPlot %>% 
  mutate(GOs_ordered = factor(row_number(), labels = Term))

min(dataToPlot$pvalue)

plot1 <- ggplot(data = dataToPlot, 
                aes(x = GeneRatio, 
                    y = GOs_ordered, 
                    color = pvalue, 
                    size = Significant)) + 
  xlim(NA, 0.4) +
  geom_point() +
  scale_size(range=c(3, 8), breaks = c(3, 10, 23)) +
  scale_color_gradient(low = "red3", 
                       high = "blue", 
                       transform = 'log10', 
                       limits = c(1e-7, NA), 
                       n.breaks = 4) +
  theme_dose(10) + 
  scale_y_discrete(labels = function(x) str_wrap(dataToPlot$Term, width = 47)) +
  ylab(NULL) + 
  xlab("Gene ratio") + 
  labs (color = "p value", 
        size = "Gene count", 
        title = "GO BP overrepresentation analysis") +
  theme(legend.position="bottom", 
        legend.box = "vertical",
        plot.title = element_text(size = 12))

# Plot effect of filgotinib on the top 15 GO pathways
min(dataToPlot$pvalFilgoEffect)
min(dataToPlot$Annotated)
max(dataToPlot$Annotated)

plot2 <- ggplot(data = dataToPlot, 
                aes(x = medianFilgoLg2FC, 
                    y = GOs_ordered, 
                    color = pvalFilgoEffect, 
                    size = Annotated)) +
    geom_point() +
    scale_size(range=c(3, 8), breaks = c(27, 400, 844)) +
    scale_color_continuous(low = "red3", 
                           high = "blue",
                           guide = guide_colorbar(reverse=FALSE), 
                           transform = 'log10',
                           limits = c(0.0006, 0.05),
                           n.breaks = 4) +
    scale_y_discrete(labels = NULL) +
    labs(size = "Gene set size", 
         color = "p value") +
    guides(color = guide_colourbar(order=1), 
           size = guide_legend(order=2)) +
    ylab(NULL) + 
    xlab("Median filgotinib induced lg2FC") +
    theme_dose(10) +
    ggtitle("Effect of filgotinib on GO term") +
    theme(legend.position="bottom", 
          legend.box = "vertical",
          plot.title = element_text(size = 12))


ggarrange(plotlist = list(plot1, plot2), ncol = 2, nrow = 1, widths = c(2, 1))


# Retrieve GO MF terms
GOdata <- new("topGOdata", 
              description = "HZ core up geneset", 
              ontology = "MF",
              allGenes = geneList,
              nodeSize = 10,
              annot = annFUN.org, 
              mapping = "org.Hs.eg.db",
              ID = "entrez")
GOdata

numSigGenes(GOdata)
# -> GO MF: 65 of the 68 genes in the HZup core gene set are feasible genes
resultFisher.elim <- runTest(GOdata, 
                             algorithm = "elim", 
                             statistic = "fisher")

sum(score(resultFisher.elim) < 0.05)
# -> 9 GO MF terms with a p-value < 0.05

fisher.elimRes <- GenTable(GOdata,
                           Fisher.elim = resultFisher.elim,
                           orderBy = "Fisher.elim", 
                           topNodes = 9,
                           numChar=1000)
fisher.elimRes$GeneRatio <- fisher.elimRes$Significant/numSigGenes(GOdata)
# In the fisher.elimRes object, the p-values are stored in the variable 'Fisher.elim' as char and not num
fisher.elimRes$pvalue <- score(resultFisher.elim, whichGO = fisher.elimRes$GO.ID)
fisher.elimRes$Fisher.elim <- NULL

# Remove terms with only 1 or 2 significant genes
fisher.elimRes <- fisher.elimRes[fisher.elimRes$Significant > 2, ]
# -> 4 terms remaining

# Save median lg2FCs and pval of filgo effect as well as gene set size of pathways to file

# Calculate adjusted p-values for filgotinib effect on pathways
# Get unadjusted p-values from gsea F vs. V.
gseaFV <- gseGO(geneList      = gene_list,
                ont           = "MF",
                OrgDb         = "org.Hs.eg.db",
                keyType       = "ENTREZID",
                minGSSize     = 10,
                maxGSSize     = 20000,
                pvalueCutoff  = 1,
                pAdjustMethod = "BH",
                verbose       = TRUE,
                by            = "fgsea")

pvals <- as.data.frame(gseaFV)[fisher.elimRes$GO.ID, "pvalue"]
# Perform BH correction for the pathways (n=4)
pvals <- p.adjust(pvals, method = "BH")

ann.genes <- genesInTerm(GOdata, fisher.elimRes$GO.ID)
medians <- sapply(ann.genes, function(X) { median(allLg2FCs[X, 'F_vs_V']) })

fisher.elimRes$medianFilgoLg2FC <- medians
fisher.elimRes$pvalFilgoEffect <- pvals
write.csv(fisher.elimRes, "Filgotinib effect on GO MF pathways relevant in HZ.csv", row.names = FALSE)

# Number of pathways on which filgotinib has a significant effect
sum(!is.na(fisher.elimRes$pvalFilgoEffect) & fisher.elimRes$pvalFilgoEffect < 0.05)
# -> 1
