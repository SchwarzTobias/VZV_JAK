# Analysis of differentially expressed genes
library(DESeq2)

# Load data
countdata <- as.matrix(read.csv("countdata.csv", row.names=1))
coldata <- read.csv("coldata.csv")
coldata$sample <- as.factor(coldata$sample)
coldata$stimulation <- as.factor(coldata$stimulation)
# Relevel factor stimulation, so the untreated samples are the reference level
coldata$stimulation <- factor(coldata$stimulation, levels = c("O", "V", "F"))

# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~ sample + stimulation)

# Prefiltering the dataset: At least 2/3 of the samples need to have a read count >5
smallestGroupSize <- 2/3*nrow(coldata)
keep <- rowSums(counts(dds) >= 5) >= smallestGroupSize
dds <- dds[keep,]
nrow(dds)
# -> 17842 genes to analyze

dds <- DESeq(dds)

# Include 'annotateGene' function from 'functions.R'

# Results VZV vs. O
resVO <- results(dds, contrast = c("stimulation", "V", "O"), independentFiltering = FALSE, alpha = 0.01)
resVO <- annotateGene(resVO)
write.csv(as.data.frame(resVO[order(resVO$log2FoldChange, decreasing = TRUE),]), file = "V_vs_O.csv")
# Genes significantly altered by VZV antigen (V vs O)
resVOsig <- resVO[!is.na(resVO$padj) & resVO$padj < 0.01 & abs(resVO$log2FoldChange) > 1, ]
summary(resVOsig, alpha = 0.01)
# 1168 up, 465 down out of 17842

# Results VZV+Filgo vs. VZV
resFV <- results(dds, contrast = c("stimulation", "F", "V"), independentFiltering = FALSE, alpha = 0.01)
resFV <- annotateGene(resFV)
write.csv(as.data.frame(resFV[order(resFV$log2FoldChange, decreasing = TRUE),]), file = "F_vs_V.csv")
# Genes significantly altered by Filgotinib (F vs V)
resFVsig <- resFV[!is.na(resFV$padj) & resFV$padj < 0.01 & abs(resFV$log2FoldChange) > 1, ]
summary(resFVsig, alpha = 0.01)
# 360 up, 1076 down out of 17842

# Results VZV+Filgo vs. O
resFO <- results(dds, contrast = c("stimulation", "F", "O"), independentFiltering = FALSE, alpha = 0.01)
resFO <- annotateGene(resFO)
write.csv(as.data.frame(resFO[order(resFO$log2FoldChange, decreasing = TRUE),]), file = "F_vs_O.csv")
resFOsig <- resFO[!is.na(resFO$padj) & resFO$padj < 0.01 & abs(resFO$log2FoldChange) > 1, ]
summary(resFOsig, alpha = 0.01)
# 1050 up, 1694 down out of 17842


# Create Volcano Plots
# Include 'plotVolcano' function from 'functions.R'
# Include function 'removeNoise' from file 'functions.R'
# Figure S1a
plotVolcano(removeNoise(resVO), title = "Effects of VZV antigen stimulation on rPBMC", pCutoff = 0.01, FCcutoff = 1)
# Figure S1b
plotVolcano(removeNoise(resFV), title = "Effects of filgotinib on VZV antigen stimulated rPBMC", pCutoff = 0.01, FCcutoff = 1)


# Effects of Filgotinib on genes significantly altered by VZV antigen
resFVinVO <- resFV[rownames(resFV) %in% rownames(resVOsig), ]

# Genes significantly altered by Filgotinib out of the genes significantly altered by VZV antigen
resFVinVOsig <- resFVinVO[resFVinVO$padj < 0.01, ]
summary(resFVinVOsig, alpha = 0.01)
# Out of the 1633 genes significantly altered by VZV antigen, 126 are upregulated and 502 are downregulated
# following addition of filgotinib (628 total, 38.5 %).

# Genes significantly altered by filgotinib out of the genes significantly altered by VZV antigen
# restricted to genes where the effects of filgotinib and VZV antigen have differing directions (counterregulation).
diffExp <- sign(resVO$log2FoldChange) != sign(resFV$log2FoldChange)
diffExp <- rownames(resVO[diffExp, ])
resCounterregulated <- resFVinVOsig[rownames(resFVinVOsig) %in% diffExp, ]
summary(resCounterregulated, alpha = 0.01)
# Out of the 628 genes significantly altered by VZV antigen stimulation and by addition of filgotinib,
# 565 are counterregulated (107 up, 458 down).

# Genes significantly up- or downregulated by VZV antigen compared to O,
# and following addition of filgotinib, are significantly up- or downregulated even more (amplification).
diffExp <- sign(resVO$log2FoldChange) == sign(resFV$log2FoldChange)
diffExp <- rownames(resVO[diffExp, ])
resAmplified <- resFVinVOsig[rownames(resFVinVOsig) %in% diffExp, ]
summary(resAmplified, alpha = 0.01)
# Out of the 628 genes significantly altered by VZV antigen stimulation and by addition of filgotinib,
# 63 are amplified (19 up, 44 down, 3.9% of the 1633 genes).

# Genes significantly altered by VZV antigen, but significant difference is abrogated upon addition of filgotinib:
# Results VZV+Filgo vs. O
resFO <- results(dds, contrast = c("stimulation", "F", "O"), independentFiltering = FALSE, alpha = 0.01)
resFO <- annotateGene(resFO)
# Genes not significantly different between F and O
resFOnosig <- resFO[is.na(resFO$padj) | resFO$padj > 0.01 | abs(resFO$log2FoldChange) < 1, ]
# Genes not significantly different between F and O 
# out of the genes significantly altered by VZV antigen
resnoFOinVO <- resFOnosig[rownames(resFOnosig) %in% rownames(resVOsig), ]
summary(resnoFOinVO, alpha = 0.01)
# -> Filgotinib abrogates the significant effect of VZV antigen stimulation in 766 out of the 1633 genes (46.9%).

# Genes significantly up- or downregulated by VZV antigen compared to O,
# but following addition of filgotinib, are significantly down- or upregulated compared to O (reversal).
# The lg2FC of V vs. O and the lg2FC of F vs. O need to have differing signs
diffExp <- sign(resVO$log2FoldChange) != sign(resFO$log2FoldChange)
resReversal <- resFO[(resFO$padj < 0.01)
                      & (rownames(resFO) %in% rownames(resVOsig))
                      & (diffExp == TRUE), ]
summary(resReversal)
# Filgotinib reverses the significant effect of VZV antigen stimulation in 19 out of the 1633 genes (1.2%).
# (5 up, 14 down)


# Plot the effects of VZV antigen stimulation vs. the effect of VZV antigen stimulation in the presence of filgotinib.

# Plot all genes excluding ncRNAs, pseudogenes, etc.
showGenes <- rownames(removeNoise(resVO))

# -> 1429 out of the 1633 genes remaining

FC_VO <- resVO[rownames(resVO) %in% showGenes, ]$log2FoldChange
FC_FO <- resFO[rownames(resFO) %in% showGenes, ]$log2FoldChange

# Create matrix with lg2FCs
ma <- cbind(FC_VO, FC_FO)
rownames(ma) <- showGenes

# Color dots according to the effects of filgotinib on VZV activated genes
# No significant effect of VZV antigen -> color 0
# Significant effect of VZV antigen -> color 1
dot.color <- as.integer(showGenes %in% rownames(resVOsig))
# Genes significantly altered by VZV antigen and by filgotinib -> color 2
dot.color2 <- as.integer(showGenes %in% rownames(resFVinVOsig))
dot.color <- dot.color + dot.color2
# Genes were the VZV effect is abrogated by addition of filgotinib -> color 3
dot.color2 <- as.integer((showGenes %in% rownames(resnoFOinVO)*3))
# Cave: filgotinib could abrogate the effect of VZV antigen, even when it has no significant effect.
dot.color <- pmax(dot.color, dot.color2)
# Genes were the VZV effect is reversed by addition of filgotinib -> color 4
dot.color2 <- as.integer((showGenes %in% rownames(resReversal)*4))
dot.color <- pmax(dot.color, dot.color2)
# Genes significantly amplified -> color 5
dot.color2 <- as.integer((showGenes %in% rownames(resAmplified)*5))
dot.color <- pmax(dot.color, dot.color2)

dot.color <- as.factor(dot.color)
levels(dot.color) <- c("No effect of VZV antigen", "No effect of filgotinib", "Significant counterregulation", "Abrogation of VZV effect",
                       "Reversal of VZV effect", "Amplification of VZV effect")

# Make dots with lower V vs O lg2FC more transparent (starting from an abs(lg2FC) of 3)
alpha <- abs(ma[, 1])
alpha[alpha > 3] <- 3
alpha <- alpha/max(alpha)

# Do not show names of genes not altered by VZV antigen
showGenes[!(showGenes %in% rownames(resVOsig))] <- NA

library("ggrepel")
# Limits of the plot
paste(min(ma[, 1]), max(ma[, 1]), min(ma[, 2]), max(ma[, 2]))

ggplot(ma, 
       aes(x = ma[, 1], y = ma[, 2])
) +
  scale_x_continuous("Effects of VZV antigen stimulation on rPBMC (lg2FC)", 
                     limits = c(-5, 8),
                     breaks = c(-4, -2, 0, 2, 4, 6, 8),
                     expand = c(0, 0.2)
  ) +
  scale_y_continuous("Effects of VZV antigen stimulation on rPBMC in the presence of filgotinib (lg2FC)", 
                     limits = c(-5, 6), 
                     breaks = c(-4, -2, 0, 2, 4, 6),
                     expand = c(0, 0.2)
  ) +
  theme(text = element_text(size=12)
  ) +
  geom_point(
    aes(color = dot.color),
    size = 1, 
    alpha = alpha
  ) + 
  scale_color_manual(
    # Use custom colors
    values = c("grey70","grey30", "green3", "blue2", "red2", "orange2")
  ) + 
  labs (
    color = "Filgotinib effect"
  ) + 
  geom_text_repel(
    aes(label = showGenes),
    size = 2.5,
    min.segment.length = 0.5, 
    segment.color = "grey60",
    seed = 42, 
    box.padding = 0.3,
    max.overlaps = 20, 
    direction = "both"
  ) +
  coord_fixed(ratio = 1)

# Print filgotinib-induced relative reduction and p-value for selected genes

genes_of_interest <- c(# Signal transducers
                      "IFNA", "IFNB", "STING1", "MAVS", "IRF3", "STAT2", "STAT3", "TYK2",
                      # Interferon stimulated genes
                      "EIF2AK2", "ISG15", "OAS1", "MX1", "CXCL10", "IFI44", "IFI44L", "IFI27", "RSAD2", "SERPING1", "SIGLEC1", "IFIT5", "OASL", "IFI35", "IFITM1", "IRF1", "IRF9",
                      # Immune system effectors
                      "IL6", "IL12B", "IFNG", "TNF", "IL1B")

idx <- which(rownames(resFV) %in% genes_of_interest)
gene_data <- data.frame(resFV[idx, ]$symbol, 
                        resFV[idx, ]$log2FoldChange,
                        # Median decrease in gene expression in %:
                        (1-2^resFV[idx, ]$log2FoldChange)*100,
                        resVO[idx, ]$padj,
                        resFV[idx, ]$padj)

colnames(gene_data) <- c("symbol", "lg2FC", "rel_red", "padj_VO", "padj_FV")

# Remove genes were change was not significant (p > 0.05)
gene_data <- gene_data[gene_data$padj_FV < 0.05, ]
# -> 20 genes remaining

# Order by lg2FC
gene_data <- gene_data[order(gene_data$lg2FC, decreasing = TRUE),]

print(gene_data)


# Plot log transformed expression of genes following rlog transformation
genes_of_interest <- c("IRF1", "IFIT5", "MX1", "TNF")

rld <- rlog(dds, blind = FALSE)
tmp <- assay(rld)
# Normalize for baseline expression
normalized <- cbind(tmp[, 2]-tmp[, 1], tmp[, 3]-tmp[, 1],
                    tmp[, 5]-tmp[, 4], tmp[, 6]-tmp[, 4],
                    tmp[, 8]-tmp[, 7], tmp[, 9]-tmp[, 7],
                    tmp[, 11]-tmp[, 10], tmp[, 12]-tmp[, 10])

coldata_norm <- data.frame(
  id = c("Ca_V", "Ca_F", "Fe_V", "Fe_F", "G_V", "G_F", "S_V", "S_F"),
  stimulation = as.factor(c("V", "F", "V", "F", "V", "F", "V", "F")),
  sample = as.factor(c("Ca", "Ca", "Fe", "Fe", "G", "G", "S", "S"))
)
coldata_norm$stimulation <- factor(coldata_norm$stimulation, levels = c("V", "F"))

library("scales")
library("ggpubr")

plotGrouped <- function(gene_of_interest) {
  data <- coldata_norm
  data$expression <- normalized[gene_of_interest, ]
  pval <- resFV$padj[rownames(resFV) == gene_of_interest]
  pval <- ifelse(pval <= 0.05, pvalue(pval, accuracy = 0.001, add_p = TRUE), "ns")
  pframe <- data.frame(group1 = "V", group2 = "F", p = pval, y.position = max(data$expression) * 1.1)
  
  ggplot(data, 
         aes(x = stimulation, 
             y = expression, 
             color = sample, 
             group = sample)) +
    geom_point(size = 2) + 
    geom_line() + 
    ggtitle(gene_of_interest) +
    theme(plot.title = element_text(size = 10, face = "bold")) +
    scale_y_continuous("lg2FC", 
                       expand = c(0, max(data$expression) * 0.35 + 0.3)) + 
    scale_x_discrete(name = NULL,
                     labels = c("VZV", "VZV+Fil"),
                     expand = c(0, 0.4)) +
    guides(color = "none") +
    stat_pvalue_manual(pframe, 
                       label = "p",
                       bracket.nudge.y = 0.15,
                       vjust = -0.5)
}

plots <- lapply(genes_of_interest, plotGrouped)
ggarrange(plotlist = plots, ncol = 4, nrow = 1)
