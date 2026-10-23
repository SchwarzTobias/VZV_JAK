# Annotate DESeqResult function
library("org.Hs.eg.db")

annotateGene <- function(DESeqResult, keytype = "SYMBOL", returnLists = FALSE) {
  # Lists of GO and KEGG ids cannot be easily saved to file, therefore are not returned by default
  DESeqResult$genetype <- unname(mapIds(org.Hs.eg.db,
                                        keys=as.character(rownames(DESeqResult)),
                                        column="GENETYPE",
                                        keytype=keytype,
                                        multiVals="first"))
  DESeqResult$genename <- unname(mapIds(org.Hs.eg.db,
                                        keys=as.character(rownames(DESeqResult)),
                                        column="GENENAME",
                                        keytype=keytype,
                                        multiVals="first"))
  DESeqResult$gene_id <- unname(mapIds(org.Hs.eg.db,
                                       keys=as.character(rownames(DESeqResult)),
                                       column="ENTREZID",
                                       keytype=keytype,
                                       multiVals="first"))
  DESeqResult$symbol <-  unname(mapIds(org.Hs.eg.db,
                                       keys=as.character(rownames(DESeqResult)),
                                       column="SYMBOL",
                                       keytype=keytype,
                                       multiVals="first"))
  if (returnLists == TRUE) {
    # Gene ontology without child terms
    DESeqResult$go <- mapIds(org.Hs.eg.db,
                             keys=as.character(rownames(DESeqResult)),
                             column="GO",
                             keytype=keytype,
                             multiVals="filter")
    # Gene ontology with child terms
    DESeqResult$go_all <- mapIds(org.Hs.eg.db,
                                 keys=as.character(rownames(DESeqResult)),
                                 column="GOALL",
                                 keytype=keytype,
                                 multiVals="filter")
    # KEGG pathway
    DESeqResult$kegg <- mapIds(org.Hs.eg.db,
                               keys=as.character(rownames(DESeqResult)),
                               column="PATH",
                               keytype=keytype,
                               multiVals="filter")
  }
  DESeqResult
}

# Create Volcano Plot
library(EnhancedVolcano)

plotVolcano <- function(res, title = "", pCutoff = 0.05, FCcutoff = 1) {
  # Lg2FC with padj != NA
  Lg2FCs <- res[!is.na(res$padj), "log2FoldChange"]
  minLg2FC <- min(Lg2FCs)
  print(paste("Min Lg2FC:", minLg2FC))
  maxLg2FC <- max(Lg2FCs)
  print(paste("Max Lg2FC:", maxLg2FC))
  minp <- min(na.omit(res$padj))
  print(paste("Min p.adj", minp))
  
  EnhancedVolcano(res,
                  lab = res$symbol,
                  title = title,
                  subtitle = NULL,
                  x = "log2FoldChange",
                  y = "padj",
                  xlim = c(minLg2FC, maxLg2FC),
                  ylim = c(0, -log10(minp)),
                  xlab = bquote(~Log[2]~Fold~Change),
                  ylab = bquote(~-Log[10]~P~Value),
                  pCutoff = pCutoff,
                  FCcutoff = FCcutoff,
                  col = c("grey30", "forestgreen", "royalblue", "red2"),
                  legendPosition = "none",
                  pointSize = 1.0,
                  labSize = 3.0,
                  drawConnectors = TRUE,
                  widthConnectors = 0.5)
}

# Function to remove uncharacterized ncRNAs and pseudogenes from DESeq2 results
library("purrr")

removeNoise <- function(DESeqResult, pseudogenes = TRUE) {
  # Remove uncharacterized ncRNAs
  remove <- which(DESeqResult$genetype == "ncRNA" & grepl('uncharacterized', DESeqResult$genename, fixed=TRUE) == TRUE)
  # Remove genes with missing annotation
  remove <- c(remove, which(is.na(DESeqResult$genetype)))
  # Remove pseudogenes
  if (pseudogenes == TRUE) {
    remove <- c(remove, which(DESeqResult$genetype == "pseudo"))
  }
  # Return DESeqResult without the genes selected for removal
  if (is_empty(remove) == TRUE) {
    DESeqResult
  } else {
    DESeqResult[-remove, ]
  }
}

# Function to load lg2FCs and padj-values from DESeq2 results saved to file
geneList <- function(filename) {
  # reading in data from Deseq2 result
  df = read.csv(filename, header=TRUE)
  
  # Remove entries with no gene symbol
  if (sum(is.na(df$symbol) | df$symbol == "") > 0) {
    print(paste("Removing entries with no gene symbol:", sum(is.na(df$symbol) | df$symbol == "")))
    df <- df[!is.na(df$symbol) | df$symbol != "", ]
  }
  
  # Remove entries with duplicate gene symbol
  if (sum(duplicated(df$symbol)) > 0) {
    print(paste("Removing entries with duplicate symbol:", sum(duplicated(df$symbol))))
    df <- df[!duplicated(df$symbol), ]
  }
  
  # We need the lg2FC and the padj-value
  gene_list <- data.frame(
    lg2FC = df$log2FoldChange,
    padj  = df$padj)
  
  # name the vector
  rownames(gene_list) <- df$symbol
  
  gene_list
}

# Create gene list for GSEA from file
geneListGSEA <- function(filename) {
  # reading in data from DESeq2
  df = read.csv(filename, header=TRUE)

  # Remove entries with no gene ID
  if (sum(is.na(df$gene_id) | df$gene_id == "") > 0) {
    print(paste("Removing entries with no gene ID:", sum(is.na(df$gene_id) | df$gene_id == "")))
    df <- df[!is.na(df$gene_id) | df$gene_id != "", ]
  }

  # Remove entries with duplicate gene ID
  if (sum(duplicated(df$gene_id)) > 0) {
    print(paste("Removing entries with duplicate gene ID:", sum(duplicated(df$gene_id))))
    df <- df[!duplicated(df$gene_id), ]
  }

  # We need the lg2FC
  gene_list <- df$log2FoldChange
  
  # name the vector
  names(gene_list) <- df[, "gene_id"]
  
  # omit NA values (such values should not exist in DESeq output)
  gene_list <- na.omit(gene_list)
  
  # sort the list in decreasing order (required for clusterProfiler)
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  gene_list
}

# Create matrix with lg2FCs of the 3 comparisons (V vs. O, F vs. V, F vs. O)
multipleComparisons <- function() {
  # Load gene lists and order by name
  V_vs_O <- geneListGSEA("V_vs_O.csv")
  V_vs_O <- V_vs_O[order(names(V_vs_O))]
  F_vs_V <- geneListGSEA("F_vs_V.csv")
  F_vs_V <- F_vs_V[order(names(F_vs_V))]
  F_vs_O <- geneListGSEA("F_vs_O.csv")
  F_vs_O <- F_vs_O[order(names(F_vs_O))]
  # Check that genes in lists are identical
  if (sum(names(V_vs_O) != names(F_vs_V)))  stop("Error: Gene lists are not identical!")
  if (sum(names(V_vs_O) != names(F_vs_O)))  stop("Error: Gene lists are not identical!")
  # Return matrix with lg2FCs
  cbind(V_vs_O, F_vs_V, F_vs_O)
}

# Plot ridgeplot
library("ggplot2")
library("DOSE")

myridgeplot <- function(coreLg2FCs, pvals, title, subtitle, vline = NULL) {
  if (!exists("minpval"))
    minpval <- min(na.omit(pvals))
  
  # nn: Vector containing the names of the different comparisons
  nn <- c("V_vs_O", "F_vs_V", "F_vs_O")
  
  # len: length of the complete core set
  len <- nrow(coreLg2FCs)
  
  # gs2val.df: colnum: sum of the genes in all gene sets (in my case 3 * length(complete_core))
  # $category: name of the gene set (in my case 'V_vs_O', 'F_vs_V', 'F_vs_O') (repeated)
  # $color: pvalue (repeated)
  # $value: lg2fcs
  gs2val.df <- data.frame(category = rep(nn, each=len),
                          p.adjust = rep(pvals, each=len),
                          value = as.vector(coreLg2FCs))
  
  # Make factor out of category
  gs2val.df$category <- factor(gs2val.df$category, levels=c("F_vs_O", "F_vs_V", "V_vs_O"))
  
  plot <- ggplot(gs2val.df, aes_string(x="value", y="category", fill="p.adjust")) +
    ggridges::geom_density_ridges() +
    scale_fill_continuous(low="red",
                          high="blue",
                          name="p.adjust",
                          guide = guide_colorbar(),
                          limits = c(minpval, 0.05),
                          transform = "log10") +
    scale_y_discrete(labels = c("VZV+Filgo vs. O", "VZV+Filgo vs. VZV", "VZV vs. O")) +
    geom_vline(xintercept=0, 
               color='grey70', 
               linetype='dashed') +
    xlab(paste("lg2FC expression (n = ", len, ")", sep = "")) + 
    ylab(NULL) + 
    labs(title = title, subtitle = subtitle) + 
    theme_dose(font.size = 10)
  
  if (!is.null(vline))
    plot <- plot + geom_vline(xintercept = vline,
                              color = 'grey8',
                              linetype = 'dotted')
  
  plot
}
