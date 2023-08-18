# Description -------------
# Identify correlations between DNA methylation and gene expression

# Get Packages --------------------------------------------------------------
# Package names
packages <- c("DMRichR", "ComplexHeatmap", "RColorBrewer", "circlize", "ggplot2", "magick", "biomaRt", "ggpubr", "gridExtra")

packages <- c("DMRichR", "dplyr", "ggplot2")
# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Load packages
invisible(lapply(packages, library, character.only = TRUE))

# Pearson functions -------------
# Evaluates pearson correlation and makes correlation scatterplots for top "n" DMRs
pearson_co <- function(methyl = methyl, reads = reads, DMRs = DMRs, project = project, n = n){
  beta <- methyl[,5:length(methyl)] %>% data.matrix() %>% t() # beta values
  counts <- reads[2:length(reads)] %>% data.matrix %>% t()
  co <- data.frame(DMR = numeric(nrow(methyl)), gene = numeric(nrow(methyl)), annot = numeric(nrow(methyl)),
                   coef = numeric(nrow(methyl)), p = numeric(nrow(methyl))) 
  for (i in 1:nrow(methyl)){
    x = cor.test(beta[,i], counts[,i], method = "pearson")
    co$p[i] = x[["p.value"]]
    co$coef[i] = cor(beta[,i], counts[,i], method = "pearson")
  }
  DMRname <- c(paste(methyl$chr, methyl$start, methyl$end, sep = ";"))
  co$DMR <- DMRname
  co$gene <- methyl$geneSymbol
  co$annot <- DMRs$annotation
  co$Direction <- DMRs$direction
  OrderedCo <- co[order(co$p),] # reorder by p value
  write.table(OrderedCo, file = glue::glue("DMR_DEG/{project}_pearson.txt"), 
              quote = FALSE, sep = '\t ', row.names = F) # save pearson 
  
  # Scatter plots
  plot_list = list()
  for (k in 1:n){
    DMR = OrderedCo[k,1]
    data = data.frame(meth = beta[,which(co$DMR == DMR)], reads = counts[,which(co$DMR == DMR)])
    p = ggpubr::ggscatter(data, x = "meth", y = "reads",
                          add = "reg.line", conf.int = TRUE, 
                          cor.coef = TRUE, cor.method = "pearson",
                          xlab = "DNA methylation", ylab = "Gene Expression",
                          title = paste(DMR, co[co$DMR == DMR, 2], sep = "-"))+
      theme(plot.title = element_text(size=12))
    plot_list[[k]] = p
  }
  
  ggsave(filename = glue::glue("DMR_DEG/{project}_TopPearson.pdf"),
         plot = gridExtra::marrangeGrob(plot_list, nrow=3, ncol=2), width = 9, height = 9) # save scatterplots
  sigco <- OrderedCo %>% plyranges::filter(p < .05)
}

# Determines percent significant correlation for each genic region
pcor <- function(DMRs = DMRs, sigco = sigco, project = project){
  DMRs$annotation <- factor(DMRs$annotation, levels=c("Promoter", "5' UTR", "Exon", "Intron", "3' UTR", "Downstream", "Distal Intergenic")) # reorder
  annots <- levels(as.factor(DMRs$annotation)) # annotations
  sigco_counts = data.frame(Annotation = c(annots, annots),
                            Correlation = rep(c("positive", "negative"), each=length(annots)), 
                            Count = numeric(length(annots)*2),
                            Percent = numeric(length(annots)*2)) # initialize dataframe
  for (j in 1:length(levels(as.factor(DMRs$annotation)))){
    annot <- levels(as.factor(DMRs$annotation))[j]
    sigco_counts[j,3] <- nrow(sigco[which(sigco$annot == annot & sigco$coef > 0),]) # positive correlation
    sigco_counts[j+length(annots),3] <- nrow(sigco[which(sigco$annot == annot & sigco$coef < 0),]) # negative correlation
  } # get counts for each annotation
  pos <- sum(sigco_counts[which(sigco_counts$Correlation == "positive"), 3])
  neg <- sum(sigco_counts[which(sigco_counts$Correlation == "negative"), 3])
  sigco_counts$Percent[1:length(annots)] <- sigco_counts$Count[1:length(annots)]/pos # positive correlation percentage
  sigco_counts$Percent[(length(annots)+1):(length(annots)*2)] <- sigco_counts$Count[(length(annots)+1):(length(annots)*2)]/neg # negative correlation percentage
  sigco_counts$Percent <- round(sigco_counts$Percent, 2)
  write.table(sigco_counts, file = glue::glue("DMR_DEG/{project}_cor_counts.txt"), quote = FALSE, sep = '\t ', row.names = F) # saves the results as a text file in the working directory
  
  # Plot
  sigco_counts$Annotation = factor(sigco_counts$Annotation, levels=c("Promoter", "5' UTR", "Exon", "Intron", "3' UTR", "Downstream", "Distal Intergenic"))
  plot <- ggplot2::ggplot(sigco_counts, aes(fill=Annotation, y=Percent, x=Correlation)) + 
    geom_bar(position="dodge", stat = "identity", color = "black")+
    scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
    labs(x ="", y = "") + ggtitle(project) +
    scale_x_discrete(labels=c("Negative Correlation", "Positive Correlation")) +
    theme_minimal() + 
    scale_fill_manual(values=wesanderson::wes_palette("Zissou1", n = 7, type = "continuous") %>% rev())
  ggsave(filename = glue::glue("DMR_DEG/{project}_sigcorr.pdf"),
         plot = plot, width = 6, height = 4) # save scatterplots
}
#### DMR correlation tests for mT vs mM #### -------------

# Get TPMs and smoothened methylation values
# TPM table
TM_countData <- read.table("DMR_DEG/mTM_TPM.txt",header=T, sep = "\t") 
rownames(TM_countData) <- TM_countData$Gene
TM_countData <- TM_countData[,2:length(TM_countData)]

# Methylation table q < .01
#DMRichR::annotationDatabases(genome = "mm9", EnsDb = FALSE)
#TM_DMRs <- read.table("mNPTM_DMRs/TM_DMR_individual_smoothed_methylation.txt",header=T, sep = "\t") %>%
#  makeGRangesFromDataFrame(., keep.extra.columns = T) %>%
#  DMRichR::annotateRegions(TxDb = TxDb, annoDb = annoDb) # 4116 DMRs
#TM_DMRs <- TM_DMRs[which(TM_DMRs$geneSymbol %in% rownames(TM_countData)),] # 3456 regions
#TM_meth <- TM_DMRs[,c(1,2,3,32,12:19)]

# Methylation table p < .05
DMRichR::annotationDatabases(genome = "mm9", EnsDb = FALSE)
TM_DMRs <- read.table("mNPTM_DMRs/TM_background_region_individual_smoothed_methylation.txt",header=T, sep = "\t") %>%
  makeGRangesFromDataFrame(., keep.extra.columns = T) %>%
  plyranges::filter(pval < 0.05) %>%
  DMRichR::annotateRegions(TxDb = TxDb, annoDb = annoDb) # 43400 DMRs
TM_DMRs <- TM_DMRs[which(TM_DMRs$geneSymbol %in% rownames(TM_countData)),] # 36300 regions
TM_DMRs <- na.omit(TM_DMRs) # 36287 regions
TM_meth <- TM_DMRs[,c(1,2,3,32,12:19)]

# TPM table for DMR genes
TM_reads <- data.frame(gene = TM_meth$geneSymbol, matrix(0,length(TM_meth$geneSymbol), ncol(TM_countData)))
for (i in 1:length(TM_meth$geneSymbol)){
  TM_reads[i,2:length(TM_reads)] = TM_countData[which(rownames(TM_countData) %in% TM_meth$geneSymbol[i]),1:ncol(TM_countData)]
}
colnames(TM_reads) = c("gene", colnames(TM_countData))

TM_sigco <- pearson_co(methyl = TM_meth, reads = TM_reads, DMRs = TM_DMRs, project = "TM", n = 30) # 8394
pcor(DMRs = TM_DMRs, sigco = TM_sigco, project = "TM")

# DMR correlation tests for early vs late -------------

# DESeq2 normalized read count table
EL_countData = read.table("DMR_DEG/PDO_readcount.txt",header=T, sep = "\t") # raw read counts
rownames(EL_countData) = EL_countData$Gene
EL_countData = EL_countData[,2:length(EL_countData)]
ELmetaData = read.table("DMR_DEG/PDO_metadata.txt",header=T, sep = "\t")
row.names(ELmetaData) = ELmetaData[,1]
ELmetaData = ELmetaData[2]
norm_read <- DESeq2::DESeqDataSetFromMatrix(countData = EL_countData, colData = ELmetaData, design = ~CombinedStage) %>%
              DESeq2::estimateSizeFactors() %>% DESeq2::counts(., normalized=TRUE) %>% as.data.frame()
write.table(norm_read, file = "DMR_DEG/PDO_normRead.txt", quote = FALSE, sep = '\t ', row.names = T) 
#norm_read <- read.table("DMR_DEG/PDO_normRead.txt", header = T, sep = '\t')

# Methylation table q < .05
#DMRichR::annotationDatabases(genome = "hg38", EnsDb = FALSE)
#EL_DMRs <- read.table("PDO_DMRs/EL_DMR_individual_smoothed_methylation.txt",header=T, sep = "\t") %>%
#          makeGRangesFromDataFrame(., keep.extra.columns = T) %>% 
#          DMRichR::annotateRegions(TxDb = TxDb, annoDb = annoDb) # 5374 DMRs
#EL_DMRs <- EL_DMRs[which(EL_DMRs$geneSymbol %in% rownames(norm_read)),] # 3940 DMRs
#EL_meth <- EL_DMRs[,c(1,2,3,51,12:42)]

# Methylation table p < .05
DMRichR::annotationDatabases(genome = "hg38", EnsDb = FALSE)
EL_DMRs <- read.table("PDO_DMRs/EL_background_region_individual_smoothed_methylation.txt",header=T, sep = "\t") %>%
  makeGRangesFromDataFrame(., keep.extra.columns = T) %>%
  plyranges::filter(pval < 0.05) %>%
  DMRichR::annotateRegions(TxDb = TxDb, annoDb = annoDb) # 48220 DMRs
EL_DMRs <- EL_DMRs[which(EL_DMRs$geneSymbol %in% rownames(norm_read)),] # 34868 DMRs
EL_meth <- EL_DMRs[,c(1,2,3,51,12:42)]

# Normalized read table for DMR genes
EL_reads = data.frame(gene = EL_meth$geneSymbol, matrix(0,length(EL_meth$geneSymbol), ncol(norm_read)))
for (i in 1:length(EL_meth$geneSymbol)){
  EL_reads[i,2:length(EL_reads)] = norm_read[which(rownames(norm_read) %in% EL_meth$geneSymbol[i]),1:ncol(norm_read)]
}
colnames(EL_reads) = c("gene", colnames(norm_read))

EL_sigco <- pearson_co(methyl = EL_meth, reads = EL_reads, DMRs = EL_DMRs, project = "EL", n = 30) # 3406
pcor(DMRs = EL_DMRs, sigco = EL_sigco, project = "EL")

# Correlation graph function -----
corrGene <- function(gene = gene, index = index, sigco = sigco, reads = reads, metadata = metadata, methyl = methyl, org = org){
  if (org == "human"){
    count <- data.frame(t(reads[which(row.names(reads) == gene),]), metadata$CombinedStage) # gene expression
    colnames(count) <- c("count", "stage")
    x <- which(methyl$geneSymbol == gene)[index]
    meth <- data.frame(meth = t(methyl[x,5:length(methyl)]), stage = metadata$CombinedStage) # methylation
    data <- data.frame(meth = meth$meth, reads = count$count)
  } else if(org == "mouse"){
    count <- data.frame(t(reads[which(row.names(reads) == gene),]), metadata$Stage) # gene expression
    colnames(count) <- c("count", "stage")
    x <- which(methyl$geneSymbol == gene)[index]
    meth <- data.frame(meth = t(methyl[x,5:length(methyl)]), stage = metadata$Stage) # methylation
    data <- data.frame(meth = meth$meth, reads = count$count)
  } else {
    print("Genome not compatible. Please choose mouse or human")
  }
  data <- data.frame(meth = meth$meth, reads = count$count)
  DMRname <- c(paste(methyl$chr[x], methyl$start[x], methyl$end[x], sep = ";"))
  sc <- ggpubr::ggscatter(data, x = "meth", y = "reads",
                          add = "reg.line", conf.int = TRUE, 
                          cor.coef = TRUE, cor.method = "pearson",
                          xlab = "DNA methylation", ylab = "Gene Expression",
                          title=paste(DMRname,gene,sep = " - "))
  rna <- ggplot2::ggplot(count, aes(x=stage, y=count, fill=stage)) + 
    geom_violin(trim=FALSE)+
    geom_boxplot(width=0.1, fill="white")+
    labs(title=gene,x="", y = "Gene Expression") +
    scale_fill_manual(values=c("blue", "red")) +
    theme_classic() + 
    theme(legend.position = "none")
  dna <- ggplot2::ggplot(meth, aes(x=stage, y=meth, fill=stage)) + 
    geom_violin(trim=FALSE)+
    geom_boxplot(width=0.1, fill="white")+
    labs(title=DMRname,x="", y = "Methylation") +
    scale_fill_manual(values=c("blue", "red")) +
    theme_classic() + 
    theme(legend.position = "none")
  cowplot::plot_grid(dna, sc, align = "h", nrow = 1, rel_widths = c(1, 1.5))
  ggsave(filename = glue::glue("DMR_DEG/{gene}_{index}.pdf"), width = 6, height = 3)
}

# Plot genes -------------
EL_sigGenes = unique(EL_sigco$gene) # 744 genes
TM_sigGenes = unique(TM_sigco$gene) # 2619 genes
int = intersect(EL_sigGenes, toupper(TM_sigGenes)) # 165 genes
print(int)

TM_prot = TM_sigco %>% plyranges::filter(annot == "Promoter") # 1047 genes
EL_prot = EL_sigco %>% plyranges::filter(annot == "Promoter") # 230 genes
int_prot = intersect(EL_prot$gene, toupper(TM_prot$gene)) # 13 genes

TMmetaData = read.table("DMR_DEG/TM_metadata.txt",header=T, sep = "\t")
row.names(TMmetaData) = TMmetaData[,1]
TMmetaData = TMmetaData[2]

hgene <- "PLAT"
hDMRs = EL_DMRs[which(EL_DMRs$geneSymbol == hgene),] # list of DMRs associated with gene
index = 1 # index of DMR to plot
corrGene(gene = hgene, index = index, sigco = EL_sigco, reads = norm_read, metadata = ELmetaData, methyl = EL_meth, org = "human")
file.rename("DMR_DEG/PLAT_1.pdf", "DMR_DEG/hPLAT_1.pdf")

mgene <- "Plat"
mDMRs = TM_DMRs[which(TM_DMRs$geneSymbol == mgene),] # list of DMRs associated with gene
index = 2
corrGene(gene = mgene, index = index, sigco = TM_sigco, reads = TM_countData, metadata = TMmetaData, methyl = TM_meth, org = "mouse")
file.rename("DMR_DEG/Plat_2.pdf", "DMR_DEG/mPlat_2.pdf")

hgene <- "LY6D"
hDMRs = EL_DMRs[which(EL_DMRs$geneSymbol == hgene),] # list of DMRs associated with gene
index = 1 # index of DMR to plot
corrGene(gene = hgene, index = index, sigco = EL_sigco, reads = norm_read, metadata = ELmetaData, methyl = EL_meth, org = "human")
file.rename("DMR_DEG/LY6D_1.pdf", "DMR_DEG/hLY6D_1.pdf")

mgene <- "Ly6d"
mDMRs = TM_DMRs[which(TM_DMRs$geneSymbol == mgene),] # list of DMRs associated with gene
index = 1
corrGene(gene = mgene, index = index, sigco = TM_sigco, reads = TM_countData, metadata = TMmetaData, methyl = TM_meth, org = "mouse")
file.rename("DMR_DEG/Ly6d_1.pdf", "DMR_DEG/mLy6d_1.pdf")