# Description
# DMR calling for mouse organoids + CpG/Genic Enrichment and machine learning for minimal DMR set
# requires bsseq objects from global.R as input

# Get Packages --------------------------------------------------------------
# Package names
packages <- c("ggplot2", "DMRichR", "eulerr", "RColorBrewer")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Load packages
invisible(lapply(packages, library, character.only = TRUE))

# Heatmap function --------------------
heatmap <- function (bs.filtered.bsseq = bs.filtered.bsseq, sigRegions = sigRegions, testCovariate = testCovariate, filename = filename, colors = colors, ...) {
  bsseq::getMeth(BSseq = bs.filtered.bsseq, regions = sigRegions, type = "smooth", what = "perRegion") %>% na.omit() %>% as.matrix() %>% 
    pheatmap::pheatmap(., scale = "row", annotation_col = bs.filtered.bsseq %>% pData() %>% as.data.frame() %>% 
                         dplyr::select_if(~nlevels(.) > 1), color = RColorBrewer::brewer.pal(11, name = "RdBu") %>% rev(), show_colnames = TRUE, 
                       border_color = "grey", main = glue::glue("Z-Scores of {length(sigRegions)} Differentially Methylated Regions"), 
                       fontsize = 16, filename = filename, width = 11, height = 8.5, cellwidth = 30, annotation_colors = colors %>% 
                         setNames(bs.filtered.bsseq %>% pData() %>% as.data.frame() %>% purrr::pluck(testCovariate) %>% unique() %>% 
                                    sort() %>% rev()) %>% list(testCovariate = .) %>% setNames(testCovariate), ...) %>% return()}

# Set global variables ------------------
genome <- as.character("mm9") # Options: hg38, hg19, mm10, mm9, rheMac10, rheMac8, rn6, danRer11, galGal6, bosTau9, panTro6, dm6, susScr11, canFam3, TAIR10, or TAIR9)
coverage <- as.integer(1) # CpG coverage cutoff minimum value is 1
perGroup <- as.double(1) # Percent of samples in all combinations of covariates meeting CpG coverage cutoff; Options: 0-1
minCpGs <- as.integer(5) # Minimum number of CpGs for a DMR
maxPerms <- as.integer(10) # Maximum number of permutations for the DMR analysis; no more than the # of samples
cutoff <- as.double(0.1) # Cutoff value [from 0 to 1] for the single CpG coefficient utilized to discover testable background regions
testCovariate <- as.character("Stage") # Test covariate 
adjustCovarite <- NULL
cores <- 20
EnsDb <- FALSE

# Setup annotation databases and load files ----------------------------------------------
load("RData/bismark_mNPTM.RData") # methylation values
load("RData/bsseq_mNPTM.RData") # smoothened methylation values
load("RData/settings_mNPTM.RData")
DMRichR::annotationDatabases(genome = genome, EnsDb = EnsDb)

# mM vs mN/mP DMRs (q < .01) ---------------------
# Modify bs.filtered
bs.filteredMPN <- bs.filtered[, which(bs.filtered$Stage != "Tumor")]
bs.filteredMPN$Diagnosis <- rep(c("Metastasis", "Control"), each = 4)

MPN_regions <- dmrseq::dmrseq(bs = bs.filteredMPN, cutoff = cutoff, minNumRegion = minCpGs, maxPerms = maxPerms, testCovariate = "Diagnosis")
MPN_regions <- MPN_regions %>% 
  plyranges::mutate(direction = dplyr::case_when(stat > 0 ~ "Hypermethylated",
                                                 stat < 0 ~ "Hypomethylated"),
                    difference = round(beta/pi *100))

if(sum(MPN_regions$qval < 0.01) >= 100){
  MPN_sigRegions <- MPN_regions %>%
    plyranges::filter(qval < 0.01)
}else if(sum(MPN_regions$qval < 0.01) == 0){
  stop(glue::glue("No significant DMRs detected in {length(MPN_regions)} background regions"))
}
dir.create("mNPTM_DMRs")
gr2bed(MPN_sigRegions, "mNPTM_DMRs/MPN01_DMRs.bed")
gr2bed(MPN_regions, "mNPTM_DMRs/MPN_backgroundRegions.bed")
save(MPN_regions, MPN_sigRegions, file = "RData/MPN_DMRs.RData")
#load("RData/MPN_DMRs.RData")

MPN_sigRegions %>%
  DMRichR::annotateRegions(TxDb = TxDb, annoDb = annoDb) %T>%
  DMRichR::DMReport(regions = MPN_regions, bs.filtered = bs.filteredMPN, coverage = coverage, name = "mNPTM_DMRs/MPN01_DMReport") %>% 
  openxlsx::write.xlsx(file = "mNPTM_DMRs/MPN01_DMRs_annotated.xlsx") # annotate DMRs

print(glue::glue("Extracting individual smoothed methylation values of DMRs..."))
bs.filtered.bsseq %>%
  DMRichR::smooth2txt(regions = MPN_sigRegions,
                      txt = "mNPTM_DMRs/MPN_DMR_individual_smoothed_methylation.txt")

print(glue::glue("Extracting individual smoothed methylation values of background regions..."))
bs.filtered.bsseq %>%
  DMRichR::smooth2txt(regions = MPN_regions,
                      txt = "mNPTM_DMRs/MPN_background_region_individual_smoothed_methylation.txt")

callback <- function(hc, mat){
  sv <- svd(t(mat))$v[,c(1)]
  dend <- reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}

MPN_sigRegions %>% heatmap(bs.filtered.bsseq = bs.filtered.bsseq, testCovariate = testCovariate, 
                           filename = "mNPTM_DMRs/MPN01_NPTM_heatmap.pdf", 
                           colors = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C"), clustering_callback = callback)
MPN_sigRegions %>% heatmap(bs.filtered.bsseq = bs.filtered.bsseq[, which(bs.filtered.bsseq$Stage != "Tumor")], testCovariate = testCovariate, 
                           filename = "mNPTM_DMRs/MPN01_heatmap.pdf", 
                           colors = c("#4DAF4A", "#377EB8", "#E41A1C"))
rm(bs.filteredMPN, callback)

# mT vs mM DMRs (q < .01) ---------------------
# Modify bs.filtered
bs.filteredTM <- bs.filtered[, which(bs.filtered$Stage == "Tumor" | bs.filtered$Stage == "Metastasis")]
bs.filteredTM$Diagnosis <- rep(c("Metastasis", "Control"), each = 4)

TM_regions <- dmrseq::dmrseq(bs = bs.filteredTM, cutoff = cutoff, minNumRegion = minCpGs, maxPerms = maxPerms, testCovariate = "Diagnosis")
TM_regions <- TM_regions %>% 
  plyranges::mutate(direction = dplyr::case_when(stat > 0 ~ "Hypermethylated",
                                                 stat < 0 ~ "Hypomethylated"),
                    difference = round(beta/pi *100))

if(sum(TM_regions$qval < 0.01) >= 100){
  TM_sigRegions <- TM_regions %>%
    plyranges::filter(qval < 0.01)
}else if(sum(TM_regions$qval < 0.01) == 0){
  stop(glue::glue("No significant DMRs detected in {length(TM_regions)} background regions"))
}
gr2bed(TM_sigRegions, "mNPTM_DMRs/TM01_DMRs.bed")
gr2bed(TM_regions, "mNPTM_DMRs/TM_backgroundRegions.bed")
save(TM_regions, TM_sigRegions, file = "RData/TM_DMRs.RData")
#load("RData/TM_DMRs.RData")

TM_sigRegions %>%
  DMRichR::annotateRegions(TxDb = TxDb, annoDb = annoDb) %T>%
  DMRichR::DMReport(regions = TM_regions, bs.filtered = bs.filteredTM, coverage = coverage, name = "mNPTM_DMRs/TM01_DMReport") %>% 
  openxlsx::write.xlsx(file = "mNPTM_DMRs/TM01_DMRs_annotated.xlsx") # annotate DMRs

print(glue::glue("Extracting individual smoothed methylation values of DMRs..."))
bs.filtered.bsseq %>%
  DMRichR::smooth2txt(regions = TM_sigRegions,
                      txt = "mNPTM_DMRs/TM_DMR_individual_smoothed_methylation.txt")

print(glue::glue("Extracting individual smoothed methylation values of background regions..."))
bs.filtered.bsseq %>%
  DMRichR::smooth2txt(regions = TM_regions,
                      txt = "mNPTM_DMRs/TM_background_region_individual_smoothed_methylation.txt")

callback <- function(hc, mat){
  sv <- svd(t(mat))$v[,c(2)]
  dend <- reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}

TM_sigRegions %>% heatmap(bs.filtered.bsseq = bs.filtered.bsseq, testCovariate = as.character("Stage"), 
                          filename = "mNPTM_DMRs/TM01_NPTM_heatmap.pdf", 
                       colors = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C"), clustering_callback = callback)
TM_sigRegions %>% heatmap(bs.filtered.bsseq = bs.filtered.bsseq[, which(bs.filtered.bsseq$Stage == "Tumor" | bs.filtered.bsseq$Stage == "Metastasis")],
                          testCovariate = as.character("Stage"), 
                          filename = "mNPTM_DMRs/TM01_heatmap.pdf",
                       colors = c("#984EA3", "#E41A1C")) # TM .01 heatmap
rm(bs.filteredTM, callback)

# mN/mP vs mT DMRs (q < .01) ---------------------
# Modify bs.filtered
bs.filteredTPN <- bs.filtered[, which(bs.filtered$Stage != "Metastasis")]
bs.filteredTPN$Diagnosis <- rep(c("Tumor", "Control"), each = 4)

TPN_regions <- dmrseq::dmrseq(bs = bs.filteredTPN,  cutoff = cutoff, minNumRegion = minCpGs, maxPerms = maxPerms, testCovariate = "Diagnosis")
TPN_regions <- TPN_regions %>% 
  plyranges::mutate(direction = dplyr::case_when(stat > 0 ~ "Hypermethylated",
                                                 stat < 0 ~ "Hypomethylated"),
                    difference = round(beta/pi *100))
if(sum(TPN_regions$qval < 0.01) >= 100){
  TPN_sigRegions <- TPN_regions %>%
    plyranges::filter(qval < 0.01)
}else if(sum(TPN_regions$qval < 0.01) == 0){
  stop(glue::glue("No significant DMRs detected in {length(TPN_regions)} background regions"))
}

gr2bed(TPN_regions, "mNPTM_DMRs/TPN_backgroundRegions.bed")
save(TPN_regions, file = "RData/TPN_BG.RData")
rm(TPN_regions, bs.filteredTPN)

# Identify overlaps (cutoff .01) ---------------
library(eulerr)

TM_hyper_sigRegions <- TM_sigRegions %>% plyranges::filter(stat > 0)
MPN_hyper_sigRegions <- MPN_sigRegions %>% plyranges::filter(stat > 0)
TM_hypo_sigRegions <- TM_sigRegions %>% plyranges::filter(stat < 0)
MPN_hypo_sigRegions <- MPN_sigRegions %>% plyranges::filter(stat < 0)
hyper_overlap <- intersect(MPN_hyper_sigRegions, TM_hyper_sigRegions, ignore.strand = T)
hypo_overlap <- intersect(MPN_hypo_sigRegions, TM_hypo_sigRegions, ignore.strand = T)
overlap <- length(hyper_overlap) + length(hypo_overlap)
TMonly <- length(TM_sigRegions) - overlap
MPNonly <- length(MPN_sigRegions) - overlap
fit = euler(c("mT vs mM" = TMonly, "mN/mP vs mM" = MPNonly, 
              "mT vs mM&mN/mP vs mM" = overlap))
pdf(file = "mNPTM_DMRs/euler01.pdf")
plot(fit, quantities = TRUE, legend = list(lables = c("mT vs mM", "mN/mP vs mM")),
     fills = list(fill = c("#66D2D6", "#E56997"), alpha = 0.8))
dev.off()
rm(hyper_overlap, hypo_overlap, overlap, TMonly, MPNonly, fit)

# CpG and genic enrichment testing for MPN ----------------------------------------
dir.create("mNPTM_DMRichments")
DMRich <- function(x){
    dmrList[x] %>% 
      DMRichR::DMRichCpG(regions = MPN_regions, genome = genome) %T>%
      openxlsx::write.xlsx(file = glue::glue("mNPTM_DMRichments/MPN_{names(dmrList)[x]}_CpG_enrichments.xlsx")) %>% 
      DMRichR::DMRichPlot(type = "CpG") %>% 
      ggplot2::ggsave(glue::glue("mNPTM_DMRichments/MPN_{names(dmrList)[x]}_CpG_enrichments.pdf"), plot = ., width = 4, height = 3)
    dmrList[x] %>% 
      DMRichR::DMRichGenic(regions = MPN_regions, TxDb = TxDb, annoDb = annoDb) %T>%
      openxlsx::write.xlsx(file = glue::glue("mNPTM_DMRichments/MPN_{names(dmrList)[x]}_genic_enrichments.xlsx")) %>% 
      DMRichR::DMRichPlot(type = "genic") %>% 
      ggplot2::ggsave(glue::glue("mNPTM_DMRichments/MPN_{names(dmrList)[x]}_genic_enrichments.pdf"), plot = ., width = 4, height = 4)
}
dmrList <- MPN_sigRegions %>% DMRichR::dmrList()
parallel::mclapply(seq_along(dmrList), DMRich, mc.cores = 1, mc.silent = TRUE)

DMparseR <- function (direction = c("All DMRs", "Hypermethylated DMRs", 
                        "Hypomethylated DMRs"), type = c("CpG", "genic")) {
  stopifnot(direction %in% c("All DMRs", "Hypermethylated DMRs", 
                             "Hypomethylated DMRs"))
  stopifnot(type %in% c("CpG", "genic"))
  print(glue::glue("Parsing {type} enrichment results for {tidyDirection}", 
                   tidyDirection = glue::glue_collapse({
                     direction
                   }, sep = ", ", last = " and ")))
  purrr::map(direction, function(direction) {
    glue::glue("mNPTM_DMRichments/MPN_{direction}_{type}_enrichments.xlsx")
  }) %>% as.vector() %>% lapply(function(file) {
    readxl::read_xlsx(file)
  }) %>% `names<-`(direction) %>% data.table::rbindlist(idcol = "Dataset") %>% 
    dplyr::as_tibble() %>% tidyr::separate(Dataset, c("Direction", "DMR")) %>% 
    dplyr::select(Direction, Annotation, OR, fdr) %>% 
    dplyr::mutate(Annotation = forcats::as_factor(Annotation)) %>% 
    return()
}

purrr::walk(dplyr::case_when(genome %in% c("hg38", "hg19", "mm10", "mm9", "rn6") ~ c("CpG", "genic"), TRUE ~ "genic") %>% unique(),
            function(type){
              DMparseR(direction =  c("All DMRs","Hypermethylated DMRs","Hypomethylated DMRs"),type = type) %>%
                DMRichR::DMRichPlot(type = type,multi = TRUE) %>% 
                ggplot2::ggsave(glue::glue("mNPTM_DMRichments/MPN_{type}_multi_plot.pdf"), plot = ., device = NULL,
                                height = dplyr::case_when(type == "genic" ~ 5, type == "CpG" ~ 3.5), width = 7)
            })

# Genic enrichment counts
MPN_hyper_annot <- MPN_hyper_sigRegions %>% 
  DMRichR::annotateRegions(TxDb = TxDb, annoDb = annoDb) %>%
  dplyr::select(annotation) %>% 
  table() %>% as.data.frame() # hypermethylated frequency table
MPN_hyper_annot$Percent <- MPN_hyper_annot[[2]]/sum(MPN_hyper_annot[[2]]) # add percent

MPN_hypo_annot <- MPN_hypo_sigRegions %>% 
  DMRichR::annotateRegions(TxDb = TxDb, annoDb = annoDb) %>%
  dplyr::select(annotation) %>% table() %>% as.data.frame() # hypomethylated frquency table# Percent graphs
MPN_hypo_annot$Percent <- MPN_hypo_annot[[2]]/sum(MPN_hypo_annot[[2]]) # add percent

MPN_counts <- data.frame(Annotation = rep(levels(MPN_hyper_annot[[1]]), 2),
                        Count = c(MPN_hyper_annot[[2]], MPN_hypo_annot[[2]]),
                        Percent = round(c(MPN_hyper_annot[[3]], MPN_hypo_annot[[3]]),2),
                        Direction = rep(c("Hypermethylated", "Hypomethylated"), each=7)) # initialize dataframe
MPN_counts$Annotation <- factor(MPN_counts$Annotation, 
                                levels=c("Promoter", "5' UTR", "Exon", "Intron", "3' UTR", "Downstream", "Distal Intergenic")) # reorder
write.table(MPN_counts, file = "mNPTM_DMRichments/MPN_genic_counts.txt", quote = FALSE, sep = '\t ', row.names = F) # saves the results as a text file in the working directory

pdf(file = "mNPTM_DMRichments/MPN_genic_counts.pdf")
ggplot(MPN_counts, aes(fill=Annotation, y=Percent, x=Direction)) + 
  geom_bar(position="dodge", stat = "identity", color = "black")+
  scale_y_continuous(labels = scales::percent) +
  labs(x ="", y = "") + ggtitle("mN/mP vs mM") +
  scale_x_discrete(labels=c("Hypermethylated", "Hypomethylated")) +
  theme_minimal() + 
  scale_fill_manual(values=wesanderson::wes_palette("Zissou1", n = 7, type = "continuous") %>% rev())
dev.off()
rm(MPN_hyper_annot, MPN_hypo_annot, MPN_counts)

# CpG enrichment counts
MPN_hyper_CpG <- MPN_hyper_sigRegions %>% 
  DMRichR::annotateRegions(TxDb = TxDb, annoDb = annoDb) %>%
  dplyr::select(CpG.Island, CpG.Shore, CpG.Shelf, Open.Sea) %>% as.data.frame() # hypermethylated CpGs
yes_hyper <- sapply(MPN_hyper_CpG,FUN = function(x){length(x[x=="Yes"])})
count_hyper <- data.frame(Count = yes_hyper, Percent = yes_hyper/length(MPN_hyper_sigRegions))

MPN_hypo_CpG <- MPN_hypo_sigRegions %>% 
  DMRichR::annotateRegions(TxDb = TxDb, annoDb = annoDb) %>%
  dplyr::select(CpG.Island, CpG.Shore, CpG.Shelf, Open.Sea) %>% as.data.frame() # hypomethylated CpGs
yes_hypo <- sapply(MPN_hypo_CpG,FUN = function(x){length(x[x=="Yes"])})
count_hypo <- data.frame(Count = yes_hypo, Percent = yes_hypo/length(MPN_hypo_sigRegions))

MPN_counts <- data.frame(CpG = rep(rownames(count_hyper),2),
                         Count = c(count_hyper[[1]], count_hypo[[1]]),
                         Percent = round(c(count_hyper[[2]], count_hypo[[2]]),2),
                         Direction = rep(c("Hypermethylated", "Hypomethylated"), each=4)) # initialize dataframe
write.table(MPN_counts, file = "mNPTM_DMRichments/MPN_CpG_counts.txt", quote = FALSE, sep = '\t ', row.names = F) # saves the results as a text file in the working directory

pdf(file = "mNPTM_DMRichments/MPN_CpG_counts.pdf")
ggplot(MPN_counts, aes(fill=CpG, y=Percent, x=Direction)) + 
  geom_bar(position="dodge", stat = "identity", color = "black")+
  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
  labs(x ="", y = "") + ggtitle("mN/mP vs mM") +
  scale_x_discrete(labels=c("Hypermethylated", "Hypomethylated")) +
  theme_minimal() + 
  scale_fill_manual(values = c("forestgreen", "goldenrod2", "dodgerblue", "blue3"))
dev.off()
rm(TM_hyper_CpG, TM_hypo_CpG, yes_hyper, yes_hypo, count_hyper, count_hypo, TM_counts)

# CpG and genic enrichment testing for TM ----------------------------------------
DMRich <- function(x){
  dmrList[x] %>% 
    DMRichR::DMRichCpG(regions = TM_regions, genome = genome) %T>%
    openxlsx::write.xlsx(file = glue::glue("mNPTM_DMRichments/TM_{names(dmrList)[x]}_CpG_enrichments.xlsx")) %>% 
    DMRichR::DMRichPlot(type = "CpG") %>% 
    ggplot2::ggsave(glue::glue("mNPTM_DMRichments/TM_{names(dmrList)[x]}_CpG_enrichments.pdf"), plot = ., width = 4, height = 3)
  dmrList[x] %>% 
    DMRichR::DMRichGenic(regions = TM_regions, TxDb = TxDb, annoDb = annoDb) %T>%
    openxlsx::write.xlsx(file = glue::glue("mNPTM_DMRichments/TM_{names(dmrList)[x]}_genic_enrichments.xlsx")) %>% 
    DMRichR::DMRichPlot(type = "genic") %>% 
    ggplot2::ggsave(glue::glue("mNPTM_DMRichments/TM_{names(dmrList)[x]}_genic_enrichments.pdf"), plot = ., width = 4, height = 4)
}
dmrList <- TM_sigRegions %>% DMRichR::dmrList()
parallel::mclapply(seq_along(dmrList), DMRich, mc.cores = 1, mc.silent = TRUE)

DMparseR <- function (direction = c("All DMRs", "Hypermethylated DMRs", 
                                    "Hypomethylated DMRs"), type = c("CpG", "genic")) {
  stopifnot(direction %in% c("All DMRs", "Hypermethylated DMRs", 
                             "Hypomethylated DMRs"))
  stopifnot(type %in% c("CpG", "genic"))
  print(glue::glue("Parsing {type} enrichment results for {tidyDirection}", 
                   tidyDirection = glue::glue_collapse({
                     direction
                   }, sep = ", ", last = " and ")))
  purrr::map(direction, function(direction) {
    glue::glue("mNPTM_DMRichments/TM_{direction}_{type}_enrichments.xlsx")
  }) %>% as.vector() %>% lapply(function(file) {
    readxl::read_xlsx(file)
  }) %>% `names<-`(direction) %>% data.table::rbindlist(idcol = "Dataset") %>% 
    dplyr::as_tibble() %>% tidyr::separate(Dataset, c("Direction", "DMR")) %>% 
    dplyr::select(Direction, Annotation, OR, fdr) %>% 
    dplyr::mutate(Annotation = forcats::as_factor(Annotation)) %>% 
    return()
}

purrr::walk(dplyr::case_when(genome %in% c("hg38", "hg19", "mm10", "mm9", "rn6") ~ c("CpG", "genic"), TRUE ~ "genic") %>% unique(),
            function(type){
              DMparseR(direction =  c("All DMRs","Hypermethylated DMRs","Hypomethylated DMRs"),type = type) %>%
                DMRichR::DMRichPlot(type = type,multi = TRUE) %>% 
                ggplot2::ggsave(glue::glue("mNPTM_DMRichments/TM_{type}_multi_plot.pdf"), plot = ., device = NULL,
                                height = dplyr::case_when(type == "genic" ~ 5, type == "CpG" ~ 3.5), width = 7)
            })

# Genic enrichment counts
TM_hyper_annot <- TM_hyper_sigRegions %>% 
  DMRichR::annotateRegions(TxDb = TxDb, annoDb = annoDb) %>%
  dplyr::select(annotation) %>% 
  table() %>% as.data.frame() # hypermethylated frequency table
TM_hyper_annot$Percent <- TM_hyper_annot[[2]]/sum(TM_hyper_annot[[2]]) # add percent

TM_hypo_annot <- TM_hypo_sigRegions %>% 
  DMRichR::annotateRegions(TxDb = TxDb, annoDb = annoDb) %>%
  dplyr::select(annotation) %>% table() %>% as.data.frame() # hypomethylated frquency table# Percent graphs
TM_hypo_annot$Percent <- TM_hypo_annot[[2]]/sum(TM_hypo_annot[[2]]) # add percent

TM_counts <- data.frame(Annotation = rep(levels(TM_hyper_annot[[1]]), 2),
                         Count = c(TM_hyper_annot[[2]], TM_hypo_annot[[2]]),
                         Percent = round(c(TM_hyper_annot[[3]], TM_hypo_annot[[3]]),2),
                         Direction = rep(c("Hypermethylated", "Hypomethylated"), each=7)) # initialize dataframe
TM_counts$Annotation <- factor(TM_counts$Annotation, 
                                levels=c("Promoter", "5' UTR", "Exon", "Intron", "3' UTR", "Downstream", "Distal Intergenic")) # reorder
write.table(TM_counts, file = "mNPTM_DMRichments/TM_genic_counts.txt", quote = FALSE, sep = '\t ', row.names = F) # saves the results as a text file in the working directory

pdf(file = "mNPTM_DMRichments/TM_genic_counts.pdf")
ggplot(TM_counts, aes(fill=Annotation, y=Percent, x=Direction)) + 
  geom_bar(position="dodge", stat = "identity", color = "black")+
  scale_y_continuous(labels = scales::percent, limits = c(0,.5)) +
  labs(x ="", y = "") + ggtitle("mT vs mM") +
  scale_x_discrete(labels=c("Hypermethylated", "Hypomethylated")) +
  theme_minimal() + 
  scale_fill_manual(values=wesanderson::wes_palette("Zissou1", n = 7, type = "continuous") %>% rev())
dev.off()
rm(TM_hyper_annot, TM_hypo_annot, TM_counts)

# CpG enrichment counts
TM_hyper_CpG <- TM_hyper_sigRegions %>% 
  DMRichR::annotateRegions(TxDb = TxDb, annoDb = annoDb) %>%
  dplyr::select(CpG.Island, CpG.Shore, CpG.Shelf, Open.Sea) %>% as.data.frame() # hypermethylated CpGs
yes_hyper <- sapply(TM_hyper_CpG,FUN = function(x){length(x[x=="Yes"])})
count_hyper <- data.frame(Count = yes_hyper, Percent = yes_hyper/length(TM_hyper_sigRegions))

TM_hypo_CpG <- TM_hypo_sigRegions %>% 
  DMRichR::annotateRegions(TxDb = TxDb, annoDb = annoDb) %>%
  dplyr::select(CpG.Island, CpG.Shore, CpG.Shelf, Open.Sea) %>% as.data.frame() # hypomethylated CpGs
yes_hypo <- sapply(TM_hypo_CpG,FUN = function(x){length(x[x=="Yes"])})
count_hypo <- data.frame(Count = yes_hypo, Percent = yes_hypo/length(TM_hypo_sigRegions))

TM_counts <- data.frame(CpG = rep(rownames(count_hyper),2),
                         Count = c(count_hyper[[1]], count_hypo[[1]]),
                         Percent = round(c(count_hyper[[2]], count_hypo[[2]]),2),
                         Direction = rep(c("Hypermethylated", "Hypomethylated"), each=4)) # initialize dataframe
write.table(TM_counts, file = "mNPTM_DMRichments/TM_CpG_counts.txt", quote = FALSE, sep = '\t ', row.names = F) # saves the results as a text file in the working directory

pdf(file = "mNPTM_DMRichments/TM_CpG_counts.pdf")
ggplot(TM_counts, aes(fill=CpG, y=Percent, x=Direction)) + 
  geom_bar(position="dodge", stat = "identity", color = "black")+
  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
  labs(x ="", y = "") + ggtitle("mT vs mM") +
  scale_x_discrete(labels=c("Hypermethylated", "Hypomethylated")) +
  theme_minimal() + 
  scale_fill_manual(values = c("forestgreen", "goldenrod2", "dodgerblue", "blue3"))
dev.off()
rm(MPN_hyper_CpG, MPN_hypo_CpG, yes_hyper, yes_hypo, count_hyper, count_hypo, MPN_counts)

# Manhattan plots -------------------------------------------------
MPN_regions %>% 
  DMRichR::annotateRegions(TxDb = TxDb, annoDb = annoDb) %>%
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>%
  sort() %>% dplyr::as_tibble() %>% dplyr::select(geneSymbol, seqnames, start, q.value)  %>% 
  dplyr::mutate(seqnames = substring(.$seqnames, 4)) %>%
  CMplot::CMplot(col = c("grey30", "grey60"), plot.type = "m", 
                 LOG10 = TRUE, ylim = NULL, threshold = 0.01, threshold.lty = 1, threshold.lwd = 1, 
                 threshold.col = "black", cex = 0.5, cex.axis = 0.7, amplify = FALSE, chr.den.col = c("darkgreen", "yellow", "red"), 
                 bin.size = 1e+06, file = "pdf", memo = "")
file.rename("Rectangular-Manhattan.q.value.pdf", "mNPTM_DMRs/MPN_manhattan.pdf")

TM_regions %>% 
  DMRichR::annotateRegions(TxDb = TxDb, annoDb = annoDb) %>%
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>%
  sort() %>% dplyr::as_tibble() %>% dplyr::select(geneSymbol, seqnames, start, q.value)  %>% 
  dplyr::mutate(seqnames = substring(.$seqnames, 4)) %>%
  CMplot::CMplot(col = c("grey30", "grey60"), plot.type = "m", 
                 LOG10 = TRUE, ylim = NULL, threshold = 0.01, threshold.lty = 1, threshold.lwd = 1, 
                 threshold.col = "black", cex = 0.5, cex.axis = 0.7, amplify = FALSE, chr.den.col = c("darkgreen", "yellow", "red"), 
                 bin.size = 1e+06, file = "pdf", memo = "")
file.rename("Rectangular-Manhattan.q.value.pdf", "mNPTM_DMRs/TM_manhattan.pdf")

# Prepare HOMER -------------------------------------------------------------------
prepHOMER <- function (sigRegions = sigRegions, regions = regions, dir.name = dir.name) 
{
  dir.create(dir.name)
  sigRegions %>% DMRichR::gr2bed(paste(dir.name,"DMRs.bed", sep="/"))
  sigRegions %>% plyranges::filter(stat > 0) %>% DMRichR::gr2bed(paste(dir.name,"DMRs_hyper.bed", sep="/"))
  sigRegions %>% plyranges::filter(stat < 0) %>% DMRichR::gr2bed(paste(dir.name,"DMRs_hypo.bed", sep="/"))
  regions %>% DMRichR::gr2bed(paste(dir.name,"background.bed", sep="/"))
}

MPN_sigRegions %>% prepHOMER(regions = MPN_regions, dir.name = "mNPTM_DMRs/MPN_HOMER")
TM_sigRegions %>% prepHOMER(regions = TM_regions, dir.name = "mNPTM_DMRs/TM_HOMER")

# Machine learning --------------------------------------------------------
setwd("mNPTM_DMRs")
bs.filtered.bsseq$Diagnosis = bs.filtered.bsseq$Stage
tryCatch({methylLearnOutput <- DMRichR::methylLearn(bs.filtered.bsseq = bs.filtered.bsseq,
                                            sigRegions = MPN_sigRegions,
                                            testCovariate = "Stage",
                                            TxDb = TxDb,
                                            annoDb = annoDb,
                                            topPercent = 1,
                                            output = "all",
                                            saveHtmlReport = TRUE)
    dir.create("./Machine_learning")
  if(length(methylLearnOutput) == 1) {
    openxlsx::write.xlsx(list(Annotations_Common_DMRs = methylLearnOutput), 
                         file = "./Machine_learning/Machine_learning_output_one.xlsx") 
  } else {
    openxlsx::write.xlsx(list(Annotations_Common_DMRs = methylLearnOutput$`Annotated common DMRs`,
                              RF_Ranking_All_DMRs = methylLearnOutput$`RF ranking`,
                              SVM_Ranking_All_DMRs = methylLearnOutput$`SVM ranking`),
                         file = "./Machine_learning/Machine_learning_output_all.xlsx") 
  }
  save(methylLearnOutput, file = "RData/machineLearning.RData")
  #load("RData/machineLearing.RData")
},
error = function(error_condition) {
  print(glue::glue("Warning: methylLearn did not finish. \\
                      You may have not had enough top DMRs across algorithms."))
})
rm(methylLearnOutput)