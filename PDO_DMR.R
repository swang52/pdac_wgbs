# Description
# DMR calling for PDOs (hN vs hT and early vs late) + enrichment analysis and machine learning

# Get Packages if haven't already --------------------------------------------------------------
packages <- c("ggplot2", "DMRichR", "eulerr", "RColorBrewer")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Load packages
invisible(lapply(packages, library, character.only = TRUE))

# Heatmap function --------------------
heatmap <- function(bs.filtered.bsseq = bs.filtered.bsseq, sigRegions = sigRegions, testCovariate = testCovariate, filename = filename, colors = colors, ...) {
  bsseq::getMeth(BSseq = bs.filtered.bsseq, regions = sigRegions, type = "smooth", what = "perRegion") %>% na.omit() %>% as.matrix() %>% 
    pheatmap::pheatmap(., scale = "row", annotation_col = bs.filtered.bsseq %>% pData() %>% as.data.frame() %>%
                         dplyr::select_if(~nlevels(.) > 1), color = RColorBrewer::brewer.pal(11, name = "RdBu") %>% rev(), show_colnames = TRUE, 
                       border_color = "grey", main = glue::glue("Z-Scores of {length(sigRegions)} Differentially Methylated Regions"), 
                       fontsize = 10, filename = filename, width = 11, height = 8.5, cellwidth = 12, annotation_colors = colors %>% 
                         setNames(bs.filtered.bsseq %>% pData() %>% as.data.frame() %>% purrr::pluck(testCovariate) %>% unique() %>% 
                                    sort() %>% rev()) %>% list(testCovariate = .) %>% setNames(testCovariate), ...) %>% return()}

# Set global variables ------------------
genome <- as.character("hg38") # Options: hg38, hg19, mm10, mm9, rheMac10, rheMac8, rn6, danRer11, galGal6, bosTau9, panTro6, dm6, susScr11, canFam3, TAIR10, or TAIR9)
coverage <- as.integer(1) # CpG coverage cutoff minimum value is 1
perGroup <- as.double(.75) # Percent of samples in all combinations of covariates meeting CpG coverage cutoff; Options: 0-1
minCpGs <- as.integer(5) # Minimum number of CpGs for a DMR
maxPerms <- as.integer(10) # Maximum number of permutations for the DMR analysis; no more than the # of samples
cutoff <- as.double(0.1) # Cutoff value [from 0 to 1] for the single CpG coefficient utilized to discover testable background regions
testCovariate <- as.character("CombinedStage") # Test covariate 
adjustCovarite <- NULL
cores <- 20
EnsDb <- FALSE

# Setup annotation databases and load files ----------------------------------------------
load("RData/bismark_NPDOs.RData") # methylation values
load("RData/bsseq_NPDOs.RData") # smoothened methylation values

DMRichR::annotationDatabases(genome = genome, EnsDb = EnsDb)

# Normal vs Tumor DMRs (p < .01) ---------------------
# Modify bs.filtered
bs.filtered$CombinedStage = bs.filtered$Stage %>% as.data.frame() %>% 
  with(., ifelse(. == "Normal", "Normal", "Tumor")) %>% as.factor()

NT_regions <- dmrseq::dmrseq(bs = bs.filtered, cutoff = cutoff, minNumRegion = minCpGs, maxPerms = maxPerms, testCovariate = testCovariate)
NT_regions <- NT_regions %>% 
  plyranges::mutate(direction = dplyr::case_when(stat > 0 ~ "Hypermethylated",
                                                 stat < 0 ~ "Hypomethylated"),
                    difference = round(beta/pi *100))

if(sum(NT_regions$qval < 0.05) < 100 & sum(NT_regions$pval < 0.01) != 0){
  NT_sigRegions <- NT_regions %>%
    plyranges::filter(pval < 0.01)
}else if(sum(NT_regions$qval < 0.05) >= 100){
  NT_sigRegions <- NT_regions %>%
    plyranges::filter(qval < 0.05)
}else if(sum(NT_regions$pval < 0.05) == 0){
  stop(glue::glue("No significant DMRs detected in {length(regions)} background regions"))
}

dir.create("PDO_DMRs")
gr2bed(NT_sigRegions, "PDO_DMRs/NT_DMRs.bed")
gr2bed(NT_regions, "PDO_DMRs/NT_backgroundRegions.bed")
save(NT_regions, NT_sigRegions, file = "RData/NT_DMRs.RData")
#load("RData/NT_DMRs.RData")

NT_sigRegions %>%
  DMRichR::annotateRegions(TxDb = TxDb, annoDb = annoDb) %>%
  DMRichR::DMReport(regions = NT_regions, bs.filtered = bs.filtered, coverage = coverage, name = "PDO_DMRs/NT_DMReport") %>% 
  openxlsx::write.xlsx(file = "PDO_DMRs/NT_DMRs_annotated.xlsx") # annotate DMRs

print(glue::glue("Extracting individual smoothed methylation values of DMRs..."))
bs.filtered.bsseq %>%
  DMRichR::smooth2txt(regions = NT_sigRegions,
                      txt = "PDO_DMRs/NT_DMR_individual_smoothed_methylation.txt")

print(glue::glue("Extracting individual smoothed methylation values of background regions..."))
bs.filtered.bsseq %>%
  DMRichR::smooth2txt(regions = NT_regions,
                      txt = "PDO_DMRs/NT_background_region_individual_smoothed_methylation.txt")

colors <- list(c("#4DAF4A", "#377EB8", "#E41A1C", "#FF7F00", "#984EA3") %>% 
                 setNames(bs.filtered.bsseq %>% pData() %>% as.data.frame() %>% purrr::pluck("Stage") %>% 
                            unique() %>% sort() %>% rev()),
               c("#7FC97F", "#BEAED4", "#666666", "#FDC086", "#FFFF99") %>% setNames(bs.filtered.bsseq %>% pData() %>% as.data.frame() %>% purrr::pluck("Subtype") %>% 
                            unique() %>% sort() %>% rev()),
                c("#8DD3C7", "#FCCDE5") %>% setNames(bs.filtered.bsseq %>% pData() %>% as.data.frame() %>% purrr::pluck("CombinedStage") %>% 
                            unique() %>% sort() %>% rev()))
names(colors) = c("Stage", "Subtype", "CombinedStage")
NT_sigRegions %>% heatmap(bs.filtered.bsseq = bs.filtered.bsseq, testCovariate = testCovariate,
                            filename = "PDO_DMRs/NT_heatmap.pdf", 
                            colors = colors)

# early vs late DMRs (p < .01) ---------------------
# Modify bs.filtered
bs.filteredPDO <- bs.filtered[, which(bs.filtered$Stage != "Normal")]
bs.filteredPDO$Stage = droplevels(bs.filteredPDO$Stage)
bs.filteredPDO$CombinedStage = bs.filteredPDO$Stage %>% as.data.frame() %>% 
  with(., ifelse(. == "Metastatic" | . == "Locally advanced", "Late", "Early")) %>% as.factor()

EL_regions <- dmrseq::dmrseq(bs = bs.filteredPDO, cutoff = cutoff, minNumRegion = minCpGs, maxPerms = maxPerms, testCovariate = testCovariate)
EL_regions <- EL_regions %>% 
  plyranges::mutate(direction = dplyr::case_when(stat > 0 ~ "Hypermethylated",
                                                 stat < 0 ~ "Hypomethylated"),
                    difference = round(beta/pi *100))

if(sum(EL_regions$qval < 0.05) < 100 & sum(EL_regions$pval < 0.01) != 0){
  EL_sigRegions <- EL_regions %>%
    plyranges::filter(pval < 0.01)
}else if(sum(EL_regions$qval < 0.05) >= 100){
  EL_sigRegions <- EL_regions %>%
    plyranges::filter(qval < 0.05)
}else if(sum(EL_regions$pval < 0.05) == 0){
  stop(glue::glue("No significant DMRs detected in {length(regions)} background regions"))
}

gr2bed(EL_sigRegions, "PDO_DMRs/EL_DMRs.bed")
gr2bed(EL_regions, "PDO_DMRs/EL_backgroundRegions.bed")
save(EL_regions, EL_sigRegions, file = "RData/EL_DMRs.RData")
#load("RData/EL_DMRs.RData")

EL_sigRegions %>%
  DMRichR::annotateRegions(TxDb = TxDb, annoDb = annoDb) %>%
  DMRichR::DMReport(regions = EL_regions, bs.filtered = bs.filteredPDO, coverage = coverage, name = "PDO_DMRs/EL_DMReport") %>% 
  openxlsx::write.xlsx(file = "PDO_DMRs/EL_DMRs_annotated.xlsx") # annotate DMRs

bs.filtered.bsseqPDO <- bs.filtered.bsseq[, which(bs.filtered.bsseq$Stage != "Normal")]
bs.filtered.bsseqPDO$Stage = droplevels(bs.filtered.bsseqPDO$Stage)
bs.filtered.bsseqPDO$CombinedStage = bs.filtered.bsseqPDO$Stage %>% as.data.frame() %>% 
  with(., ifelse(. == "Metastatic" | . == "Locally advanced", "Late", "Early")) %>% as.factor()

print(glue::glue("Extracting individual smoothed methylation values of DMRs..."))
bs.filtered.bsseqPDO %>%
  DMRichR::smooth2txt(regions = EL_sigRegions,
                      txt = "PDO_DMRs/EL_DMR_individual_smoothed_methylation.txt")

print(glue::glue("Extracting individual smoothed methylation values of background regions..."))
bs.filtered.bsseqPDO %>%
  DMRichR::smooth2txt(regions = EL_regions,
                      txt = "PDO_DMRs/EL_background_region_individual_smoothed_methylation.txt")

colors <- list(c("#4DAF4A", "#E41A1C", "#FF7F00", "#984EA3") %>% 
                 setNames(bs.filtered.bsseqPDO %>% pData() %>% as.data.frame() %>% purrr::pluck("Stage") %>% 
                            unique() %>% sort() %>% rev()),
               c("#7FC97F", "#BEAED4", "#FDC086", "#FFFF99") %>% setNames(bs.filtered.bsseqPDO %>% pData() %>% as.data.frame() %>% purrr::pluck("Subtype") %>% 
                                                                            unique() %>% sort() %>% rev()),
               c("#8DD3C7", "#FCCDE5") %>% setNames(bs.filtered.bsseqPDO %>% pData() %>% as.data.frame() %>% purrr::pluck("CombinedStage") %>% 
                                                      unique() %>% sort() %>% rev()))
names(colors) = c("Stage", "Subtype", "CombinedStage")
EL_sigRegions %>% heatmap(bs.filtered.bsseq = bs.filtered.bsseqPDO, testCovariate = as.character("Stage"), 
                          filename = "PDO_DMRs/EL_heatmap.pdf", 
                          colors = colors)
EL_sigRegions %>% heatmap(bs.filtered.bsseq = bs.filtered.bsseqPDO, testCovariate = testCovariate, 
                           sigRegions = EL_sigRegions, filename = "PDO_DMRs/EL_heatmap.pdf", 
                           colors = colors)

# Identify overlaps (cutoff .01) ---------------
EL_hyper_sigRegions <- EL_sigRegions %>% plyranges::filter(stat > 0)
NT_hyper_sigRegions <- NT_sigRegions %>% plyranges::filter(stat > 0)
EL_hypo_sigRegions <- EL_sigRegions %>% plyranges::filter(stat < 0)
NT_hypo_sigRegions <- NT_sigRegions %>% plyranges::filter(stat < 0)
hyper_overlap <- intersect(NT_hyper_sigRegions, EL_hyper_sigRegions, ignore.strand = T)
hypo_overlap <- intersect(NT_hypo_sigRegions, EL_hypo_sigRegions, ignore.strand = T)
overlap <- length(hyper_overlap) + length(hypo_overlap)
ELonly <- length(EL_sigRegions) - overlap
NTonly <- length(NT_sigRegions) - overlap
fit = euler(c("Early vs Late" = ELonly, "hN vs hT" = NTonly, 
              "Early vs Late&hN vs hT" = overlap))
pdf(file = "PDO_DMRs/euler01.pdf")
plot(fit, quantities = TRUE, legend = list(lables = c("Early vs Late", "hN vs hT")),
     fills = list(fill = c("#66D2D6", "#E56997"), alpha = 0.8))
dev.off()
rm(hyper_overlap, hypo_overlap, overlap, ELonly, NTonly, fit)

# CpG and genic enrichment testing for Early vs Late ----------------------------------------
dir.create("PDO_DMRichments")
DMRich <- function(x){
  dmrList[x] %>% 
    DMRichR::DMRichCpG(regions = EL_regions, genome = genome) %T>%
    openxlsx::write.xlsx(file = glue::glue("PDO_DMRichments/EL_{names(dmrList)[x]}_CpG_enrichments.xlsx")) %>% 
    DMRichR::DMRichPlot(type = "CpG") %>% 
    ggplot2::ggsave(glue::glue("PDO_DMRichments/EL_{names(dmrList)[x]}_CpG_enrichments.pdf"),
                    plot = ., width = 4, height = 3)
  dmrList[x] %>% 
    DMRichR::DMRichGenic(regions = EL_regions, TxDb = TxDb, annoDb = annoDb) %T>%
    openxlsx::write.xlsx(file = glue::glue("PDO_DMRichments/EL_{names(dmrList)[x]}_genic_enrichments.xlsx")) %>% 
    DMRichR::DMRichPlot(type = "genic") %>% 
    ggplot2::ggsave(glue::glue("PDO_DMRichments/EL_{names(dmrList)[x]}_genic_enrichments.pdf"), plot = ., width = 4, height = 4)
}
dmrList <- EL_sigRegions %>% DMRichR::dmrList()
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
    glue::glue("PDO_DMRichments/EL_{direction}_{type}_enrichments.xlsx")
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
                ggplot2::ggsave(glue::glue("PDO_DMRichments/EL_{type}_multi_plot.pdf"), plot = ., device = NULL,
                                height = dplyr::case_when(type == "genic" ~ 5, type == "CpG" ~ 3.5), width = 7)
            })

# Genic enrichment counts
EL_hyper_annot <- EL_hyper_sigRegions %>% 
  DMRichR::annotateRegions(TxDb = TxDb, annoDb = annoDb) %>%
  dplyr::select(annotation) %>% 
  table() %>% as.data.frame() # hypermethylated frequency table
EL_hyper_annot$Percent <- EL_hyper_annot[[2]]/sum(EL_hyper_annot[[2]]) # add percent

EL_hypo_annot <- EL_hypo_sigRegions %>% 
  DMRichR::annotateRegions(TxDb = TxDb, annoDb = annoDb) %>%
  dplyr::select(annotation) %>% table() %>% as.data.frame() # hypomethylated frquency table# Percent graphs
EL_hypo_annot$Percent <- EL_hypo_annot[[2]]/sum(EL_hypo_annot[[2]]) # add percent

EL_counts <- data.frame(Annotation = c(levels(EL_hyper_annot[[1]]), levels(EL_hypo_annot[[1]])),
                         Count = c(EL_hyper_annot[[2]], EL_hypo_annot[[2]]),
                         Percent = round(c(EL_hyper_annot[[3]], EL_hypo_annot[[3]]),2),
                         Direction = c(rep("Hypermethylated",length(EL_hyper_annot[[1]])), 
                                       rep("Hypomethylated", length(EL_hypo_annot[[1]])))) # initialize dataframe
EL_counts$Annotation <- factor(EL_counts$Annotation, 
                                levels=c("Promoter", "5' UTR", "Exon", "Intron", "3' UTR", "Downstream", "Distal Intergenic")) # reorder
write.table(EL_counts, file = "PDO_DMRichments/EL_genic_counts.txt", quote = FALSE, sep = '\t ', row.names = F) # saves the results as a text file in the working directory

pdf(file = "PDO_DMRichments/EL_genic_counts.pdf")
ggplot(EL_counts, aes(fill=Annotation, y=Percent, x=Direction)) + 
  geom_bar(position="dodge", stat = "identity", color = "black")+
  scale_y_continuous(labels = scales::percent) +
  labs(x ="", y = "") + ggtitle("Early vs Late") +
  scale_x_discrete(labels=c("Hypermethylated", "Hypomethylated")) +
  theme_minimal() + 
  scale_fill_manual(values=wesanderson::wes_palette("Zissou1", n = 7, type = "continuous") %>% rev())
dev.off()
rm(EL_hyper_annot, EL_hypo_annot, EL_counts)

# CpG enrichment counts
EL_hyper_CpG <- EL_hyper_sigRegions %>% 
  DMRichR::annotateRegions(TxDb = TxDb, annoDb = annoDb) %>%
  dplyr::select(CpG.Island, CpG.Shore, CpG.Shelf, Open.Sea) %>% as.data.frame() # hypermethylated CpGs
yes_hyper <- sapply(EL_hyper_CpG,FUN = function(x){length(x[x=="Yes"])})
count_hyper <- data.frame(Count = yes_hyper, Percent = yes_hyper/length(EL_hyper_sigRegions))

EL_hypo_CpG <- EL_hypo_sigRegions %>% 
  DMRichR::annotateRegions(TxDb = TxDb, annoDb = annoDb) %>%
  dplyr::select(CpG.Island, CpG.Shore, CpG.Shelf, Open.Sea) %>% as.data.frame() # hypomethylated CpGs
yes_hypo <- sapply(EL_hypo_CpG,FUN = function(x){length(x[x=="Yes"])})
count_hypo <- data.frame(Count = yes_hypo, Percent = yes_hypo/length(EL_hypo_sigRegions))

EL_counts <- data.frame(CpG = rep(rownames(count_hyper),2),
                         Count = c(count_hyper[[1]], count_hypo[[1]]),
                         Percent = round(c(count_hyper[[2]], count_hypo[[2]]),2),
                         Direction = rep(c("Hypermethylated", "Hypomethylated"), each=4)) # initialize dataframe
write.table(EL_counts, file = "PDO_DMRichments/EL_CpG_counts.txt", quote = FALSE, sep = '\t ', row.names = F) # saves the results as a text file in the working directory

pdf(file = "PDO_DMRichments/EL_CpG_counts.pdf")
ggplot(EL_counts, aes(fill=CpG, y=Percent, x=Direction)) + 
  geom_bar(position="dodge", stat = "identity", color = "black")+
  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
  labs(x ="", y = "") + ggtitle("Early vs Late") +
  scale_x_discrete(labels=c("Hypermethylated", "Hypomethylated")) +
  theme_minimal() + 
  scale_fill_manual(values = c("forestgreen", "goldenrod2", "dodgerblue", "blue3"))
dev.off()

# CpG and genic enrichment testing for hN vs hT ----------------------------------------
DMRich <- function(x){
  dmrList[x] %>% 
    DMRichR::DMRichCpG(regions = NT_regions, genome = genome) %T>%
    openxlsx::write.xlsx(file = glue::glue("PDO_DMRichments/NT_{names(dmrList)[x]}_CpG_enrichments.xlsx")) %>% 
    DMRichR::DMRichPlot(type = "CpG") %>% 
    ggplot2::ggsave(glue::glue("PDO_DMRichments/NT_{names(dmrList)[x]}_CpG_enrichments.pdf"),
                    plot = ., width = 4, height = 3)
  dmrList[x] %>% 
    DMRichR::DMRichGenic(regions = NT_regions, TxDb = TxDb, annoDb = annoDb) %T>%
    openxlsx::write.xlsx(file = glue::glue("PDO_DMRichments/NT_{names(dmrList)[x]}_genic_enrichments.xlsx")) %>% 
    DMRichR::DMRichPlot(type = "genic") %>% 
    ggplot2::ggsave(glue::glue("PDO_DMRichments/NT_{names(dmrList)[x]}_genic_enrichments.pdf"), plot = ., width = 4, height = 4)
}
dmrList <- NT_sigRegions %>% DMRichR::dmrList()
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
    glue::glue("PDO_DMRichments/NT_{direction}_{type}_enrichments.xlsx")
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
                ggplot2::ggsave(glue::glue("PDO_DMRichments/NT_{type}_multi_plot.pdf"), plot = ., device = NULL,
                                height = dplyr::case_when(type == "genic" ~ 5, type == "CpG" ~ 3.5), width = 7)
            })

# Genic enrichment counts
NT_hyper_annot <- NT_hyper_sigRegions %>% 
  DMRichR::annotateRegions(TxDb = TxDb, annoDb = annoDb) %>%
  dplyr::select(annotation) %>% 
  table() %>% as.data.frame() # hypermethylated frequency table
NT_hyper_annot$Percent <- NT_hyper_annot[[2]]/sum(NT_hyper_annot[[2]]) # add percent

NT_hypo_annot <- NT_hypo_sigRegions %>% 
  DMRichR::annotateRegions(TxDb = TxDb, annoDb = annoDb) %>%
  dplyr::select(annotation) %>% table() %>% as.data.frame() # hypomethylated frquency table# Percent graphs
NT_hypo_annot$Percent <- NT_hypo_annot[[2]]/sum(NT_hypo_annot[[2]]) # add percent

NT_counts <- data.frame(Annotation = rep(levels(NT_hyper_annot[[1]]), 2),
                        Count = c(NT_hyper_annot[[2]], NT_hypo_annot[[2]]),
                        Percent = round(c(NT_hyper_annot[[3]], NT_hypo_annot[[3]]),2),
                        Direction = rep(c("Hypermethylated", "Hypomethylated"), each=7)) # initialize dataframe
NT_counts$Annotation <- factor(NT_counts$Annotation, 
                               levels=c("Promoter", "5' UTR", "Exon", "Intron", "3' UTR", "Downstream", "Distal Intergenic")) # reorder
write.table(NT_counts, file = "PDO_DMRichments/NT_genic_counts.txt", quote = FALSE, sep = '\t ', row.names = F) # saves the results as a text file in the working directory

pdf(file = "PDO_DMRichments/NT_genic_counts.pdf")
ggplot(NT_counts, aes(fill=Annotation, y=Percent, x=Direction)) + 
  geom_bar(position="dodge", stat = "identity", color = "black")+
  scale_y_continuous(labels = scales::percent, limits = c(0,.5)) +
  labs(x ="", y = "") + ggtitle("hN vs hT") +
  scale_x_discrete(labels=c("Hypermethylated", "Hypomethylated")) +
  theme_minimal() + 
  scale_fill_manual(values=wesanderson::wes_palette("Zissou1", n = 7, type = "continuous") %>% rev())
dev.off()
rm(NT_hyper_annot, NT_hypo_annot, NT_counts)

# CpG enrichment counts
NT_hyper_CpG <- NT_hyper_sigRegions %>% 
  DMRichR::annotateRegions(TxDb = TxDb, annoDb = annoDb) %>%
  dplyr::select(CpG.Island, CpG.Shore, CpG.Shelf, Open.Sea) %>% as.data.frame() # hypermethylated CpGs
yes_hyper <- sapply(NT_hyper_CpG,FUN = function(x){length(x[x=="Yes"])})
count_hyper <- data.frame(Count = yes_hyper, Percent = yes_hyper/length(NT_hyper_sigRegions))

NT_hypo_CpG <- NT_hypo_sigRegions %>% 
  DMRichR::annotateRegions(TxDb = TxDb, annoDb = annoDb) %>%
  dplyr::select(CpG.Island, CpG.Shore, CpG.Shelf, Open.Sea) %>% as.data.frame() # hypomethylated CpGs
yes_hypo <- sapply(NT_hypo_CpG,FUN = function(x){length(x[x=="Yes"])})
count_hypo <- data.frame(Count = yes_hypo, Percent = yes_hypo/length(NT_hypo_sigRegions))

NT_counts <- data.frame(CpG = rep(rownames(count_hyper),2),
                        Count = c(count_hyper[[1]], count_hypo[[1]]),
                        Percent = round(c(count_hyper[[2]], count_hypo[[2]]),2),
                        Direction = rep(c("Hypermethylated", "Hypomethylated"), each=4)) # initialize dataframe
write.table(NT_counts, file = "PDO_DMRichments/NT_CpG_counts.txt", quote = FALSE, sep = '\t ', row.names = F) # saves the results as a text file in the working directory

pdf(file = "PDO_DMRichments/NT_CpG_counts.pdf")
ggplot(NT_counts, aes(fill=CpG, y=Percent, x=Direction)) + 
  geom_bar(position="dodge", stat = "identity", color = "black")+
  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
  labs(x ="", y = "") + ggtitle("hN vs hT") +
  scale_x_discrete(labels=c("Hypermethylated", "Hypomethylated")) +
  theme_minimal() + 
  scale_fill_manual(values = c("forestgreen", "goldenrod2", "dodgerblue", "blue3"))
dev.off()
rm(NT_hyper_CpG, NT_hypo_CpG, yes_hyper, yes_hypo, count_hyper, count_hypo, NT_counts)

# Manhattan plots -------------------------------------------------
EL_regions %>% 
  DMRichR::annotateRegions(TxDb = TxDb, annoDb = annoDb) %>%
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>%
  sort() %>% dplyr::as_tibble() %>% dplyr::select(geneSymbol, seqnames, start, q.value)  %>% 
  dplyr::mutate(seqnames = substring(.$seqnames, 4)) %>%
  CMplot::CMplot(col = c("grey30", "grey60"), plot.type = "m", 
                 LOG10 = TRUE, ylim = NULL, threshold = 0.01, threshold.lty = 1, threshold.lwd = 1, 
                 threshold.col = "black", cex = 0.5, cex.axis = 0.7, amplify = FALSE, chr.den.col = c("darkgreen", "yellow", "red"), 
                 bin.size = 1e+06, file = "pdf", memo = "")
file.rename("Rect_Manhtn.q.value.pdf", "PDO_DMRs/EL_manhattan.pdf")

NT_regions %>% 
  DMRichR::annotateRegions(TxDb = TxDb, annoDb = annoDb) %>%
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>%
  sort() %>% dplyr::as_tibble() %>% dplyr::select(geneSymbol, seqnames, start, q.value)  %>% 
  dplyr::mutate(seqnames = substring(.$seqnames, 4)) %>%
  CMplot::CMplot(col = c("grey30", "grey60"), plot.type = "m", 
                 LOG10 = TRUE, ylim = NULL, threshold = 0.01, threshold.lty = 1, threshold.lwd = 1, 
                 threshold.col = "black", cex = 0.5, cex.axis = 0.7, amplify = FALSE, chr.den.col = c("darkgreen", "yellow", "red"), 
                 bin.size = 1e+06, file = "pdf", memo = "")
file.rename("Rect_Manhtn.q.value.pdf", "PDO_DMRs/NT_manhattan.pdf")

# Prepare HOMER -------------------------------------------------------------------
prepHOMER <- function (sigRegions = sigRegions, regions = regions, dir.name = dir.name) 
{
  dir.create(dir.name)
  sigRegions %>% DMRichR::gr2bed(paste(dir.name,"DMRs.bed", sep="/"))
  sigRegions %>% plyranges::filter(stat > 0) %>% DMRichR::gr2bed(paste(dir.name,"DMRs_hyper.bed", sep="/"))
  sigRegions %>% plyranges::filter(stat < 0) %>% DMRichR::gr2bed(paste(dir.name,"DMRs_hypo.bed", sep="/"))
  regions %>% DMRichR::gr2bed(paste(dir.name,"background.bed", sep="/"))
}

EL_sigRegions %>% prepHOMER(regions = EL_regions, dir.name = "PDO_DMRs/EL_HOMER")
NT_sigRegions %>% prepHOMER(regions = NT_regions, dir.name = "PDO_DMRs/NT_HOMER")

# ChromHMM and Reference Epigenomes ---------------------------------------
if(length(grep("genomecenter.ucdavis.edu", .libPaths())) > 0 & genome == "hg38"){
  LOLA <- function(x){
    dir.create(names(dmrList)[x])
    setwd(names(dmrList)[x])
    dmrList[x] %>%
      DMRichR::chromHMM(regions = regions,
                        cores = floor(cores/3)) %>% 
      DMRichR::chromHMM_heatmap()
    if(file.exists("Rplots.pdf")){file.remove("Rplots.pdf")}
  }
  
  cores = 3
  dir.create("LOLA_NT")
  setwd("LOLA_NT")
  dmrList <- NT_sigRegions %>% DMRichR::dmrList()
  regions = NT_regions
  parallel::mclapply(seq_along(dmrList), LOLA, mc.cores = 3, mc.silent = TRUE)
  setwd("..")
  
  dir.create("LOLA_EL")
  setwd("LOLA_EL")
  dmrList <- EL_sigRegions %>% DMRichR::dmrList()
  regions = EL_regions
  parallel::mclapply(seq_along(dmrList), LOLA, mc.cores = 3, mc.silent = TRUE)
  setwd("..")
}

# Machine learning --------------------------------------------------------
setwd("PDO_DMRs")
tryCatch({methylLearnOutput <- DMRichR::methylLearn(bs.filtered.bsseq = bs.filtered.bsseqPDO,
                                                    sigRegions = EL_sigRegions,
                                                    testCovariate = testCovariate,
                                                    TxDb = TxDb,
                                                    annoDb = annoDb,
                                                    topPercent = 1,
                                                    output = "all",
                                                    saveHtmlReport = TRUE)
dir.create("Machine_learning")
if(length(methylLearnOutput) == 1) {
  openxlsx::write.xlsx(list(Annotations_Common_DMRs = methylLearnOutput), 
                       file = "Machine_learning/Machine_learning_output_one.xlsx") 
} else {
  openxlsx::write.xlsx(list(Annotations_Common_DMRs = methylLearnOutput$`Annotated common DMRs`,
                            RF_Ranking_All_DMRs = methylLearnOutput$`RF ranking`,
                            SVM_Ranking_All_DMRs = methylLearnOutput$`SVM ranking`),
                       file = "Machine_learning/Machine_learning_output_all.xlsx") 
}
save(methylLearnOutput, file = "RData/machineLearning.RData")
#load("RData/machineLearing.RData")
},
error = function(error_condition) {
  print(glue::glue("Warning: methylLearn did not finish. \\
                      You may have not had enough top DMRs across algorithms."))
})
rm(methylLearnOutput)
file.rename("Machine_learning", "EL_Machine_Learning")

tryCatch({methylLearnOutput <- DMRichR::methylLearn(bs.filtered.bsseq = bs.filtered.bsseq,
                                                    sigRegions = NT_sigRegions,
                                                    testCovariate = testCovariate,
                                                    TxDb = TxDb,
                                                    annoDb = annoDb,
                                                    topPercent = 1,
                                                    output = "all",
                                                    saveHtmlReport = TRUE)
dir.create("Machine_learning")
if(length(methylLearnOutput) == 1) {
  openxlsx::write.xlsx(list(Annotations_Common_DMRs = methylLearnOutput), 
                       file = "Machine_learning/Machine_learning_output_one.xlsx") 
} else {
  openxlsx::write.xlsx(list(Annotations_Common_DMRs = methylLearnOutput$`Annotated common DMRs`,
                            RF_Ranking_All_DMRs = methylLearnOutput$`RF ranking`,
                            SVM_Ranking_All_DMRs = methylLearnOutput$`SVM ranking`),
                       file = "Machine_learning/Machine_learning_output_all.xlsx") 
}
save(methylLearnOutput, file = "RData/machineLearning.RData")
#load("RData/machineLearing.RData")
},
error = function(error_condition) {
  print(glue::glue("Warning: methylLearn did not finish. \\
                      You may have not had enough top DMRs across algorithms."))
})
rm(methylLearnOutput)
file.rename("Machine_learning", "NT_Machine_Learning")
