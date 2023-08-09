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
                       fontsize = 10, filename = filename, width = 11, height = 8.5, cellwidth = 12, annotation_colors = colors, ...) %>% return()}

genicCount <- function(sigRegions = sigRegions, project = c("TM", "MPN", "EL", "NT", "sub"), ...){
  hyper <- sigRegions %>% plyranges::filter(stat > 0) %>% 
    DMRichR::annotateRegions(TxDb = TxDb, annoDb = annoDb) %>%
    dplyr::select(annotation) %>% table() %>% as.data.frame() %>%
    dplyr::mutate(Percent = Freq/sum(Freq)) %>%
    stats::setNames(c("Annotation","Frequency","Percentage"))
  hypo <- sigRegions %>% plyranges::filter(stat < 0) %>% 
    DMRichR::annotateRegions(TxDb = TxDb, annoDb = annoDb) %>%
    dplyr::select(annotation) %>% table() %>% as.data.frame() %>% 
    dplyr::mutate(Percent = Freq/sum(Freq)) %>%
    stats::setNames(c("Annotation","Frequency","Percentage"))
  data.frame(Annotation = c(hyper$Annotation, hypo$Annotation),
             Count = c(hyper$Freq, hypo$Freq),
             Percent = round(c(hyper$Percent, hypo$Percent),2),
             Direction = c(rep("Hypermethylated", nrow(hyper)), rep("Hypomethylated", nrow(hypo)))) %>% 
    dplyr::mutate(Annotation = stringr::str_replace(Annotation, "3' UTR", "3UTR")) %>%
    dplyr::mutate(Annotation = stringr::str_replace(Annotation, "5' UTR", "5UTR")) %>%
    write.table(., file = glue::glue("DMRichments/{project}_genic_counts.txt"), quote = FALSE, sep = '\t ', row.names = F)
}

cpgCount <- function(sigRegions = sigRegions, project = c("TM", "MPN", "EL", "NT", "sub"), ...){
  hyper <- sigRegions %>% plyranges::filter(stat > 0) %>% 
    DMRichR::annotateRegions(TxDb = TxDb, annoDb = annoDb) %>%
    dplyr::select(CpG.Island, CpG.Shore, CpG.Shelf, Open.Sea) %>% as.data.frame() %>%
    sapply(., FUN = function(x){length(x[x=="Yes"])}) %>%
    data.frame(Count = ., Percent = ./length(sigRegions%>% plyranges::filter(stat > 0)))
  hypo <- sigRegions %>% plyranges::filter(stat < 0) %>% 
    DMRichR::annotateRegions(TxDb = TxDb, annoDb = annoDb) %>%
    dplyr::select(CpG.Island, CpG.Shore, CpG.Shelf, Open.Sea) %>% as.data.frame() %>%
    sapply(., FUN = function(x){length(x[x=="Yes"])}) %>%
    data.frame(Count = ., Percent = ./length(sigRegions%>% plyranges::filter(stat < 0)))
  data.frame(CpG = rep(rownames(hyper),2),
             Count = c(hyper$Count, hypo$Count),
             Percent = round(c(hyper$Percent, hypo$Percent),2),
             Direction = rep(c("Hypermethylated", "Hypomethylated"), each = 4)) %>%
    write.table(., file = glue::glue("DMRichments/{project}_CpG_counts.txt"), quote = FALSE, sep = '\t ', row.names = F)
}
# Setup annotation databases and load files ----------------------------------------------
load("RData/bismark_NPDOs.RData") # methylation values
load("RData/bsseq_NPDOs.RData") # smoothened methylation values
load("RData/settings_NPDOs.RData")
DMRichR::annotationDatabases(genome = genome, EnsDb = EnsDb)

# Normal vs Tumor DMRs ---------------------
NT_regions <- dmrseq::dmrseq(bs = bs.filtered, cutoff = cutoff, minNumRegion = minCpGs, maxPerms = maxPerms, testCovariate = testCovariate)
NT_regions <- NT_regions %>% 
  plyranges::mutate(direction = dplyr::case_when(stat > 0 ~ "Hypermethylated",
                                                 stat < 0 ~ "Hypomethylated"),
                    difference = round(beta/pi *100))

if(sum(NT_regions$qval < 0.05) < 100 & sum(NT_regions$pval < 0.01) != 0){
  NT_sigRegions <- NT_regions %>%
    plyranges::filter(pval < 0.05)
}else if(sum(NT_regions$qval < 0.05) >= 100){
  NT_sigRegions <- NT_regions %>%
    plyranges::filter(qval < 0.05)
}else if(sum(regions$pval < 0.05) == 0){
  stop(glue::glue("No significant DMRs detected in {length(regions)} background regions"))
}
NT_sigRegions <- NT_sigRegions %>% plyranges::filter(seqnames != "chrX")
NT_sigRegions <- NT_sigRegions %>% plyranges::filter(seqnames != "chrY")

dir.create("PDO_DMRs")
gr2bed(NT_sigRegions, "PDO_DMRs/NT_DMRs.bed")
gr2bed(NT_regions, "PDO_DMRs/NT_backgroundRegions.bed")
save(NT_regions, NT_sigRegions, file = "RData/NT_DMRs.RData")
#load("RData/NT_DMRs.RData")

NT_sigRegions %>%
  DMRichR::annotateRegions(TxDb = TxDb, annoDb = annoDb) %T>%
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

# early vs late DMRs ---------------------
# Modify bs.filtered
bs.filteredPDO <- bs.filtered[, which(bs.filtered$Stage != "Normal")]
bs.filteredPDO$Stage <- droplevels(bs.filteredPDO$Stage) # drop normal from levels
bs.filteredPDO$CombinedStage = bs.filteredPDO$Stage %>% as.data.frame() %>% 
  with(., ifelse(. == "Metastatic" | . == "Locally advanced", "Late", "Early")) %>% as.factor()
save(bs.filteredPDO, file = "RData/bismark_PDOs.RData")
#load("RData/bismark_PDOs.RData")

EL_regions <- dmrseq::dmrseq(bs = bs.filteredPDO, cutoff = cutoff, minNumRegion = minCpGs, maxPerms = maxPerms, testCovariate = testCovariate)
EL_regions <- EL_regions %>% 
  plyranges::mutate(direction = dplyr::case_when(stat > 0 ~ "Hypermethylated",
                                                 stat < 0 ~ "Hypomethylated"),
                    difference = round(beta/pi *100))

if(sum(EL_regions$qval < 0.05) < 100 & sum(EL_regions$pval < 0.01) != 0){
  EL_sigRegions <- EL_regions %>%
    plyranges::filter(pval < 0.05)
}else if(sum(EL_regions$qval < 0.05) >= 100){
  EL_sigRegions <- EL_regions %>%
    plyranges::filter(qval < 0.05)
}else if(sum(regions$pval < 0.05) == 0){
  stop(glue::glue("No significant DMRs detected in {length(regions)} background regions"))
}
EL_sigRegions <- EL_sigRegions %>% plyranges::filter(seqnames != "chrX")
EL_sigRegions <- EL_sigRegions %>% plyranges::filter(seqnames != "chrY")

gr2bed(EL_sigRegions, "PDO_DMRs/EL_DMRs.bed")
gr2bed(EL_regions, "PDO_DMRs/EL_backgroundRegions.bed")
save(EL_regions, EL_sigRegions, file = "RData/EL_DMRs.RData")
#load("RData/EL_DMRs.RData")

EL_sigRegions %>%
  DMRichR::annotateRegions(TxDb = TxDb, annoDb = annoDb) %T>%
  DMRichR::DMReport(regions = EL_regions, bs.filtered = bs.filteredPDO, coverage = coverage, name = "PDO_DMRs/EL_DMReport") %>% 
  openxlsx::write.xlsx(file = "PDO_DMRs/EL_DMRs_annotated.xlsx") # annotate DMRs

bs.filtered.bsseqPDO <- bs.filtered.bsseq[, which(bs.filtered.bsseq$Stage != "Normal")]
bs.filtered.bsseqPDO$Stage = droplevels(bs.filtered.bsseqPDO$Stage) # drop normal from levels
bs.filtered.bsseqPDO$CombinedStage = bs.filtered.bsseqPDO$Stage %>% as.data.frame() %>% 
  with(., ifelse(. == "Metastatic" | . == "Locally advanced", "Late", "Early")) %>% as.factor()
save(bs.filtered.bsseqPDO, file = "RData/bsseq_PDOs.RData")
#load("RData/bsseq_PDOs.RData")

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

# Identify overlaps ---------------
EL_hyper_sigRegions <- EL_sigRegions %>% plyranges::filter(stat > 0)
NT_hyper_sigRegions <- NT_sigRegions %>% plyranges::filter(stat > 0)
EL_hypo_sigRegions <- EL_sigRegions %>% plyranges::filter(stat < 0)
NT_hypo_sigRegions <- NT_sigRegions %>% plyranges::filter(stat < 0)
hyper_overlap <- GenomicRanges::intersect(NT_hyper_sigRegions, EL_hyper_sigRegions, ignore.strand = T)
hypo_overlap <- GenomicRanges::intersect(NT_hypo_sigRegions, EL_hypo_sigRegions, ignore.strand = T)
overlap <- length(hyper_overlap) + length(hypo_overlap)
ELonly <- length(EL_sigRegions) - overlap
NTonly <- length(NT_sigRegions) - overlap
fit = euler(c("Early vs Late" = ELonly, "hN vs hT" = NTonly, 
              "Early vs Late&hN vs hT" = overlap))
pdf(file = "PDO_DMRs/euler.pdf")
plot(fit, quantities = TRUE, legend = list(lables = c("Early vs Late", "hN vs hT")),
     fills = list(fill = c("#66D2D6", "#E56997"), alpha = 0.8))
dev.off()

# CpG and genic enrichment testing for hN vs hT ----------------------------------------
dir.create("DMRichments")
DMRich <- function(x){
  dmrList[x] %>% 
    DMRichR::DMRichCpG(regions = NT_regions, genome = genome) %T>%
    openxlsx::write.xlsx(file = glue::glue("DMRichments/NT_{names(dmrList)[x]}_CpG_enrichments.xlsx")) %>% 
    DMRichR::DMRichPlot(type = "CpG") %>% 
    ggplot2::ggsave(glue::glue("DMRichments/NT_{names(dmrList)[x]}_CpG_enrichments.pdf"),
                    plot = ., width = 4, height = 3)
  dmrList[x] %>% 
    DMRichR::DMRichGenic(regions = NT_regions, TxDb = TxDb, annoDb = annoDb) %T>%
    openxlsx::write.xlsx(file = glue::glue("DMRichments/NT_{names(dmrList)[x]}_genic_enrichments.xlsx")) %>% 
    DMRichR::DMRichPlot(type = "genic") %>% 
    ggplot2::ggsave(glue::glue("DMRichments/NT_{names(dmrList)[x]}_genic_enrichments.pdf"), plot = ., width = 4, height = 4)
}
dmrList <- NT_sigRegions %>% DMRichR::dmrList()
parallel::mclapply(seq_along(dmrList), DMRich, mc.cores = 1, mc.silent = TRUE)

genicCount(sigRegions = NT_sigRegions, project = "NT")
cpgCount(sigRegions = NT_sigRegions, project = "NT")

# CpG and genic enrichment testing for Early vs Late ----------------------------------------
DMRich <- function(x){
  dmrList[x] %>% 
    DMRichR::DMRichCpG(regions = EL_regions, genome = genome) %T>%
    openxlsx::write.xlsx(file = glue::glue("DMRichments/EL_{names(dmrList)[x]}_CpG_enrichments.xlsx")) %>% 
    DMRichR::DMRichPlot(type = "CpG") %>% 
    ggplot2::ggsave(glue::glue("DMRichments/EL_{names(dmrList)[x]}_CpG_enrichments.pdf"),
                    plot = ., width = 4, height = 3)
  dmrList[x] %>% 
    DMRichR::DMRichGenic(regions = EL_regions, TxDb = TxDb, annoDb = annoDb) %T>%
    openxlsx::write.xlsx(file = glue::glue("DMRichments/EL_{names(dmrList)[x]}_genic_enrichments.xlsx")) %>% 
    DMRichR::DMRichPlot(type = "genic") %>% 
    ggplot2::ggsave(glue::glue("DMRichments/EL_{names(dmrList)[x]}_genic_enrichments.pdf"), plot = ., width = 4, height = 4)
}
dmrList <- EL_sigRegions %>% DMRichR::dmrList()
parallel::mclapply(seq_along(dmrList), DMRich, mc.cores = 1, mc.silent = TRUE)

genicCount(sigRegions = EL_sigRegions, project = "EL")
cpgCount(sigRegions = EL_sigRegions, project = "EL")

# Manhattan plots -------------------------------------------------
EL_regions %>% 
  DMRichR::annotateRegions(TxDb = TxDb, annoDb = annoDb) %>%
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>%
  sort() %>% dplyr::as_tibble() %>% dplyr::select(geneSymbol, seqnames, start, q.value)  %>% 
  dplyr::mutate(seqnames = substring(.$seqnames, 4)) %>%
  CMplot::CMplot(col = c("grey30", "grey60"), plot.type = "m", 
                 LOG10 = TRUE, ylim = NULL, threshold = 0.05, threshold.lty = 1, threshold.lwd = 1, 
                 threshold.col = "black", cex = 0.5, cex.axis = 0.7, amplify = FALSE, chr.den.col = c("darkgreen", "yellow", "red"), 
                 bin.size = 1e+06, file = "pdf", memo = "")
file.rename("Rect-Manhtn.q.value.pdf", "PDO_DMRs/EL_manhattan.pdf")

NT_regions %>% 
  DMRichR::annotateRegions(TxDb = TxDb, annoDb = annoDb) %>%
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>%
  sort() %>% dplyr::as_tibble() %>% dplyr::select(geneSymbol, seqnames, start, q.value)  %>% 
  dplyr::mutate(seqnames = substring(.$seqnames, 4)) %>%
  CMplot::CMplot(col = c("grey30", "grey60"), plot.type = "m", 
                 LOG10 = TRUE, ylim = NULL, threshold = 0.05, threshold.lty = 1, threshold.lwd = 1, 
                 threshold.col = "black", cex = 0.5, cex.axis = 0.7, amplify = FALSE, chr.den.col = c("darkgreen", "yellow", "red"), 
                 bin.size = 1e+06, file = "pdf", memo = "")
file.rename("Rect-Manhtn.q.value.pdf", "PDO_DMRs/NT_manhattan.pdf")

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

# Machine learning --------------------------------------------------------
dir.create("NT_Machine_learning")
setwd("NT_Machine_learning")
regions <- NT_regions %>% plyranges::filter(pval < 0.01)
tryCatch({methylLearnOutput <- DMRichR::methylLearn(bs.filtered.bsseq = bs.filtered.bsseq,
                                                    sigRegions = regions,
                                                    testCovariate = testCovariate,
                                                    TxDb = TxDb,
                                                    annoDb = annoDb,
                                                    topPercent = 1,
                                                    output = "all",
                                                    saveHtmlReport = TRUE)
if(length(methylLearnOutput) == 1) {
  openxlsx::write.xlsx(list(Annotations_Common_DMRs = methylLearnOutput), 
                       file = "./Machine_learning/Machine_learning_output_one.xlsx") 
} else {
  openxlsx::write.xlsx(list(Annotations_Common_DMRs = methylLearnOutput$`Annotated common DMRs`,
                            RF_Ranking_All_DMRs = methylLearnOutput$`RF ranking`,
                            SVM_Ranking_All_DMRs = methylLearnOutput$`SVM ranking`),
                       file = "./Machine_learning/Machine_learning_output_all.xlsx") 
}
save(methylLearnOutput, file = "RData/NT_machineLearning.RData")
#load("RData/machineLearing.RData")
},
error = function(error_condition) {
  print(glue::glue("Warning: methylLearn did not finish. \\
                      You may have not had enough top DMRs across algorithms."))
})
rm(methylLearnOutput)

dir.create("EL_Machine_learning")
setwd("EL_Machine_learning")
regions <- EL_regions %>% plyranges::filter(pval < 0.01)
tryCatch({methylLearnOutput <- DMRichR::methylLearn(bs.filtered.bsseq = bs.filtered.bsseq,
                                                    sigRegions = regions,
                                                    testCovariate = testCovariate,
                                                    TxDb = TxDb,
                                                    annoDb = annoDb,
                                                    topPercent = 1,
                                                    output = "all",
                                                    saveHtmlReport = TRUE)
if(length(methylLearnOutput) == 1) {
  openxlsx::write.xlsx(list(Annotations_Common_DMRs = methylLearnOutput), 
                       file = "./Machine_learning/Machine_learning_output_one.xlsx") 
} else {
  openxlsx::write.xlsx(list(Annotations_Common_DMRs = methylLearnOutput$`Annotated common DMRs`,
                            RF_Ranking_All_DMRs = methylLearnOutput$`RF ranking`,
                            SVM_Ranking_All_DMRs = methylLearnOutput$`SVM ranking`),
                       file = "./Machine_learning/Machine_learning_output_all.xlsx") 
}
save(methylLearnOutput, file = "RData/EL_machineLearning.RData")
#load("RData/machineLearing.RData")
},
error = function(error_condition) {
  print(glue::glue("Warning: methylLearn did not finish. \\
                      You may have not had enough top DMRs across algorithms."))
})
rm(methylLearnOutput)