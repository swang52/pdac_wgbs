# Get Packages --------------------------------------------------------------
cat("\n[DMRichR] Initializing \t\t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")

if(length(grep("genomecenter.ucdavis.edu", .libPaths())) > 0){
  .libPaths("/share/lasallelab/programs/DMRichR/R_4.1")
  AnnotationHub::setAnnotationHubOption("CACHE", "/share/lasallelab/programs/DMRichR/R_4.1")
  ExperimentHub::setExperimentHubOption("CACHE", "/share/lasallelab/programs/DMRichR/R_4.1")
}

if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if(!requireNamespace("remotes", quietly = TRUE))
  install.packages("remotes")
if(suppressPackageStartupMessages(!requireNamespace("DMRichR", quietly = TRUE))){
  Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = TRUE)
  BiocManager::install("ben-laufer/DMRichR")
}
suppressPackageStartupMessages(library(DMRichR))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(RColorBrewer))

# Functions ----------------
heatmap <- function(bs.filtered.bsseq = bs.filtered.bsseq, sigRegions = sigRegions, testCovariate = testCovariate, filename = filename, colors = colors, ...) {
  bsseq::getMeth(BSseq = bs.filtered.bsseq, regions = sigRegions, type = "smooth", what = "perRegion") %>% na.omit() %>% as.matrix() %>% 
    pheatmap::pheatmap(., scale = "row", annotation_col = bs.filtered.bsseq %>% pData() %>% as.data.frame() %>%
                         dplyr::select_if(~nlevels(.) > 1), color = RColorBrewer::brewer.pal(11, name = "RdBu") %>% rev(), show_colnames = TRUE, 
                       border_color = "grey", main = glue::glue("Z-Scores of {length(sigRegions)} Differentially Methylated Regions"), 
                       fontsize = 8, filename = filename, width = 11, height = 8.5, cellwidth = 12, annotation_colors = colors, ...) %>% return()}

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

## Load and process samples ## ------------------
load("RData/bismark_PDOs.RData")
bs.filteredSub <- bs.filteredPDO[, which(bs.filteredPDO$Subtype != "ADEX")]
bs.filteredSub <- bs.filteredSub[, which(bs.filteredSub$Subtype != "Immunogenic")]
bs.filteredSub$Subtype <- droplevels(bs.filteredSub$Subtype) # drop normal from levels
bs.filteredSub$CombinedStage <- c()
save(bs.filteredSub, file = "RData/bismark_Sub.RData")
#load("RData/bismark_Sub.RData")
rm(bs.filteredPDO)

load("bsseq_PDOs.RData")
bs.filtered.bsseqSub <- bs.filtered.bsseqPDO[, which(bs.filtered.bsseqPDO$Subtype != "ADEX")]
bs.filtered.bsseqSub <- bs.filtered.bsseqSub[,which(bs.filtered.bsseqSub$Subtype != "Immunogenic")]
bs.filtered.bsseqSub$Subtype = droplevels(bs.filtered.bsseqSub$Subtype) # drop normal from levels
bs.filtered.bsseqSub$CombinedStage = c()
save(bs.filtered.bsseqSub, file = "RData/bsseq_sub.RData")
#load("RData/bsseq_sub.RData")
rm(bs.filtered.bsseqPDO)

load("RData/settings_NPDOs.RData")
testCovariate <- as.character("Subtype") # Test covariate 
DMRichR::annotationDatabases(genome = genome, EnsDb = EnsDb) # setup annotation databases

## Make global density plots ## ------------------
dir.create("Subtype_Global")
bs.filtered.bsseqSub %>% DMRichR::globalStats(genome = genome,
                                           testCovariate = testCovariate) %>%
  openxlsx::write.xlsx("Subtype_Global/smoothed_globalStats.xlsx") # Smoothened statistics

windows <- bs.filtered.bsseqSub %>% DMRichR::windows(goi = goi)
CpGs <- bs.filtered.bsseqSub %>% DMRichR::CpGs()
CGi <- bs.filtered.bsseqSub %>% DMRichR::CGi(genome = genome)  
plots <- c("windows", "CpGs", "CGi")
group =  bs.filtered.bsseqSub %>% pData() %>% dplyr::as_tibble() %>%
  dplyr::pull(!!"Subtype") %>% forcats::fct_rev()

plot = function (matrix = matrix, group = NA) {print(glue::glue("Density plot of {length(matrix)} sites"))
  matrix %>% dplyr::as_tibble() %>% magrittr::set_colnames(paste(group, seq_along(1:length(group)))) %>% dplyr::transmute(
    Group1 = dplyr::select(., dplyr::contains(levels(group)[1])) %>% rowMeans() * 100, 
    Group2 = dplyr::select(., dplyr::contains(levels(group)[2])) %>% rowMeans() * 100) %>% 
    magrittr::set_colnames(c(levels(group)[1],levels(group)[2])) %>% 
    tidyr::gather(key = "variable", value = "value") %>% dplyr::mutate(variable = factor(.$variable, 
                                                                                         levels = levels(group))) %>% ggplot(aes(value, color = variable)) + 
    geom_density(size = 1.2) + labs(x = "Percent Methylation", y = "Density", color = "Group") + theme_classic() +
    scale_color_manual(values = c("#7FC97F", "#BEAED4")) + 
    scale_x_continuous(expand = c(0.05, 0.05), breaks = c(0, 25, 50, 75, 100)) + scale_y_continuous(expand = c(0, 0.001)) + 
    theme(axis.text = element_text(size = 8), axis.title = element_text(size = 8), strip.text = element_text(size = 8), 
          legend.text = element_text(size = 8), legend.position = "bottom", legend.title = element_text(size = 8)) %>% return()
}

purrr::walk(plots,
            function(plotMatrix, group =  bs.filtered.bsseqSub %>% pData() %>% dplyr::as_tibble() %>%
                       dplyr::pull(!!"Subtype") %>% forcats::fct_rev()){
              title <- dplyr::case_when(plotMatrix == "windows" ~ "20Kb Windows", plotMatrix == "CpGs" ~ "Single CpG", plotMatrix == "CGi" ~ "CpG Island")
              plotMatrix %>% get() %>% plot(group = group) %>% 
                ggplot2::ggsave(glue::glue("Subtype_Global/{title} Density Plot.pdf"), plot = ., device = NULL, width = 11, height = 4)
              })
              
# DMRs ---------------------
sub_regions <- dmrseq::dmrseq(bs = bs.filteredSub, cutoff = cutoff, minNumRegion = minCpGs, maxPerms = maxPerms, testCovariate = testCovariate)
sub_regions <- sub_regions %>% 
  plyranges::mutate(direction = dplyr::case_when(stat > 0 ~ "Hypermethylated",
                                                 stat < 0 ~ "Hypomethylated"),
                    difference = round(beta/pi *100))

if(sum(sub_regions$qval < 0.05) < 100 & sum(sub_regions$pval < 0.05) != 0){
  sub_sigRegions <- sub_regions %>%
    plyranges::filter(pval < 0.05)
}else if(sum(sub_regions$qval < 0.05) >= 100){
  sub_sigRegions <- sub_regions %>%
    plyranges::filter(qval < 0.05)
}else if(sum(sub_regions$pval < 0.05) == 0){
  stop(glue::glue("No significant DMRs detected in {length(regions)} background regions"))
}
sub_sigRegions <- sub_sigRegions %>% plyranges::filter(seqnames != "chrX") # remove X chromosome DMRs
sub_sigRegions <- sub_sigRegions %>% plyranges::filter(seqnames != "chrY") # remove Y chromosome DMRs

dir.create("Subtype_DMRs")
gr2bed(sub_sigRegions, "Subtype_DMRs/sub_DMRs.bed")
gr2bed(sub_regions, "Subtype_DMRs/sub_backgroundRegions.bed")
save(sub_regions, sub_sigRegions, file = "RData/sub_DMRs.RData")
#load("RData/sub_DMRs.RData")

sub_sigRegions %>%
  DMRichR::annotateRegions(TxDb = TxDb, annoDb = annoDb) %T>%
  DMRichR::DMReport(regions = sub_regions, bs.filtered = bs.filteredSub, coverage = coverage, name = "Subtype_DMRs/sub_DMReport") %>% 
  openxlsx::write.xlsx(file = "Subtype_DMRs/sub_DMRs_annotated.xlsx") # annotate DMRs

print(glue::glue("Extracting individual smoothed methylation values of DMRs..."))
bs.filtered.bsseqSub %>%
  DMRichR::smooth2txt(regions = sub_sigRegions,
                      txt = "Subtype_DMRs/sub_DMR_individual_smoothed_methylation.txt")

print(glue::glue("Extracting individual smoothed methylation values of background regions..."))
bs.filtered.bsseqSub %>%
  DMRichR::smooth2txt(regions = sub_regions,
                      txt = "Subtype_DMRs/sub_background_region_individual_smoothed_methylation.txt")

colors <- list(c("#4DAF4A", "#E41A1C", "#FF7F00", "#984EA3") %>% 
                 setNames(bs.filtered.bsseqSub %>% pData() %>% as.data.frame() %>% purrr::pluck("Stage") %>% 
                            unique() %>% sort() %>% rev()),
               c("#7FC97F", "#BEAED4") %>% setNames(bs.filtered.bsseqSub %>% pData() %>% as.data.frame() %>% purrr::pluck("Subtype") %>% 
                                                                            unique() %>% sort() %>% rev()))
names(colors) = c("Stage", "Subtype")
sub_sigRegions %>% heatmap(bs.filtered.bsseq = bs.filtered.bsseqSub, testCovariate = testCovariate, 
                           filename = "Subtype_DMRs/sub_heatmap.pdf", 
                           colors = colors)

# CpG and genic enrichment testing for Subtype  ----------------------------------------
DMRich <- function(x){
  dmrList[x] %>% 
    DMRichR::DMRichCpG(regions = sub_regions, genome = genome) %T>%
    openxlsx::write.xlsx(file = glue::glue("DMRichments/sub_{names(dmrList)[x]}_CpG_enrichments.xlsx")) %>% 
    DMRichR::DMRichPlot(type = "CpG") %>% 
    ggplot2::ggsave(glue::glue("DMRichments/sub_{names(dmrList)[x]}_CpG_enrichments.pdf"),
                    plot = ., width = 4, height = 3)
  dmrList[x] %>% 
    DMRichR::DMRichGenic(regions = sub_regions, TxDb = TxDb, annoDb = annoDb) %T>%
    openxlsx::write.xlsx(file = glue::glue("DMRichments/sub_{names(dmrList)[x]}_genic_enrichments.xlsx")) %>% 
    DMRichR::DMRichPlot(type = "genic") %>% 
    ggplot2::ggsave(glue::glue("DMRichments/sub_{names(dmrList)[x]}_genic_enrichments.pdf"), plot = ., width = 4, height = 4)
}
dmrList <- sub_sigRegions %>% DMRichR::dmrList()
parallel::mclapply(seq_along(dmrList), DMRich, mc.cores = 1, mc.silent = TRUE)

genicCount(sigRegions = sub_sigRegions, project = "sub")
cpgCount(sigRegions = sub_sigRegions, project = "sub")

# Manhattan plots -------------------------------------------------
sub_regions %>% 
  DMRichR::annotateRegions(TxDb = TxDb, annoDb = annoDb) %>%
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>%
  sort() %>% dplyr::as_tibble() %>% dplyr::select(geneSymbol, seqnames, start, p.value)  %>% 
  dplyr::mutate(seqnames = substring(.$seqnames, 4)) %>%
  CMplot::CMplot(col = c("grey30", "grey60"), plot.type = "m", 
                 LOG10 = TRUE, ylim = NULL, threshold = 0.05, threshold.lty = 1, threshold.lwd = 1, 
                 threshold.col = "black", cex = 0.5, cex.axis = 0.7, amplify = FALSE, chr.den.col = c("darkgreen", "yellow", "red"), 
                 bin.size = 1e+06, file = "pdf")
file.rename("Rect_Manhtn.p.value.pdf", "Subtype_DMRs/Manhattan.pdf")

# Prepare HOMER -------------------------------------------------------------------
prepHOMER <- function (sigRegions = sigRegions, regions = regions, dir.name = dir.name) 
{
  dir.create(dir.name)
  sigRegions %>% DMRichR::gr2bed(paste(dir.name,"DMRs.bed", sep="/"))
  sigRegions %>% plyranges::filter(stat > 0) %>% DMRichR::gr2bed(paste(dir.name,"DMRs_hyper.bed", sep="/"))
  sigRegions %>% plyranges::filter(stat < 0) %>% DMRichR::gr2bed(paste(dir.name,"DMRs_hypo.bed", sep="/"))
  regions %>% DMRichR::gr2bed(paste(dir.name,"background.bed", sep="/"))
}

sub_sigRegions %>% prepHOMER(regions = sub_regions, dir.name = "Subtype_DMRs/Sub_HOMER")

# Machine Learning -------------------------------------------------------------------
sub_sigRegions <- sub_regions %>% plyranges::filter(pval < 0.05) %>% plyranges::filter(seqnames != "chrX") %>% plyranges::filter(seqnames != "chrY")

dir.create("Sub")
setwd("Sub")
tryCatch({methylLearnOutput <- DMRichR::methylLearn(bs.filtered.bsseq = bs.filtered.bsseqSub,
                                                    sigRegions = sub_sigRegions,
                                                    testCovariate = "Subtype",
                                                    TxDb = TxDb,
                                                    annoDb = annoDb,
                                                    topPercent = 1,
                                                    output = "all",
                                                    saveHtmlReport = TRUE)
if(length(methylLearnOutput) == 1) {
  openxlsx::write.xlsx(list(Annotations_Common_DMRs = methylLearnOutput), 
                       file = "./Machine_learning_output_one.xlsx") 
} else {
  openxlsx::write.xlsx(list(Annotations_Common_DMRs = methylLearnOutput$`Annotated common DMRs`,
                            RF_Ranking_All_DMRs = methylLearnOutput$`RF ranking`,
                            SVM_Ranking_All_DMRs = methylLearnOutput$`SVM ranking`),
                       file = "./Machine_learning_output_all.xlsx") 
}
save(methylLearnOutput, file = "Sub_machineLearning.RData")
#load("machineLearing.RData")
},
error = function(error_condition) {
  print(glue::glue("Warning: methylLearn did not finish. \\
                      You may have not had enough top DMRs across algorithms."))
})
rm(methylLearnOutput)

# DMR plots -----------------
annoTrack <- dmrseq::getAnnot("hg38")
avg_bs = collapseBSseq(bs.filtered, group = c("A", "B", "A", "A", "B", "A", "A", "B"))
avg_bs = collapseBSseq(bs.filtered, group = c("A", "A", "A", "A", "A", "A", "A", "A", 
                                              "B", "B", "B", "A", "A", "B", "A", "B", "A", "B"))
gene_name = "BMP4"
coor = data.frame(chr = c("chr14"), start = c(53949736),
                   end = c(53956825)) %>% makeGRangesFromDataFrame(.)
DMRs = sub_sigRegions %>%
  DMRichR::annotateRegions(TxDb = TxDb, annoDb = annoDb) %>%
  plyranges::filter(geneSymbol == gene_name) %>%
  makeGRangesFromDataFrame(., keep.extra.columns = TRUE)
pdf(paste(gene_name,".pdf", sep = ""), height = 4, width = 8)
dmrseq::plotDMRs(avg_bs_n, regions = coor, testCovariate = "Subtype",
                 extend = 1000, main = gene_name,
                 annoTrack = annoTrack, addRegions = DMRs,
                 regionCol = "#FF00001A",lwd = 1,
                 addPoints = FALSE, qval = FALSE, stat = FALSE)
dev.off()
