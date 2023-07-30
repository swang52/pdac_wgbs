# Description -------------
# Perform ChromHMM enrichment on all DMR sets
# Requires bed files of mT vs mM and mN/mP vs mM DMRs liftover to hg38

# Get Packages --------------------------------------------------------------
cat("\n[DMRichR] Initializing \t\t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")

if(length(grep("genomecenter.ucdavis.edu", .libPaths())) > 0){
  .libPaths("/share/hwanglab/programs/DMRichR/R_4.1")
  AnnotationHub::setAnnotationHubOption("CACHE", "/share/hwanglab/programs/DMRichR/R_4.1")
  ExperimentHub::setExperimentHubOption("CACHE", "/share/hwanglab/programs/DMRichR/R_4.1")
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

# Make GRanges objects for mNPTM -------------
cat("\n[DMRichR] Making GRanges objects \t\t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
hyper_MPN = read.table("MPN_hyper_hg38.bed", sep = "\t", col.names = c("chr", "start", "end"))
hypo_MPN = read.table("MPN_hypo_hg38.bed", sep = "\t", col.names = c("chr", "start", "end"))
hyper_TM = read.table("TM_hyper_hg38.bed", sep = "\t", col.names = c("chr", "start", "end"))
hypo_TM = read.table("TM_hypo_hg38.bed", sep = "\t", col.names = c("chr", "start", "end"))
TM_background = read.table("TM_bg_hg38.bed", sep = "\t", col.names = c("chr", "start", "end"))
MPN_background = read.table("MPN_bg_hg38.bed", sep = "\t", col.names = c("chr", "start", "end"))

hyper_MPN = makeGRangesFromDataFrame(hyper_MPN)
hyper_TM = makeGRangesFromDataFrame(hyper_TM)
hypo_MPN = makeGRangesFromDataFrame(hypo_MPN)
hypo_TM = makeGRangesFromDataFrame(hypo_TM)
TM = c(hyper_TM, hypo_TM)
MPN = c(hyper_MPN, hypo_MPN)
TM_background = makeGRangesFromDataFrame(TM_background)
MPN_background = makeGRangesFromDataFrame(MPN_background)

# Get PDO DMRs
load("RData/NT_DMRs.RData")
load("RData/EL_DMRs.RData")
load("RData/sub_DMRs.RData")

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
  dir.create("LOLA_MPN")
  setwd("LOLA_MPN")
  dmrList <- GenomicRanges::GRangesList(`All DMRs` = MPN, `Hypermethylated DMRs` = hyper_MPN, `Hypomethylated DMRs` = hypo_MPN)
  regions = MPN_background
  parallel::mclapply(seq_along(dmrList), LOLA, mc.cores = 3, mc.silent = TRUE)
  setwd("..")
  
  dir.create("LOLA_TM")
  setwd("LOLA_TM")
  dmrList <- GenomicRanges::GRangesList(`All DMRs` = TM, `Hypermethylated DMRs` = hyper_TM, `Hypomethylated DMRs` = hypo_TM)
  regions = TM_background
  parallel::mclapply(seq_along(dmrList), LOLA, mc.cores = 3, mc.silent = TRUE)
  setwd("..")

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
  
  dir.create("LOLA_Subtype")
  setwd("LOLA_Subtype")
  dmrList <- sub_sigRegions %>% DMRichR::dmrList()
  regions = sub_regions
  parallel::mclapply(seq_along(dmrList), LOLA, mc.cores = 3, mc.silent = TRUE)
  setwd("..")
}