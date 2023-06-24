# Get Packages
library(DMRichR)

# Make GRanges objects for mNPTM -------------
cat("\n[DMRichR] Making GRanges objects \t\t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
hyper_MPN = read.table("MPN_hyper_hg38.bed", sep = "\t", col.names = c("chr", "start", "end"))
hypo_MPN = read.table("MPN_hypo_hg38.bed", sep = "\t", col.names = c("chr", "start", "end"))
hyper_TM = read.table("TM_hyper_hg38.bed", sep = "\t", col.names = c("chr", "start", "end"))
hypo_TM = read.table("TM_hypo_hg38.bed", sep = "\t", col.names = c("chr", "start", "end"))
TM_background = read.table("TM_bg_hg38.bed", sep = "\t", col.names = c("chr", "start", "end"))
MPN_background = read.table("MPN_bg_hg38.bed", sep = "\t", col.names = c("chr", "start", "end"))

hyper_MPN = makeGRangesFromDataFrame(hyper_MPN) # 4753 regions
hyper_TM = makeGRangesFromDataFrame(hyper_TM) # 2392 regions
hypo_MPN = makeGRangesFromDataFrame(hypo_MPN) # 12520 regions
hypo_TM = makeGRangesFromDataFrame(hypo_TM) # 1379 regions
TM = c(hyper_TM, hypo_TM)
MPN = c(hyper_MPN, hypo_MPN)
TM_background = makeGRangesFromDataFrame(TM_background)
MPN_background = makeGRangesFromDataFrame(MPN_background)

# Load DMRs for PDOs
load("RData/NT_DMRs.RData")
load("RData/EL_DMRs.RData")

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
}