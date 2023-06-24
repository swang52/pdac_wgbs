# Description -------------
# Processes cytosine reports, calculates smoothened individual methylation values, and generates global methylation methylation 
# plots/statistics (i.e. PCA, global density plots, average CpG methylation statistics) for PDOs and murine organoids

# Get Packages --------------------------------------------------------------
# Package names
packages <- c("ggplot2", "DMRichR", "RColorBrewer")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Load packages
invisible(lapply(packages, library, character.only = TRUE))

# Make mNPTM PCA and global density plots -----------------------------------------------------------

## Global variables ##
genome <- as.character("mm9") # Options: hg38, hg19, mm10, mm9, rheMac10, rheMac8, rn6, danRer11, galGal6, bosTau9, panTro6, dm6, susScr11, canFam3, TAIR10, or TAIR9)
coverage <- as.integer(1) # CpG coverage cutoff minimum value is 1
perGroup <- as.double(1) # Percent of samples in all combinations of covariates meeting CpG coverage cutoff; Options: 0-1
minCpGs <- as.integer(5) # Minimum number of CpGs for a DMR
maxPerms <- as.integer(8) # Maximum number of permutations for the DMR analysis; no more than the # of samples
cutoff <- as.double(0.1) # Cutoff value [from 0 to 1] for the single CpG coefficient utilized to discover testable background regions
testCovariate <- as.character("Stage") # Test covariate 
EnsDb <- FALSE

## Load and process samples ##
bs.filtered <- processBismark(files = list.files(path = getwd(), pattern = "*.txt.gz"),
                              meta = openxlsx::read.xlsx("sample_info.xlsx",colNames = TRUE) %>% dplyr::mutate_if(is.character, as.factor),
                              testCovariate = testCovariate, adjustCovariate = adjustCovariate, matchCovariate = matchCovariate,
                              coverage = coverage, cores = cores, perGroup = perGroup, sexCheck = sexCheck)
bs.filtered = chrSelectBSseq(bs.filtered, seqnames = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11",
                                                       "chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19"), order = TRUE) # removes sex chromosomes
save(bs.filtered, file = "RData/bismark_mNPTM.RData")
#load("RData/bismark_mNPTM.RData")

bs.filtered.bsseq <- bsseq::BSmooth(bs.filtered, BPPARAM = BiocParallel::SnowParam(workers = 1, progressbar = TRUE)) # get individual smoothened methylation values
save(bs.filtered.bsseq,file = "RData/bsseq_mNPTM.RData")
#load("RData/bsseq_mNPTM.RData")

DMRichR::annotationDatabases(genome = genome, EnsDb = EnsDb) # setup annotation databases
settings_env <- ls(all = TRUE)
save(list = settings_env, file = "RData/settings_mNPTM.RData")
#load("RData/settings_mNPTM.RData")

## Perform PCA ##
group =  bs.filtered.bsseq %>% pData() %>% dplyr::as_tibble() %>%
  dplyr::pull(!!"Stage") %>% forcats::fct_rev()
beta = t(bs.filtered.bsseq %>% DMRichR::CpGs()) # get methylation values
beta = beta[,order(apply(beta,2,var))]
pca = prcomp(beta) # conducts principle component analysis
PCbeta = data.frame(pca$x, Group=group)
pov1 = summary(pca)$importance[2,1] # PC1 proportion of variance
pov1 = round(pov1*100,digits=2)
pov2 = summary(pca)$importance[2,2] # PC2 proportion of variance
pov2 = round(pov2*100,digits=2)

## Create PCA Plot ##
dir.create("mNPTM_Global")
name = "mNPTM_PCA.pdf"
ggplot(PCbeta,aes(x=PC1,y=PC2,col=group))+
  geom_point(size=3)+ # sets point size
  geom_text(aes(label=rownames(pData(bs.filtered.bsseq))),size=2,nudge_x=0,nudge_y=10)+ # sets label size and position
  theme_classic()+ # gets theme
  theme(axis.text=element_text(size=3), axis.title=element_text(size=8),plot.title = element_text(size=10))+ # sets size of axis labels and title
  scale_color_manual(values = brewer.pal(length(levels(as.factor(group))), "Set1"))+ # sets color of dots
  ggtitle(paste("PCA of Methylation at", ncol(beta) ,"CpG Sites",sep=" "))+ # adds title
  xlab(paste("PC 1 (", pov1, "%)", sep = ""))+ # labels PC1
  ylab(paste("PC 2 (", pov2, "%)", sep = "")) # labels PC2
ggsave(filename = paste("mNPTM_Global/", name), plot = last_plot(), device = pdf(), width = 5.5, height = 4, path = getwd()) # saves PCA plot
rm(beta,pca, PCbeta, group)

## Make global density plots ##
bs.filtered.bsseq %>% DMRichR::globalStats(genome = genome,
                                           testCovariate = testCovariate,
                                           adjustCovariate = adjustCovariate,
                                           matchCovariate = matchCovariate) %>%
  openxlsx::write.xlsx("mNPTM_Global/smoothed_globalStats.xlsx") # Smoothened statistics

windows <- bs.filtered.bsseq %>% DMRichR::windows(goi = goi)
CpGs <- bs.filtered.bsseq %>% DMRichR::CpGs()
CGi <- bs.filtered.bsseq %>% DMRichR::CGi(genome = genome)
plots <- c("windows", "CpGs", "CGi")
group =  bs.filtered.bsseq %>% pData() %>% dplyr::as_tibble() %>%
  dplyr::pull(!!"Stage") %>% forcats::fct_rev()

plot = function (matrix = matrix, group = NA) {print(glue::glue("Density plot of {length(matrix)} sites"))
  matrix %>% dplyr::as_tibble() %>% magrittr::set_colnames(paste(group, seq_along(1:length(group)))) %>% dplyr::transmute(
    Group1 = dplyr::select(., dplyr::contains(levels(group)[1])) %>% rowMeans() * 100, 
    Group2 = dplyr::select(., dplyr::contains(levels(group)[2])) %>% rowMeans() * 100,
    Group3 = dplyr::select(., dplyr::contains(levels(group)[3])) %>% rowMeans() * 100, 
    Group4 = dplyr::select(., dplyr::contains(levels(group)[4])) %>% rowMeans() * 100,) %>% 
    magrittr::set_colnames(c(levels(group)[1],levels(group)[2], levels(group)[3], levels(group)[4])) %>% 
    tidyr::gather(key = "variable", value = "value") %>% dplyr::mutate(variable = factor(.$variable, 
                                                                                         levels = levels(group))) %>% ggplot(aes(value, color = variable)) + 
    geom_density(size = 1.2) + labs(x = "Percent Methylation", y = "Density", color = "Group") + theme_classic() +
    scale_color_manual(values = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C")) + 
    scale_x_continuous(expand = c(0.05, 0.05), breaks = c(0, 25, 50, 75, 100)) + scale_y_continuous(expand = c(0, 0.001)) + 
    theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20), strip.text = element_text(size = 20), 
          legend.text = element_text(size = 20), legend.position = "bottom", legend.title = element_text(size = 20)) %>% return()
}

purrr::walk(plots,
            function(plotMatrix, group =  bs.filtered.bsseq %>% pData() %>% dplyr::as_tibble() %>%
                       dplyr::pull(!!"Stage") %>% forcats::fct_rev()){
              title <- dplyr::case_when(plotMatrix == "windows" ~ "20Kb Windows", plotMatrix == "CpGs" ~ "Single CpG", plotMatrix == "CGi" ~ "CpG Island")
              plotMatrix %>% get() %>% plot(group = group) %>% 
                ggplot2::ggsave(glue::glue("mNPTM_Global/{title} Density Plot.pdf"), plot = ., device = NULL, width = 11, height = 4)
            })

# Make PDO PCA and global density plots -----------------------------------------------------------
## Global variables ##
genome <- as.character("hg38") # Options: hg38, hg19, mm10, mm9, rheMac10, rheMac8, rn6, danRer11, galGal6, bosTau9, panTro6, dm6, susScr11, canFam3, TAIR10, or TAIR9)
coverage <- as.integer(1) # CpG coverage cutoff minimum value is 1
perGroup <- as.double(.75) # Percent of samples in all combinations of covariates meeting CpG coverage cutoff; Options: 0-1
minCpGs <- as.integer(5) # Minimum number of CpGs for a DMR
maxPerms <- as.integer(10) # Maximum number of permutations for the DMR analysis; no more than the # of samples
cutoff <- as.double(0.1) # Cutoff value [from 0 to 1] for the single CpG coefficient utilized to discover testable background regions
testCovariate <- as.character("CombinedStage") # Test covariate 
EnsDb <- FALSE

## Load and process samples ##
bs.filtered <- processBismark(files = list.files(path = getwd(), pattern = "*.txt.gz"),
                              meta = openxlsx::read.xlsx("sample_info.xlsx",colNames = TRUE) %>% dplyr::mutate_if(is.character, as.factor),
                              testCovariate = testCovariate, adjustCovariate = adjustCovariate, matchCovariate = matchCovariate,
                              coverage = coverage, cores = 1, perGroup = perGroup, sexCheck = sexCheck)
save(bs.filtered, file = "RData/bismark_NPDOs.RData")
#load("RData/bismark_NPDOs.RData")

bs.filtered.bsseq <- bsseq::BSmooth(bs.filtered, BPPARAM = BiocParallel::SnowParam(workers = 1, progressbar = TRUE)) # get individual smoothened methylation values
save(bs.filtered.bsseq,file = "RData/bsseq_NPDOs.RData")
#load("RData/bsseq_NPDOs.RData")

DMRichR::annotationDatabases(genome = genome, EnsDb = EnsDb) # setup annotation databases
settings_env <- ls(all = TRUE)
save(list = settings_env, file = "RData/settings_NPDOs.RData")
#load("RData/settings_NPDOs.RData")

## Perform PCA ##
group =  bs.filtered.bsseq %>% pData() %>% dplyr::as_tibble() %>%
  dplyr::pull(!!"Stage") %>% forcats::fct_rev()
beta = t(bs.filtered.bsseq %>% DMRichR::CpGs()) # get methylation values
beta = beta[,order(apply(beta,2,var))]
pca = prcomp(beta) # conducts principle component analysis
PCbeta = data.frame(pca$x, Group=group)
pov1 = summary(pca)$importance[2,1] # PC1 proportion of variance
pov1 = round(pov1*100,digits=2)
pov2 = summary(pca)$importance[2,2] # PC2 proportion of variance
pov2 = round(pov2*100,digits=2)

## Create PCA Plot ##
name = "PDO_labeled_PCA.pdf"
ggplot(PCbeta,aes(x=PC1,y=PC2,col=group))+
  geom_point(size=3)+ # sets point size
  geom_text(aes(label=rownames(pData(bs.filtered.bsseq))),size=2,nudge_x=0,nudge_y=5)+ # sets label size and position
  theme_classic()+ # gets theme
  theme(axis.text=element_text(size=3), axis.title=element_text(size=8),plot.title = element_text(size=10))+ # sets size of axis labels and title
  scale_color_manual(values = c("#4DAF4A", "#377EB8", "#E41A1C", "#FF7F00", "#984EA3")) + 
  ggtitle(paste("PCA of Methylation at", ncol(beta) ,"CpG Sites",sep=" "))+ # adds title
  xlab(paste("PC 1 (", pov1, "%)", sep = ""))+ # labels PC1
  ylab(paste("PC 2 (", pov2, "%)", sep = "")) # labels PC2
ggsave(filename = paste("mNPTM_Global/", name), plot = last_plot(), device = pdf(), width = 5.5, height = 4, path = getwd()) # saves PCA plot

## Make global density plots ##
bs.filtered.bsseq %>% DMRichR::globalStats(genome = genome,
                                           testCovariate = testCovariate,
                                           adjustCovariate = adjustCovariate,
                                           matchCovariate = matchCovariate) %>%
  openxlsx::write.xlsx("PDO_Global/smoothed_globalStats.xlsx") # Smoothened statistics

windows <- bs.filtered.bsseq %>% DMRichR::windows(goi = goi)
CpGs <- bs.filtered.bsseq %>% DMRichR::CpGs()
CGi <- bs.filtered.bsseq %>% DMRichR::CGi(genome = genome)  
plots <- c("windows", "CpGs", "CGi")
group =  bs.filtered.bsseq %>% pData() %>% dplyr::as_tibble() %>%
  dplyr::pull(!!"Stage") %>% forcats::fct_rev()

plot = function (matrix = matrix, group = NA) {print(glue::glue("Density plot of {length(matrix)} sites"))
  matrix %>% dplyr::as_tibble() %>% magrittr::set_colnames(paste(group, seq_along(1:length(group)))) %>% dplyr::transmute(
    Group1 = dplyr::select(., dplyr::contains(levels(group)[1])) %>% rowMeans() * 100, 
    Group2 = dplyr::select(., dplyr::contains(levels(group)[2])) %>% rowMeans() * 100,
    Group3 = dplyr::select(., dplyr::contains(levels(group)[3])) %>% rowMeans() * 100,
    Group4 = dplyr::select(., dplyr::contains(levels(group)[4])) %>% rowMeans() * 100,
    Group5 = dplyr::select(., dplyr::contains(levels(group)[5])) %>% rowMeans() * 100,) %>% 
    magrittr::set_colnames(c(levels(group)[1],levels(group)[2], levels(group)[3], levels(group)[4], levels(group)[5])) %>% 
    tidyr::gather(key = "variable", value = "value") %>% dplyr::mutate(variable = factor(.$variable, 
                                                                                         levels = levels(group))) %>% ggplot(aes(value, color = variable)) + 
    geom_density(size = 1.2) + labs(x = "Percent Methylation", y = "Density", color = "Group") + theme_classic() +
    scale_color_manual(values = c("#4DAF4A", "#377EB8", "#E41A1C", "#FF7F00", "#984EA3")) + 
    scale_x_continuous(expand = c(0.05, 0.05), breaks = c(0, 25, 50, 75, 100)) + scale_y_continuous(expand = c(0, 0.001)) + 
    theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20), strip.text = element_text(size = 20), 
          legend.text = element_text(size = 20), legend.position = "bottom", legend.title = element_text(size = 20)) %>% return()
}

dir.create("Global")
purrr::walk(plots,
            function(plotMatrix, group =  bs.filtered.bsseq %>% pData() %>% dplyr::as_tibble() %>%
                       dplyr::pull(!!"Stage") %>% forcats::fct_rev()){
              title <- dplyr::case_when(plotMatrix == "windows" ~ "20Kb Windows", plotMatrix == "CpGs" ~ "Single CpG", plotMatrix == "CGi" ~ "CpG Island")
              plotMatrix %>% get() %>% plot(group = group) %>% 
                ggplot2::ggsave(glue::glue("PDO_Global/{title} Density Plot.pdf"), plot = ., device = NULL, width = 11, height = 4)
              })
