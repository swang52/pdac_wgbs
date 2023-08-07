# Description
# Performs gene ontology enrichment analysis for all DMR sets

# Get Packages --------------------------------------------------------------
# Package names
packages <- c("enrichR", "enrichplot", "dplyr", "ggplot2", "cowplot")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Load packages
invisible(lapply(packages, library, character.only = TRUE))

# Upload DMR files and make into list -------
load("RData/TM_DMRs.RData")
load("RData/MPN_DMRs.RData")
load("RData/EL_DMRs.RData")
load("RData/NT_DMRs.RData")
load("RData/sub_DMRs.RData")

# GO functions -------
GO <- function(sigRegions = sigRegions, dbs = dbs, TxDb, annoDb){
  sigRegions %>%
    DMRichR::annotateRegions(TxDb = TxDb,annoDb = annoDb) %>%
    plyranges::filter(annotation != "Distal Intergenic") %>%
    dplyr::select(geneSymbol) %>%
    purrr::flatten() %>%
    enrichR::enrichr(dbs) %>%
    purrr::map(~ dplyr::arrange(., Adjusted.P.value))
}

GOplot <- function(data = data, dbs = dbs, top = top){
  go <- data[[1]][1:top,]
  for (i in 2:length(dbs)){
    go <- rbind(go, data[[i]][1:top,])
  }
  go <- go %>%
    dplyr::mutate(Term = gsub("\\s*\\([^\\)]+\\)","",as.character(go[,1]))) %>%
    data.frame(., Count = as.numeric(sapply(strsplit(go[,2], "/"), `[`, 1)),
               GeneRatio = sapply(strsplit(go[,2], "/"), function(x) as.numeric(x[1])/as.numeric(x[2])),
               Category = rep(dbs, each=top))
  ggplot(data = go, aes(x = reorder(Term, -log10(P.value)), y = -log10(P.value), color = -log10(Adjusted.P.value), size = Count)) + 
    geom_point() +
    scale_color_gradient(low = "blue", high = "red") +
    theme_bw() + 
    theme(text = element_text(size = 5),
          panel.grid.minor = element_line(size = 0.2), panel.grid.major = element_line(size = .2),
          axis.ticks=element_line(size=.2)) +
    ylab("-log10(P value)") + 
    xlab("") +
    coord_flip() + 
    facet_grid(rows = vars(Category), scales = "free_y")
}

# Human Ontologies --------
dir.create("PDO_Ontologies")
DMRichR::annotationDatabases(genome = "hg38", EnsDb = TRUE)
dbs <- c("GO_Biological_Process_2023",
         "GO_Cellular_Component_2023",
         "GO_Molecular_Function_2023",
         "KEGG_2021_Human")

NT <- GO(sigRegions = NT_sigRegions, dbs = dbs, TxDb = TxDb, annoDb = annoDb)
NT %>% openxlsx::write.xlsx(file = "PDO_Ontologies/NT_enrichr.xlsx")

EL <- GO(sigRegions = EL_sigRegions, dbs = dbs, TxDb = TxDb, annoDb = annoDb)
EL %>% openxlsx::write.xlsx(file = "PDO_Ontologies/EL_enrichr.xlsx")

sub <- GO(sigRegions = sub_sigRegions, dbs = dbs, TxDb = TxDb, annoDb = annoDb)
sub %>% openxlsx::write.xlsx(file = "PDO_Ontologies/sub_enrichr.xlsx")

# GO dot plot for PDOs ------
top <- 8
pNT <- GOplot(data = NT, dbs = dbs, top = top) +
  scale_size(breaks = c(25, 50, 75), range = c(.5, 2.5)) 
ggsave(pNT, width = 5, height = 5, file = "PDO_Ontologies/NT.pdf")

pEL <- GOplot(data = EL, dbs = dbs, top = top) +
  scale_size(breaks = c(10, 20, 30), range = c(.5, 2.5)) 
ggsave(pEL, width = 5, height = 5, file = "PDO_Ontologies/EL.pdf")

pSub <- GOplot(data = sub, dbs = dbs, top = top)
ggsave(pSub, width = 5, height = 5, file = "PDO_Ontologies/sub.pdf")

# Mouse Ontologies ------------------------------------
dir.create("mNPTM_Ontologies")
DMRichR::annotationDatabases(genome = "mm9", EnsDb = TRUE)
dbs <- c("GO_Biological_Process_2023",
         "GO_Cellular_Component_2023",
         "GO_Molecular_Function_2023",
         "KEGG_2019_Mouse")

TM <- GO(sigRegions = TM_sigRegions, dbs = dbs, TxDb = TxDb, annoDb = annoDb)
TM %>% openxlsx::write.xlsx(file = "mNPTM_Ontologies/TM_enrichr.xlsx")
MPN <- GO(sigRegions = MPN_sigRegions, dbs = dbs, TxDb = TxDb, annoDb = annoDb)
MPN %>% openxlsx::write.xlsx(file = "mNPTM_Ontologies/MPN_enrichr.xlsx")

# GO dot plot for MPN and TM ------
top <- 8
pMPN <- GOplot(data = MPN, dbs = dbs, top = top) +
  scale_size(breaks = c(50, 100, 150), range = c(.5, 2.5)) +
  scale_y_continuous(breaks = seq(0, 25, by = 5), limits = c(0,23)) +
  scale_color_gradient(low = "blue", high = "red", limits = c(3,21))
ggsave(pMPN, width = 5, height = 5, file = "mNPTM_Ontologies/MPN.pdf")

pTM <- GOplot(data = TM, dbs = dbs, top = top) +
  scale_size(breaks = c(50, 100, 150), range = c(.5, 2.5)) +
  scale_y_continuous(breaks = seq(0, 25, by = 5), limits = c(0,23)) +
  scale_color_gradient(low = "blue", high = "red", limits = c(3,21))
ggsave(pTM, width = 5, height = 5, file = "mNPTM_Ontologies/TM.pdf")

ggsave(plot_grid(pMPN, pTM, ncol=1, align="v"), width = 5, height = 7.5, file = "mNPTM_Ontologies/NPTM.pdf")