# load packages
library(enrichR)
library(enrichplot)
library(dplyr)
library(ggplot2)

# Upload DMR files and make into list
load("RData/TM_DMRs.RData")
load("RData/MPN_DMRs.RData")
load("RData/EL_DMRs.RData")
load("RData/NT_DMRs.RData")
load("RData/sub_DMRs.RData")

# Human Ontologies --------
dir.create("PDO_Ontologies")
genome = "hg38"
DMRichR::annotationDatabases(genome = genome, EnsDb = TRUE)
# listEnrichrDbs()
dbs <- c("GO_Biological_Process_2023",
         "GO_Cellular_Component_2023",
         "GO_Molecular_Function_2023",
         "KEGG_2021_Human")

# hN vs hT
NT = NT_sigRegions %>%
  DMRichR::annotateRegions(TxDb = TxDb,annoDb = annoDb) %>%
  plyranges::filter(annotation != "Distal Intergenic") %>%
  dplyr::select(geneSymbol) %>%
  purrr::flatten() %>%
  enrichR::enrichr(dbs) %>%
  purrr::map(~ dplyr::arrange(., Adjusted.P.value))
NT %>%
  openxlsx::write.xlsx(file = glue::glue("PDO_Ontologies/NT_enrichr.xlsx"))

# Early vs Late
EL = EL_sigRegions %>%
  DMRichR::annotateRegions(TxDb = TxDb,annoDb = annoDb) %>%
  plyranges::filter(annotation != "Distal Intergenic") %>%
  dplyr::select(geneSymbol) %>%
  purrr::flatten() %>%
  enrichR::enrichr(dbs) %>%
  purrr::map(~ dplyr::arrange(., Adjusted.P.value))
EL %>%
  openxlsx::write.xlsx(file = glue::glue("PDO_Ontologies/EL_enrichr.xlsx"))

# Subtype
sub = sub_sigRegions %>%
  DMRichR::annotateRegions(TxDb = TxDb,annoDb = annoDb) %>%  
  plyranges::filter(annotation != "Distal Intergenic") %>%
  dplyr::select(geneSymbol) %>%
  purrr::flatten() %>%
  enrichR::enrichr(dbs) %>%
  purrr::map(~ dplyr::arrange(., Adjusted.P.value))
sub %>%
  openxlsx::write.xlsx(file = glue::glue("PDO_Ontologies/subtype_enrichr.xlsx"))

# Mouse Ontologies ------------------------------------
dir.create("mNPTM_Ontologies")
genome = "mm9"
DMRichR::annotationDatabases(genome = genome, EnsDb = TRUE)
dbs <- c("GO_Biological_Process_2023",
         "GO_Cellular_Component_2023",
         "GO_Molecular_Function_2023",
         "KEGG_2019_Mouse")

TM = TM_sigRegions %>%
  DMRichR::annotateRegions(TxDb = TxDb,annoDb = annoDb) %>%  
  plyranges::filter(annotation != "Distal Intergenic") %>%
  dplyr::select(geneSymbol) %>%
  purrr::flatten() %>%
  enrichR::enrichr(dbs) %>%
  purrr::map(~ dplyr::arrange(., Adjusted.P.value))
TM %>%
  openxlsx::write.xlsx(file = glue::glue("mNPTM_Ontologies/TM_enrichr.xlsx"))

MPN = MPN_sigRegions %>%
  DMRichR::annotateRegions(TxDb = TxDb,annoDb = annoDb) %>%  
  plyranges::filter(annotation != "Distal Intergenic") %>%
  dplyr::select(geneSymbol) %>%
  purrr::flatten() %>%
  enrichR::enrichr(dbs) %>%
  purrr::map(~ dplyr::arrange(., Adjusted.P.value))
MPN %>%
  openxlsx::write.xlsx(file = glue::glue("mNPTM_Ontologies/MPN_enrichr.xlsx"))

# GO dot plot for MPN and TM ------
top = 10
go = TM$KEGG_2019_Mouse[1:top,]
go = rbind(go, TM$GO_Biological_Process_2023[1:top,])
go = rbind(go, TM$GO_Cellular_Component_2023[1:top,])
go = rbind(go, TM$GO_Molecular_Function_2023[1:top,])
go = data.frame(go, Count = as.numeric(sapply(strsplit(go[,2], "/"), `[`, 1)))
go = data.frame(go, GeneRatio = sapply(strsplit(go[,2], "/"), function(x) as.numeric(x[1])/as.numeric(x[2])))
go = data.frame(go, Category = rep(c("1-KEGG", "2-Biological Process", "3-Cellular Component", "4-Molecular Function"),each=top))

p = ggplot(data = go, aes(x = reorder(Term, -Adjusted.P.value), y = Adjusted.P.value, color = GeneRatio, size = Count)) + 
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw() + 
  ylab("Adjusted p-value") + 
  xlab("") +
  coord_flip()
plot = p + facet_grid(rows = vars(Category), scales = "free_y")
ggsave(plot, width = 10, height = 11, file = "mNPTM_Ontologies/TM.pdf")

top = 10
go = MPN$KEGG_2019_Mouse[1:top,]
go = rbind(go, MPN$GO_Biological_Process_2023[1:top,])
go = rbind(go, MPN$GO_Cellular_Component_2023[1:top,])
go = rbind(go, MPN$GO_Molecular_Function_2023[1:top,])
go = data.frame(go, Count = as.numeric(sapply(strsplit(go[,2], "/"), `[`, 1)))
go = data.frame(go, GeneRatio = sapply(strsplit(go[,2], "/"), function(x) as.numeric(x[1])/as.numeric(x[2])))
go = data.frame(go, Category = rep(c("1-KEGG", "2-Biological Process", "3-Cellular Component", "4-Molecular Function"),each=top))

p = ggplot(data = go, aes(x = reorder(Term, -Adjusted.P.value), y = Adjusted.P.value, color = GeneRatio, size = Count)) + 
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw() + 
  ylab("Adjusted p-value") + 
  xlab("") +
  coord_flip()
plot = p + facet_grid(rows = vars(Category), scales = "free_y")
ggsave(plot, width = 10, height = 11, file = "mNPTM_Ontologies/MPN.pdf")

# GO dot plot for PDOs ------
top = 10
go = EL$KEGG_2021_Human[1:top,]
go = rbind(go, EL$GO_Biological_Process_2023[1:top,])
go = rbind(go, EL$GO_Cellular_Component_2023[1:top,])
go = rbind(go, EL$GO_Molecular_Function_2023[1:top,])
go = data.frame(go, Count = as.numeric(sapply(strsplit(go[,2], "/"), `[`, 1)))
go = data.frame(go, GeneRatio = sapply(strsplit(go[,2], "/"), function(x) as.numeric(x[1])/as.numeric(x[2])))
go = data.frame(go, Category = rep(c("1-KEGG", "2-Biological Process", "3-Cellular Component", "4-Molecular Function"),each=top))

p = ggplot(data = go, aes(x = reorder(Term, -Adjusted.P.value), y = Adjusted.P.value, color = GeneRatio, size = Count)) + 
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw() + 
  ylab("Adjusted p-value") + 
  xlab("") +
  coord_flip()
plot = p + facet_grid(rows = vars(Category), scales = "free_y")
ggsave(plot, width = 10, height = 11, file = "PDO_Ontologies/EL.pdf")

top = 10
go = NT$KEGG_2021_Human[1:top,]
go = rbind(go, NT$GO_Biological_Process_2023[1:top,])
go = rbind(go, NT$GO_Cellular_Component_2023[1:top,])
go = rbind(go, NT$GO_Molecular_Function_2023[1:top,])
go = data.frame(go, Count = as.numeric(sapply(strsplit(go[,2], "/"), `[`, 1)))
go = data.frame(go, GeneRatio = sapply(strsplit(go[,2], "/"), function(x) as.numeric(x[1])/as.numeric(x[2])))
go = data.frame(go, Category = rep(c("1-KEGG", "2-Biological Process", "3-Cellular Component", "4-Molecular Function"),each=top))

p = ggplot(data = go, aes(x = reorder(Term, -Adjusted.P.value), y = Adjusted.P.value, color = GeneRatio, size = Count)) + 
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw() + 
  ylab("Adjusted p-value") + 
  xlab("") +
  coord_flip()
plot = p + facet_grid(rows = vars(Category), scales = "free_y")
ggsave(plot, width = 10, height = 11, file = "PDO_Ontologies/NT.pdf")

top = 10
go = sub$KEGG_2021_Human[1:top,]
go = rbind(go, sub$GO_Biological_Process_2023[1:top,])
go = rbind(go, sub$GO_Cellular_Component_2023[1:top,])
go = rbind(go, sub$GO_Molecular_Function_2023[1:top,])
go = data.frame(go, Count = as.numeric(sapply(strsplit(go[,2], "/"), `[`, 1)))
go = data.frame(go, GeneRatio = sapply(strsplit(go[,2], "/"), function(x) as.numeric(x[1])/as.numeric(x[2])))
go = data.frame(go, Category = rep(c("1-KEGG", "2-Biological Process", "3-Cellular Component", "4-Molecular Function"),each=top))

p = ggplot(data = go, aes(x = reorder(Term, -Adjusted.P.value), y = Adjusted.P.value, color = GeneRatio, size = Count)) + 
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw() + 
  ylab("Adjusted p-value") + 
  xlab("") +
  coord_flip()
plot = p + facet_grid(rows = vars(Category), scales = "free_y")
ggsave(plot, width = 10, height = 11, file = "PDO_Ontologies/Subtype.pdf")

