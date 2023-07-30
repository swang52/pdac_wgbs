# Description
# CpG and genic enrichment analysis for mouse organoid DMRs
# requires DMR.RData objects from mNPTM_DMR.R as input

# Get Packages --------------------------------------------------------------
# Package names
packages <- c("DMRichR", "cowplot", "ggplot2")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Load packages
invisible(lapply(packages, library, character.only = TRUE))

# Upload DMR files and make into list ------
load("RData/TM_DMRs.RData")
load("RData/MPN_DMRs.RData")
load("RData/NT_DMRs.RData")
load("RData/EL_DMRs.RData")

dir.create("DMRichments")
# Enrichment multiplot functions ----------------------------------------
DMparseR <- function (direction = c("All DMRs", "Hypermethylated DMRs", "Hypomethylated DMRs"), 
                      type = c("CpG", "genic"),
                      project = c("TM", "MPN", "EL", "NT", "sub")) {
  stopifnot(direction %in% c("All DMRs", "Hypermethylated DMRs", 
                             "Hypomethylated DMRs"))
  stopifnot(type %in% c("CpG", "genic"))
  print(glue::glue("Parsing {type} enrichment results for {tidyDirection}", 
                   tidyDirection = glue::glue_collapse({
                     direction
                   }, sep = ", ", last = " and ")))
  purrr::map(direction, function(direction) {
    glue::glue("DMRichments/{project}_{direction}_{type}_enrichments.xlsx")
  }) %>% as.vector() %>% lapply(function(file) {
    readxl::read_xlsx(file)
  }) %>% `names<-`(direction) %>% data.table::rbindlist(idcol = "Dataset") %>% 
    dplyr::as_tibble() %>% tidyr::separate(Dataset, c("Direction", "DMR")) %>% 
    dplyr::select(Direction, Annotation, OR, fdr) %>% 
    dplyr::mutate(Annotation = forcats::as_factor(Annotation)) %>% 
    return()
}

countEnrich <- function (project = c("TM", "MPN", "EL", "NT", "sub"), type = c("CpG", "genic")) {
  data <- read.table(file = glue::glue("DMRichments/{project}_{type}_counts.txt"), header = T, sep = "\t")
    if(project == "TM"){title <- "mT vs mM"
      }else if(project == "MPN"){title <- "mN/mP vs mM"
        }else if(project == "EL"){title <- "Early vs Late"
          }else if(project == "NT"){title <- "mN vs mT"
            }else if(project == "sub"){title <-"Progenitor vs Squamous"
              }
    if(type == "CpG"){
      data$CpG <- factor(data$CpG, 
            levels=c("CpG.Island", "CpG.Shelf", "CpG.Shore", "Open.Sea"))
      ggplot(data, aes(fill=CpG, y=Percent, x=Direction)) + 
        geom_bar(position="dodge", stat = "identity", color = "black")+
        scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
        labs(x ="", y = "") + ggtitle(title) +
        scale_x_discrete(labels=c("Hypermethylated", "Hypomethylated")) +
        theme_minimal() + 
        scale_fill_manual(values = c("forestgreen", "goldenrod2", "dodgerblue", "blue3"))
    }else if(type == "genic"){
      data$Annotation <- factor(data$Annotation, 
            levels=c("Promoter", "5UTR", "Exon", "Intron", "3UTR", "Downstream", "Distal Intergenic"))
      ggplot(data, aes(fill=Annotation, y=Percent, x=Direction)) + 
        geom_bar(position="dodge", stat = "identity", color = "black")+
        scale_y_continuous(labels = scales::percent) +
        labs(x ="", y = "") + ggtitle(title) +
        scale_x_discrete(labels=c("Hypermethylated", "Hypomethylated")) +
        theme_minimal() + 
        scale_fill_manual(values=wesanderson::wes_palette("Zissou1", n = 7, type = "continuous") %>% rev())
    }
  }

# CpG and genic enrichment testing plots ----------------------------------------
theme = theme(text = element_text(size = 6), axis.text = element_text(size = 6), 
              axis.title = element_text(size = 6), strip.text = element_text(size = 6), 
              legend.text = element_text(size = 6)) 

theme2 = theme(axis.text.y = element_blank(),
               axis.ticks.y = element_blank(),
               axis.title.y = element_blank())
# CpG
MPN_CpG <- DMparseR(direction =  c("Hypermethylated DMRs","Hypomethylated DMRs"),type = "CpG", project = "MPN") %>%
                DMRichR::DMRichPlot(type = "CpG",multi = TRUE) + theme + ggtitle("mN/mP vs mM")
TM_CpG <- DMparseR(direction =  c("Hypermethylated DMRs","Hypomethylated DMRs"),type = "CpG", project = "TM") %>%
  DMRichR::DMRichPlot(type = "CpG",multi = TRUE) + theme + ggtitle("mT vs mM")
NT_CpG <- DMparseR(direction =  c("Hypermethylated DMRs","Hypomethylated DMRs"),type = "CpG", project = "NT") %>%
  DMRichR::DMRichPlot(type = "CpG",multi = TRUE) + theme + ggtitle("hN vs hT")
EL_CpG <- DMparseR(direction =  c("Hypermethylated DMRs","Hypomethylated DMRs"),type = "CpG", project = "EL") %>%
  DMRichR::DMRichPlot(type = "CpG",multi = TRUE) + theme + ggtitle("Early vs Late")
cpgPlot <- egg::ggarrange(MPN_CpG, TM_CpG + theme2, 
                            NT_CpG + theme2, EL_CpG + theme2, nrow=1)
ggplot2::ggsave(cpgPlot, width = 8, height = 2, filename = "DMRichments/CpG_multi_plot.pdf")

#Genic
MPN_genic <- DMparseR(direction =  c("Hypermethylated DMRs","Hypomethylated DMRs"),type = "genic", project = "MPN") %>%
  DMRichR::DMRichPlot(type = "genic",multi = TRUE) + ggtitle("mN/mP vs mM") + theme
TM_genic <- DMparseR(direction =  c("Hypermethylated DMRs","Hypomethylated DMRs"),type = "genic", project = "TM") %>%
  DMRichR::DMRichPlot(type = "genic",multi = TRUE) + ggtitle("mT vs mM") + theme
NT_genic <- DMparseR(direction =  c("Hypermethylated DMRs","Hypomethylated DMRs"),type = "genic", project = "NT") %>%
  DMRichR::DMRichPlot(type = "genic",multi = TRUE) + ggtitle("hN vs hT") + theme
EL_genic <- DMparseR(direction =  c("Hypermethylated DMRs","Hypomethylated DMRs"),type = "genic", project = "EL") %>%
  DMRichR::DMRichPlot(type = "genic",multi = TRUE) + ggtitle("Early vs Late") + theme
genicPlot <- egg::ggarrange(MPN_genic, TM_genic + theme2, 
                            NT_genic + theme2, EL_genic + theme2, nrow=1)
ggplot2::ggsave(genicPlot, width = 8, height = 3, filename = "DMRichments/genic_multi_plot.pdf")

genicCpG <- egg::ggarrange(MPN_CpG, TM_CpG + theme2, 
               NT_CpG + theme2, EL_CpG + theme2,
               MPN_genic, TM_genic + theme2, 
               NT_genic + theme2, EL_genic + theme2,
               nrow = 2, heights = c(1,2))
save(genicCpG, file = "RData/genicCpG_fig.RData")
ggplot2::ggsave(genicCpG, width = 8, height = 5, filename = "DMRichments/genic_CpG_multi_plot.pdf")

rm(MPN_CpG, TM_CpG, NT_CpG, EL_CpG, 
   MPN_genic, TM_genic, NT_genic, EL_genic)

# CpG and genic enrichment counts ----------------
MPN_CpG <- countEnrich(project = "MPN", type = "CpG")
TM_CpG <- countEnrich(project = "TM", type = "CpG")
NT_CpG <- countEnrich(project = "NT", type = "CpG")
EL_CpG <- countEnrich(project = "EL", type = "CpG")
ggsave(plot_grid(MPN_CpG + theme(legend.position = "none"), 
                 TM_CpG + theme(legend.position = "none"), 
                 NT_CpG + theme(legend.position = "none"), 
                 EL_CpG, nrow = 1), 
       width = 8, height = 2, file = "DMRichments/CpG_count_multiPlot.pdf")

MPN_genic <- countEnrich(project = "MPN", type = "genic")
TM_genic <- countEnrich(project = "TM", type = "genic")
NT_genic <- countEnrich(project = "NT", type = "genic")
EL_genic <- countEnrich(project = "EL", type = "genic")

ggsave(plot_grid(MPN_genic + theme(legend.position = "none"), 
                 TM_genic + theme(legend.position = "none"),
                 NT_genic + theme(legend.position = "none"), 
                 EL_genic, nrow = 1), 
       width = 8, height = 2, file = "DMRichments/genic_count_multiPlot.pdf")

ggsave(plot_grid(MPN_CpG + theme(legend.position = "none"), 
                 TM_CpG + theme(legend.position = "none"), 
                 NT_CpG + theme(legend.position = "none"), 
                 EL_CpG, 
                 MPN_genic + theme(legend.position = "none"), 
                 TM_genic + theme(legend.position = "none"),
                 NT_genic + theme(legend.position = "none"), 
                 EL_genic, nrow = 2, align = "v"), 
       width = 8, height = 6, file = "DMRichments/genic_CpG_multi.pdf")

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