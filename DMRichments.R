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
load("RData/sub_DMRs.RData")

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
        theme(text = element_text(size = 6), axis.text = element_text(size = 6), 
              axis.title = element_text(size = 6), legend.text = element_text(size = 6)) +
        scale_fill_manual(values = c("forestgreen", "goldenrod2", "dodgerblue", "blue3"))
    }else if(type == "genic"){
      data$Annotation <- factor(data$Annotation, 
            levels=c("Promoter", "5UTR", "Exon", "Intron", "3UTR", "Downstream", "Distal Intergenic"))
      ggplot(data, aes(fill=Annotation, y=Percent, x=Direction)) + 
        geom_bar(position="dodge", stat = "identity", color = "black")+
        scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
        labs(x ="", y = "") + ggtitle(title) +
        scale_x_discrete(labels=c("Hypermethylated", "Hypomethylated")) +
        theme_minimal() + 
        theme(text = element_text(size = 6), axis.text = element_text(size = 6), 
              axis.title = element_text(size = 6), legend.text = element_text(size = 6)) +
        scale_fill_manual(values=wesanderson::wes_palette("Zissou1", n = 7, type = "continuous") %>% rev())
    }
  }

# CpG and genic enrichment testing plots ----------------------------------------
theme = theme(text = element_text(size = 6), axis.text = element_text(size = 6), 
              axis.title = element_text(size = 6), strip.text = element_text(size = 6), 
              legend.text = element_text(size = 6)) + geom_bar(stat = "identity", color = "Black", size = .5)

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
sub_CpG <- DMparseR(direction =  c("Hypermethylated DMRs","Hypomethylated DMRs"),type = "CpG", project = "sub") %>%
  DMRichR::DMRichPlot(type = "CpG",multi = TRUE) + theme + ggtitle("Progenitor vs Squamous")

#Genic
MPN_genic <- DMparseR(direction =  c("Hypermethylated DMRs","Hypomethylated DMRs"),type = "genic", project = "MPN") %>%
  DMRichR::DMRichPlot(type = "genic",multi = TRUE) + ggtitle("mN/mP vs mM") + theme
TM_genic <- DMparseR(direction =  c("Hypermethylated DMRs","Hypomethylated DMRs"),type = "genic", project = "TM") %>%
  DMRichR::DMRichPlot(type = "genic",multi = TRUE) + ggtitle("mT vs mM") + theme
NT_genic <- DMparseR(direction =  c("Hypermethylated DMRs","Hypomethylated DMRs"),type = "genic", project = "NT") %>%
  DMRichR::DMRichPlot(type = "genic",multi = TRUE) + ggtitle("hN vs hT") + theme
EL_genic <- DMparseR(direction =  c("Hypermethylated DMRs","Hypomethylated DMRs"),type = "genic", project = "EL") %>%
  DMRichR::DMRichPlot(type = "genic",multi = TRUE) + ggtitle("Early vs Late") + theme
sub_genic <- DMparseR(direction =  c("Hypermethylated DMRs","Hypomethylated DMRs"),type = "genic", project = "sub") %>%
  DMRichR::DMRichPlot(type = "genic",multi = TRUE) + theme + ggtitle("Progenitor vs Squamous")

genicCpG <- egg::ggarrange(MPN_CpG, TM_CpG + theme2, 
               NT_CpG + theme2, EL_CpG + theme2,
               MPN_genic, TM_genic + theme2, 
               NT_genic + theme2, EL_genic + theme2,
               nrow = 2, heights = c(1,2))
ggplot2::ggsave(genicCpG, width = 8, height = 5, filename = "DMRichments/genic_CpG_multi_plot.pdf")

SubgenicCpG <- egg::ggarrange(sub_CpG, sub_genic, nrow=2)
ggplot2::ggsave(SubgenicCpG, width = 8, height = 5, file = "DMRichments/sub_genic_CpG_multiPlot.pdf")

rm(MPN_CpG, TM_CpG, NT_CpG, EL_CpG, 
   MPN_genic, TM_genic, NT_genic, EL_genic,
   sub_CpG, sub_genic)

# CpG and genic enrichment counts ----------------
theme = theme(legend.position = "none")
theme2 = theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.title.y = element_blank())
MPN_CpG <- countEnrich(project = "MPN", type = "CpG")
TM_CpG <- countEnrich(project = "TM", type = "CpG")
NT_CpG <- countEnrich(project = "NT", type = "CpG")
EL_CpG <- countEnrich(project = "EL", type = "CpG")

MPN_genic <- countEnrich(project = "MPN", type = "genic")
TM_genic <- countEnrich(project = "TM", type = "genic")
NT_genic <- countEnrich(project = "NT", type = "genic")
EL_genic <- countEnrich(project = "EL", type = "genic")

genicCpG <- egg::ggarrange(MPN_CpG + theme, 
                           TM_CpG + theme + theme2, 
                           NT_CpG + theme + theme2, 
                           EL_CpG + theme2, 
                           MPN_genic + theme, 
                           TM_genic + theme + theme2, 
                           NT_genic + theme + theme2, 
                           EL_genic + theme2, nrow=2)
ggplot2::ggsave(genicCpG, width = 8, height = 5, file = "DMRichments/genic_CpG_count_multiPlot.pdf")

sub_CpG <- countEnrich(project = "sub", type = "CpG")
sub_genic <- countEnrich(project = "sub", type = "genic")
SubgenicCpG <- egg::ggarrange(sub_CpG, sub_genic, nrow=1)
ggplot2::ggsave(SubgenicCpG, width = 8, height = 5, file = "DMRichments/sub_genic_CpG_count_multiPlot.pdf")

