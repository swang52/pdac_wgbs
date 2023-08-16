# Description
# Performs makes HOMER enrichment graphs for all DMR sets

# Get Packages --------------------------------------------------------------
# Package names
packages <- c("ggplot2", "dplyr", "enrichplot", "cowplot")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Load packages
invisible(lapply(packages, library, character.only = TRUE))

# HOMER graph function -----------------
HOMER <- function(type = c("TM", "MPN", "EL", "NT", "sub"), top = top){
  read.table(file = glue::glue("HOMER/{type}_hyper_output/knownResults.txt"), row.names = NULL) %>% 
    dplyr::slice(1:top) %>%
    rbind(read.table(file = glue::glue("HOMER/{type}_hypo_output/knownResults.txt"), row.names = NULL) %>%
            dplyr::slice(1:top)) %>%
    stats::setNames(c("Motif", "Consensus","p.value", "logP", "q.value","countTarget","percentTarget","countBG","percentBG")) %>%
    dplyr::mutate(Motif = sapply(strsplit(Motif, "/"), function(x) x[1])) %>%
    dplyr::mutate(Category = rep(c("Hypermethylated", "Hypomethylated"),each=top)) %>%
    ggplot2::ggplot(ggplot2::aes(x = reorder(Motif, -logP/2.303), y = -logP/2.303)) +  
    ggplot2::geom_bar(stat="identity") +
    ggplot2::labs(y = "-log(P value)", x = "") +
    ggplot2::coord_flip() +
    ggplot2::theme_bw() + 
    ggplot2::theme(text = element_text(size = 5),
                   panel.grid.minor = element_line(size = 0.2), panel.grid.major = element_line(size = .2),
                   axis.ticks=element_line(size=.2)) +  ggplot2::ggtitle(type) + 
    ggplot2::facet_grid(rows = vars(Category), scales = "free_y")
}

# Make plots --------------------
top <- 10
TMplot <- HOMER(type = "TM", top)
MPNplot <- HOMER(type = "MPN", top)
ggsave(plot_grid(MPNplot, TMplot, ncol=1, align="v"), width = 1.7, height = 8.5, file = "HOMER/NPTM_HOMER.pdf")

ELplot <- HOMER(type = "EL", top)
NTplot <- HOMER(type = "NT", top)
ggsave(plot_grid(ELplot, NTplot, ncol=1, align="v"), width = 1.7, height = 8.5, file = "HOMER/NTEL_HOMER.pdf")

Subplot <- HOMER(type = "sub", top)
ggsave(Subplot, width = 5, height = 4, file = "HOMER/Sub.pdf")