# Description -------------
# Isolates ChromHMM enrichment from pancreas for all DMR sets and plots heatmap

# Get Packages --------------------------------------------------------------
# Package names
packages <- c("readr", "ComplexHeatmap", "circlize", "dplyr")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Load packages
invisible(lapply(packages, library, character.only = TRUE))

# Extract chromHMM enrichment from pancreas --------
chromTable <- function(project = project, celltype = celltype){
  dataset <- vector(mode='list', length=length(project))
  names(dataset) <- project
  for (i in 1:length(project)){
    dataset[[i]] <- readr::read_tsv(glue::glue("LOLA_{project[i]}/Hypermethylated DMRs/ChromHMM/allEnrichments.tsv")) %>%
    as.data.frame() %>% dplyr::filter(cellType == celltype) %>% arrange(description) %>%
    rbind(readr::read_tsv(glue::glue("LOLA_{project[i]}/Hypomethylated DMRs/ChromHMM/allEnrichments.tsv")) %>%
    as.data.frame() %>% dplyr::filter(cellType == celltype) %>% arrange(description))
  }
  dataset
}

project <- c("MPN", "TM", "EL", "NT", "Subtype")
celltype <- "Pancreas"
data <- chromTable(project = project, celltype = celltype)
openxlsx::write.xlsx(data, file = 'ChromHMM_pancreas.xlsx') 

# Read ChromHMM enrichment from pancreas ---------
#sheets <- readxl::excel_sheets("ChromHMM_pancreas.xlsx")
#data <- lapply(sheets, function(x) readxl::read_excel("ChromHMM_pancreas.xlsx", sheet = x))
#data <- lapply(data, as.data.frame)
#names(data) <- sheets

# Make ChromHMM fold enrichment heatmaps -------------
states <- c("Active TSS","Flanking Active TSS", "Transcription at Gene 5' and 3'", "Strong Transcription", 
            "Weak Transcription", "Genic Enhancers", "Enhancers", "ZNF Genes & Repeats", "Heterochromatin",
            "Bivalent/Poised TSS", "Flanking Bivalent TSS/Enhancer", "Bivalent Enhancer", 
            "Repressed PolyComb", "Weak Repressed PolyComb", "Quiescent/Low")

# odds ratio matrix
data <- lapply(data, function(x) split(x, f = x$userSet))
odds <- lapply(data, function(x) lapply(x, function (x) x$oddsRatio))
odds <- t(as.data.frame(odds))
colnames(odds) <- states

# fold enrichment matrix
index <- which(odds < 1)
fold <- odds
for (x in index){
  fold[x] = -1/odds[x]
}
fold[fold == -Inf] = 0
colnames(fold) <- states

# q value matrix
q <- lapply(data, function(x) lapply(x, function (x) as.numeric(x$qValue)))
q <- t(as.data.frame(q))
colnames(q) <- states

# Stage heatmap
stage_fold <- fold[1:8,]
stage_q <- q[1:8,]
pdf(file = "stage_chromHMM.pdf", width = 11, height = 8)
stage_heat <- Heatmap(stage_fold, cell_fun = function(j, i, x, y, width, height, fill) {
  if(stage_q[i, j] < 0.001) {
    grid.text("***", x, y)
  } else if(stage_q[i, j] < 0.01) {
    grid.text("**", x, y)
  } else if(stage_q[i, j] < 0.05) {
    grid.text("*", x, y)
  }
}, row_order = c(1:8), column_order = c(1:15),
row_split = rep(c("A","B","C","D"), each = 2),
col=colorRamp2(c(-5, 0, 15), c("blue", "white", "red")),
rect_gp = gpar(col = "black", lwd = .5),
width = ncol(stage_fold)*unit(10, "mm"), height = nrow(stage_fold)*unit(10, "mm"),
name = "Fold Enrichment")
dev.off()

# Subtype heatmap ------------
sub_fold <- fold[9:10,]
sub_q <- q[9:10,]
rownames(sub_fold) <- c("Hyper", "Hypo")
rownames(sub_q) <- c("Hyper", "Hypo")

pdf(file = "subtype_chromHMM.pdf", width = 11, height = 8)
Heatmap(sub_fold, cell_fun = function(j, i, x, y, width, height, fill){
  if(sub_q[i, j] < 0.001) {
    grid.text("***", x, y)
  } else if(sub_q[i, j] < 0.01) {
    grid.text("**", x, y)
  } else if(sub_q[i, j] < 0.05) {
    grid.text("*", x, y)
  }
}, row_order = c(1:2), column_order = c(1:15), 
  width = ncol(sub_fold)*unit(10, "mm"), height = nrow(sub_fold)*unit(10, "mm"),
  rect_gp = gpar(col = "black", lwd = .5),
  col=colorRamp2(c(min(sub_fold), 0, max(sub_fold)), c("blue", "white", "red")),
  name = "Fold Enrichment")
dev.off()
