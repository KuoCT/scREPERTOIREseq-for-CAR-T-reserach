# Set up environment ----
# Set working directory
working_dir <- file.path(this.path::here(),"..")
setwd(working_dir)

# Load libraries
# GEX
library(Seurat)
options(Seurat.object.assay.version = "v3")
# GSEA/ORA
library(WebGestaltR)
# Plots
library(ggplot2)
library(ggprism)
library(ComplexHeatmap)
library(scales)
library(patchwork)
# Tool
library(dplyr)
library(tibble)
# library(tidyr)
# library(stringr)
# library(purrr)
library(svglite)

# Create folders
exp_name <- tools::file_path_sans_ext(basename(this.path::this.path()))
exp_date <- strsplit(exp_name, "_")[[1]][2]

seu_obj_out <- file.path(getwd(), "seuratObjects")
GEX_out <- file.path(getwd(), "out", exp_name, "GEX")
GSEA_out <- file.path(getwd(), "out", exp_name, "GSEA")

dir.create(seu_obj_out, showWarnings = F, recursive = T)
dir.create(GEX_out, showWarnings = F, recursive = T)
dir.create(GSEA_out, showWarnings = F, recursive = T)

# Seed
set.seed(1234)

# Load functions
source("scripts/seuratFunctions.R")
source("scripts/toolFunctions.R")

# Import data ----
h5 <- Read10X_h5("cellrangerData/GEX_PDL1MoCARJNG/outs/filtered_feature_bc_matrix.h5")
print(h5[1:10, 1:5])

# Create Seurat Object ----
seu <- CreateSeuratObject(h5, assay = "RNA", project = "PDL1MoCARJNG")
rm(h5)

# Pre-processing ----
seu <- seurat_qc_flow(
  obj = seu, 
  out_path = GEX_out, 
  plot_limit_percent_mt = c(0, 35) 
)

# Normalize (SCTransform) ----
options(future.globals.maxSize = 4 * 1024^3) 
seu <- SCTransform(seu)

# PCA/UMAP ----
seu <- RunPCA(seu)
p <- ElbowPlot(seu, ndims = 50)
print(p)
saveSVG(p, file.path(GEX_out, "ElbowPlot.svg"), width = 6, height = 4.65)

seu <- RunUMAP(seu , dims = 1:50) %>%
  FindNeighbors(dims = 1:50) %>%
  FindClusters(resolution = 0.3)

improveUMAP(DimPlot(seu, pt.size = 1, shuffle = T) + labs(title = "UMAP: Clusters") + ggPreset("DimPlot"))

seu <- RenameSeuratClusters(
  seu, 
  paste0("C", 1:9), 
  new_column_name = "clusters_v1"
)

seu_sub <- subset(seu, idents = "C9", invert = T)

newPlotSet(TRUE)
p <- improveUMAP(DimPlot(seu_sub, pt.size = 1, group.by = "clusters_v1", shuffle = T) + labs(title = "UMAP: Clusters") + ggPreset("DimPlot"))
p %>% addPlot("UMAP")
saveSVG(p, file.path(GEX_out, "UMAP_Clusters.svg"), width = 4.5, height = 5.2, skip = T)

p <- FeaturePlot2(seu_sub, features = "EGFP") %>% improveUMAP()
p %>% addPlot("EGFP")
saveSVG(p, file.path(GEX_out, "FeaturePlot_EGFP.svg"), width = 4.5, height = 5.05, skip = T)

p <- FeaturePlot2(seu_sub, features = "CARTag") %>% improveUMAP()
p %>% addPlot("CARTag")
saveSVG(p, file.path(GEX_out, "FeaturePlot_CARTag.svg"), width = 4.5, height = 5.05, skip = T)

p <- FeaturePlot2(seu_sub, features = "CD4") %>% improveUMAP()
p %>% addPlot("CD4")
saveSVG(p, file.path(GEX_out, "FeaturePlot_CD4.svg"), width = 4.5, height = 5.05, skip = T)

p <- FeaturePlot2(seu_sub, features = "CD8A") %>% improveUMAP()
p %>% addPlot("CD8A")
saveSVG(p, file.path(GEX_out, "FeaturePlot_CD8A.svg"), width = 4.5, height = 5.05, skip = T)

p <- FeaturePlot2(seu_sub, features = "NFATC1") %>% improveUMAP()
p %>% addPlot("NFATC1")
saveSVG(p, file.path(GEX_out, "FeaturePlot_NFATC1.svg"), width = 4.5, height = 5.05, skip = T)

listPlotSet()
p <- wrap_plots(PlotSet(), ncol = 3)
savePNG(p, file.path(GEX_out, "FeaturePlots_Compile.png"),width = 8, height = 8, skip = T)
# closePlotSet()

# Module score ----
geneset_JAKSTAT <- read.csv("other/JAK_STAT_genes.csv") %>% pull(1)
seu_sub <- AddModuleScore2(seu_sub, geneset_JAKSTAT, ScoreName = "Gene set JAKSTAT")

newPlotSet(TRUE)
p <- improveUMAP(DimPlot(seu_sub, pt.size = 1, group.by = "clusters_v1", shuffle = T) + labs(title = "UMAP: Clusters") + ggPreset("DimPlot"))
p %>% addPlot("UMAP")

p <- FeaturePlot2(seu_sub, features = "Gene set JAKSTAT") %>% improveUMAP()
p %>% addPlot("Gene set JAKSTAT")
saveSVG(p, file.path(GEX_out, "ModuleScore_Gene_set_JAKSTAT.svg"), width = 4.5, height = 5.05, skip = T)

p <- FeaturePlot2(seu_sub, features = "EGFP") %>% improveUMAP()
p %>% addPlot("EGFP")

listPlotSet()
p <- wrap_plots(PlotSet(), ncol = 3)
savePNG(p, file.path(GEX_out, "ModuleScore_Compile.png"),width = 8, height = 4.2, skip = F)
# closePlotSet()

# Save RDS ----
saveRDS2(seu_sub, file.path(seu_obj_out, "seu_sub.rds"), overwrite = T)

# GSEA ----
DEG_all_clusters <- FindAllMarkers(seu_sub, logfc.threshold = 0) %>% select(gene, everything())
write.csv(DEG_all_clusters, file.path(GEX_out, "DEG_all_clusters.csv"), row.names = F)

for (i in unique(DEG_all_clusters$cluster)) {
  project_name <- paste0("Project_", i)
  project_dir <- file.path(GSEA_out, project_name)
  gsea_dir <- file.path(project_dir, paste0(project_name, "_GSEA"))
  score_path <- file.path(project_dir, paste0("interestingID_mappingTable_", i, ".txt"))
  
  # Run WebGestaltR
  if (dir.exists(project_dir)) {
    message("Skip WebGestaltR ", i, ": existing output found.")
  } else {
    WebGestaltR(
      enrichMethod = "GSEA",
      organism = "hsapiens",
      enrichDatabase = "pathway_KEGG",
      interestGene = DEG_all_clusters %>% filter(cluster == i) %>% select(gene, avg_log2FC),
      interestGeneType = "genesymbol",
      referenceSet = "genome_protein-coding",
      minNum = 5,
      maxNum = 2000,
      sigMethod = "top",
      fdrMethod = "BH",
      useWeightedSetCover = TRUE,
      usekMedoid = FALSE,
      fdrThr = 0.05,
      topThr = 25,
      setCoverNum = 10,
      kMedoid_k = 10,
      reportNum = 40,
      isOutput = TRUE,
      outputDirectory = GSEA_out,
      projectName = i
    )
  }
  
  # GSEA color bar
  svg_files <- list.files(path = gsea_dir, pattern = "\\.svg$", full.names = T, recursive = F)
  svg_files <- svg_files[!grepl("_color\\.svg$", basename(svg_files))]
  color_files <- file.path(gsea_dir, paste0(tools::file_path_sans_ext(basename(svg_files)), "_color.svg"))
  
  if (TRUE) svg_files <- svg_files[!file.exists(color_files)] # Set "FALSE" to redraw color bar
  if (length(svg_files) == 0) {
    message("Skip color bar ", i, ": all SVG files already processed.")
    next
  }
  
  for (svg in svg_files) {
    improveGSEAplot(
      svg_path = svg,
      score_path = score_path,
      output_folder = gsea_dir,
      hit_height_scale = 0.5,
      border_color = "black",
      border_width = 1.5
    )
  }
}

# PseudoBulk GEX ----
pseudoBulk <- AverageExpression(seu_sub)$RNA %>% as.matrix()
write.csv(pseudoBulk %>% mat2df("gene_symbol"), file.path(GEX_out, "PseudoBulk_all.csv"), row.names = F)

## Filtering
GRCh38_p14_gene_info <- read.csv("other/GRCh38_p14_gene_info.csv")
protein_coding_genes <- GRCh38_p14_gene_info %>%
  filter(Gene_type == "protein_coding") %>%
  select(Gene_stable_ID, Gene_name) %>%
  unlist(use.names = FALSE) %>%
  as.character() %>%
  unique()

filter_candidates <- c(protein_coding_genes, "CARTag")
pseudoBulk_filtered <- pseudoBulk %>%
  mat2df("gene_symbol") %>%
  mutate(
    gene_without_suffix = sub("\\.[0-9]+$", "", gene_symbol)
  ) %>%
  filter(
    gene_symbol %in% filter_candidates |
      gene_without_suffix %in% filter_candidates
  ) %>%
  select(-gene_without_suffix) %>%
  filter(rowSums(across(-gene_symbol) != 0) > 0) %>%
  column_to_rownames("gene_symbol") %>%
  df2mat()

nrow(pseudoBulk)
nrow(pseudoBulk_filtered)

pseudoBulk_filtered_Z <- pseudoBulk_filtered %>% as.matrix() %>% t() %>% scale() %>% t()
if (anyNA(pseudoBulk_filtered_Z)) pseudoBulk_filtered_Z[is.na(pseudoBulk_filtered_Z)] <- 0

write.csv(pseudoBulk_filtered, file.path(GEX_out, "pseudoBulk_filtered.csv"))
write.csv(pseudoBulk_filtered_Z,file.path(GEX_out, "pseudoBulk_filtered_Z.csv"))

capture.output(sessionInfo(), file = "sessionInfo_RNA.txt")
