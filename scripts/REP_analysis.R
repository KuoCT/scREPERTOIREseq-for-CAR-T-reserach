# Set up environment ----
# Set working directory
working_dir <- file.path(this.path::here(),"..")
setwd(working_dir)

# Load libraries
# GEX
library(Seurat)
options(Seurat.object.assay.version = "v3")
# Plots
library(ComplexHeatmap)
library(ggplot2)
library(ggprism)
library(scales)
library(patchwork)
# Tool
library(dplyr)
library(tibble)
library(tidyr)
library(readr)
library(stringr)
library(purrr)
library(svglite)

# Create folders
exp_name <- tools::file_path_sans_ext(basename(this.path::this.path()))
exp_date <- strsplit(exp_name, "_")[[1]][2]

seu_obj_out <- file.path(getwd(), "seuratObjects")
REP_out <- file.path(getwd(), "out", exp_name, "REP")
out <- file.path(getwd(), "out", exp_name)
trainingData_out <- file.path(getwd(), "trainingData")

dir.create(seu_obj_out, showWarnings = F, recursive = T)
dir.create(REP_out, showWarnings = F, recursive = T)
dir.create(out, showWarnings = F, recursive = T)
dir.create(trainingData_out, showWarnings = F, recursive = T)

# Seed
set.seed(1234)

# Load functions
source("scripts/seuratFunctions.R")
source("scripts/toolFunctions.R")
source("scripts/scRepertoireTool.R")

# Import data ----
if (!exists("contig_1")) contig_1 <- readContig("cellrangerData/scFv_PDL1MoCARJNG_1/", sample_id = "scFv_PDL1MoCARJNG")
if (!exists("contig_2")) contig_2 <- readContig("cellrangerData/scFv_PDL1MoCARJNG_2/", sample_id = "scFv_PDL1MoCARJNG")

length(contig_1$contig_id)
length(contig_2$contig_id)

# write.csv(contig_1, file.path(REP_out, "contig_1.csv"), row.names = F)
# write.csv(contig_2, file.path(REP_out, "contig_2.csv"), row.names = F)

# UMI QC (filter 1) ----
contig_1_filter1 <- contig_1 %>% filter(umis >= 3)
contig_2_filter1 <- contig_2 %>% filter(umis >= 3)

write.csv(contig_1_filter1, file.path(REP_out, "contig_1_filter1_UMI.csv"), row.names = F)
write.csv(contig_2_filter1, file.path(REP_out, "contig_2_filter1_UMI.csv"), row.names = F)

# Franking site QC (filter 2) ----
franking_table <- data.frame(
  query_id = c("Myc_tag", "3_tail"),
  sequence = c(
    "GAGCAGAAACTTATTTCAGAGGAAGACCTG",
    "GGACCGACCACGACGCCAGCG"
  )
)

contig_1_frankingSite <- findFrankingSite(
  df = contig_1_filter1,
  franking_table = franking_table,
  sequence_col = "sequence",
  id_col = "contig_id",
  blast_program = "blastn",
  min_score = 90, 
  min_pident = 90,
  keep_reverse_strand = F,
  output_dir = file.path(REP_out, "BLAST_contig_1_frankingSite")
)

contig_2_frankingSite <- findFrankingSite(
  df = contig_2_filter1,
  franking_table = franking_table,
  sequence_col = "sequence",
  id_col = "contig_id",
  blast_program = "blastn",
  min_score = 90, 
  min_pident = 90,
  keep_reverse_strand = F,
  output_dir = file.path(REP_out, "BLAST_contig_2_frankingSite")
)

# contig_1_franking_site_hits <- attr(contig_1_frankingSite, "blast_hits")
# contig_2_franking_site_hits <- attr(contig_2_frankingSite, "blast_hits")

contig_1_filter2 <- cutFrankingSite(contig_1_frankingSite)
contig_2_filter2 <- cutFrankingSite(contig_2_frankingSite)

write.csv(contig_1_filter2, file.path(REP_out, "contig_1_filter2_frankingSite.csv"), row.names = F)
write.csv(contig_2_filter2, file.path(REP_out, "contig_2_filter2_frankingSite.csv"), row.names = F)

# scFv length QC (filter 3) ----
newPlotSet(TRUE)
p <- seqLengthPlot(contig_1_filter2, colors = c("#F8766D", "black", "#F8766D", "#7B1400")) + 
  coord_cartesian(xlim = c(0, 1000)) + 
  geom_vline(xintercept = c(600, 900), linetype = "dashed") +
  labs(title = "scFv Lengths", subtitle = "Batch 1") + 
  ggPreset("default")
p %>% addPlot("Batch 1")
saveSVG(p, file.path(out, "contig_1_lengths.svg"), width = 5, height = 4)

p <- seqLengthPlot(contig_1_filter2 %>% filter(insert_seq_length > 100), colors = c("#F8766D", "black", "#F8766D", "#7B1400")) + 
  coord_cartesian(xlim = c(0, 1000)) + 
  geom_vline(xintercept = c(600, 900), linetype = "dashed") +
  labs(title = "scFv Lengths", subtitle = "Batch 1 narrowed") + 
  ggPreset("default")
p %>% addPlot("Batch 1 narrowed")
saveSVG(p, file.path(out, "contig_1_lengths_v2.svg"), width = 5, height = 4)

p <- seqLengthPlot(contig_2_filter2, colors = c("#00BFC4", "black", "#00BFC4", "#044A4D")) + 
  coord_cartesian(xlim = c(0, 1000)) + 
  geom_vline(xintercept = c(600, 900), linetype = "dashed") +
  labs(title = "scFv Lengths", subtitle = "Batch 2") + 
  ggPreset("default")
p %>% addPlot("Batch 2")
saveSVG(p, file.path(out, "contig_2_lengths.svg"), width = 5, height = 4)

p <- seqLengthPlot(contig_2_filter2 %>% filter(insert_seq_length > 100), colors = c("#00BFC4", "black", "#00BFC4", "#044A4D")) + 
  coord_cartesian(xlim = c(0, 1000)) + 
  geom_vline(xintercept = c(600, 900), linetype = "dashed") +
  labs(title = "scFv Lengths", subtitle = "Batch 2 narrowed") + 
  ggPreset("default")
p %>% addPlot("Batch 2 narrowed")
saveSVG(p, file.path(out, "contig_2_lengths_v2.svg"), width = 5, height = 4)

listPlotSet()
p <- wrap_plots(PlotSet(), ncol = 2)
savePNG(p, file.path(out, "contig_lengthPlots.png"), width = 8, height = 9)

cols_to_remove <- c("Myc_tag_count", "Myc_tag_position", "Myc_tag_score", 
                    "3_tail_count", "3_tail_position", "3_tail_score", 
                    "left_hit_start", "left_hit_end", "left_hit_score", 
                    "right_hit_start", "right_hit_end", "right_hit_score", "sequence_length", "valid_franking_order")

contig_1_filter3 <- contig_1_filter2 %>% filter(insert_seq_length <= 900, insert_seq_length >= 600) %>% select(-all_of(cols_to_remove))
contig_2_filter3 <- contig_2_filter2 %>% filter(insert_seq_length <= 900, insert_seq_length >= 600) %>% select(-all_of(cols_to_remove))

write.csv(contig_1_filter3, file.path(REP_out, "contig_1_filter3_scFvLength.csv"), row.names = F)
write.csv(contig_2_filter3, file.path(REP_out, "contig_2_filter3_scFvLength.csv"), row.names = F)

# Find linker QC (filter 4)----
linkers <- read.csv("other/linker_list.csv")

contig_1_linker <- findLinker(
  df = contig_1_filter3,
  linker_table = linkers,
  sequence_col = "insert_seq",
  id_col = "contig_id",
  min_score = 90,
  min_pident = 90,
  keep_reverse_strand = FALSE,
  blast_program = "blastn",
  task = "blastn",
  dust = "no",
  output_dir = file.path(REP_out, "BLAST_contig_1_linker")
)

contig_2_linker <- findLinker(
  df = contig_2_filter3,
  linker_table = linkers,
  sequence_col = "insert_seq",
  id_col = "contig_id",
  min_score = 90,
  min_pident = 90,
  keep_reverse_strand = FALSE,
  blast_program = "blastn",
  task = "blastn",
  dust = "no",
  output_dir = file.path(REP_out, "BLAST_contig_2_linker")
)

# contig_1_linker_hits <- attr(contig_1_linker, "blast_hits")
# contig_2_linker_hits <- attr(contig_2_linker, "blast_hits")

contig_1_filter4 <- contig_1_linker %>% dplyr::filter(!is.na(linker_name)) %>% cutLinker()
contig_2_filter4 <- contig_2_linker %>% dplyr::filter(!is.na(linker_name)) %>% cutLinker()

write.csv(contig_1_filter4, file.path(REP_out, "contig_1_filter4_linker.csv"), row.names = F)
write.csv(contig_2_filter4, file.path(REP_out, "contig_2_filter4_linker.csv"), row.names = F)

## (No linker) ----
contig_1_filter4_nolinker <- contig_1_linker %>% dplyr::filter(is.na(linker_name))
contig_2_filter4_nolinker <- contig_2_linker %>% dplyr::filter(is.na(linker_name))

write.csv(contig_1_filter4_nolinker, file.path(REP_out, "contig_1_filter4_nolinker.csv"), row.names = F)
write.csv(contig_2_filter4_nolinker, file.path(REP_out, "contig_2_filter4_nolinker.csv"), row.names = F)

# Find known sequences ----
known_sequences <- read.csv("other/known_CAR_sequence.csv") %>%
  dplyr::mutate(sequence = paste0(V_left, linker, V_right))

cols_to_remove <- c("linker_position", "linker_score", "linker_start", "linker_end")

scFv_1_annotation <- findKnownSequence(
  df = contig_1_filter4 %>% select(-all_of(cols_to_remove)),
  known_table = known_sequences,
  sequence_col = "insert_seq",
  id_col = "contig_id",
  min_score = 90,
  min_pident = 90,
  keep_reverse_strand = F,
  blast_program = "blastn",
  task = "blastn",
  output_dir = file.path(REP_out, "BLAST_scFv_1_known_sequence")
)

scFv_2_annotation <- findKnownSequence(
  df = contig_2_filter4 %>% select(-all_of(cols_to_remove)),
  known_table = known_sequences,
  sequence_col = "insert_seq",
  id_col = "contig_id",
  min_score = 90,
  min_pident = 90,
  keep_reverse_strand = F,
  blast_program = "blastn",
  task = "blastn",
  output_dir = file.path(REP_out, "BLAST_scFv_2_known_sequence")
)

write.csv(scFv_1_annotation, file.path(REP_out, "scFv_1_annotation.csv"), row.names = F)
write.csv(scFv_2_annotation, file.path(REP_out, "scFv_2_annotation.csv"), row.names = F)

## (No linker) ----
scFv_1_annotation_nolinker <- findKnownSequence(
  df = contig_1_filter4_nolinker %>% select(-any_of(cols_to_remove)),
  known_table = known_sequences,
  sequence_col = "insert_seq",
  id_col = "contig_id",
  min_score = 90,
  min_pident = 90,
  keep_reverse_strand = F,
  blast_program = "blastn",
  task = "blastn",
  output_dir = file.path(REP_out, "BLAST_scFv_1_nolinker_known_sequence")
)

scFv_2_annotation_nolinker <- findKnownSequence(
  df = contig_2_filter4_nolinker %>% select(-any_of(cols_to_remove)),
  known_table = known_sequences,
  sequence_col = "insert_seq",
  id_col = "contig_id",
  min_score = 90,
  min_pident = 90,
  keep_reverse_strand = F,
  blast_program = "blastn",
  task = "blastn",
  output_dir = file.path(REP_out, "BLAST_scFv_2_nolinker_known_sequence")
)

write.csv(scFv_1_annotation_nolinker, file.path(REP_out, "scFv_1_annotation_nolinker.csv"), row.names = F)
write.csv(scFv_2_annotation_nolinker, file.path(REP_out, "scFv_2_annotation_nolinker.csv"), row.names = F)

# Summary ----
scFv_1_summary <- repSummarize(
  raw = contig_1,
  UMI_filter = contig_1_filter1,
  franking_site = contig_1_filter2,
  length_filter = contig_1_filter3,
  linker_filter = contig_1_filter4,
  annotation = scFv_1_annotation
)

scFv_2_summary <- repSummarize(
  raw = contig_2,
  UMI_filter = contig_2_filter1,
  franking_site = contig_2_filter2,
  length_filter = contig_2_filter3,
  linker_filter = contig_2_filter4,
  annotation = scFv_2_annotation
)

write.csv(scFv_1_summary, file.path(REP_out, "scFv_1_summary.csv"), row.names = F)
write.csv(scFv_2_summary, file.path(REP_out, "scFv_2_summary.csv"), row.names = F)

# Combine ----
scFv_comb_annotation <- bind_rows(
  list(
    scFv_1_annotation %>% mutate(batch_id = "batch_1"),
    scFv_2_annotation %>% mutate(batch_id = "batch_2")
  )) %>% select(
    "sample_id", "batch_id", "contig_id", "barcode",
    "is_cell", "high_confidence", "reads", "umis",
    "insert_seq_length", "linker_name", "scFv_type",
    "sequence",
    "insert_seq", "insert_seq_quals",
    "V1_seq", "linker_seq", "V2_seq"
  )
scFv_comb_annotation %>%
  dplyr::count(sample_id, batch_id, contig_id, name = "n") %>%
  dplyr::filter(n > 1)
as_tibble(scFv_comb_annotation)

# IgBLAST ---- 
scFv_comb_annotation_ig <- annotate_scfv_igblast(
  scFvs = scFv_comb_annotation,
  ids = c("sample_id", "batch_id", "contig_id"),
  work_dir = file.path(REP_out, "IgBLAST_scFv_comb"),
  threads = 4
)

# names(scFv_comb_annotation_ig)
# table(scFv_comb_annotation_ig$is_cell)
# table(scFv_comb_annotation_ig$high_confidence)
# unique(scFv_comb_annotation_ig$V1_V_call)
# unique(scFv_comb_annotation_ig$V2_V_call)
# table(nchar(scFv_comb_annotation_ig$V1_CDR3))
# table(nchar(scFv_comb_annotation_ig$V2_CDR3))

## Post-IgBLAST cleaning ----
### Step1: Internal QC ----
scFv_comb_annotation_ig_clean_1 <- igblast_qc(
  scFvs = scFv_comb_annotation_ig,
  keep_qc = F,
  keep_fail = F
)
names(scFv_comb_annotation_ig_clean_1)

### Step2: Sequence deduplication ----
scFv_comb_annotation_ig_clean_2 <- deduplicate_exact_scfv(scFv_comb_annotation_ig_clean_1, min_reads = 100)
scFv_comb_annotation_ig_clean_2 %>% # Check duplicates
  dplyr::count(sample_id, barcode, insert_seq, name = "n") %>%
  dplyr::filter(n > 1)

### Step3: NGS error correction ----
scFv_comb_annotation_ig_clean_3 <- scfv_error_correction(
  scFvs = scFv_comb_annotation_ig_clean_2,
  max_nt_mismatches = 10,
  min_reads = 100,
  min_major_mismatch_q = 50,
  max_minor_mismatch_q = 40
)
correction_report <- attr(scFv_comb_annotation_ig_clean_3, "correction_report")
print(correction_report, width = 2000)

### Step4: Low Q filter ----
# filter low min_q sequences
scFv_comb_annotation_ig_clean_4 <- scFv_comb_annotation_ig_clean_3 %>% filter(insert_min_q >= 20)

### Step5: Clonal analysis ----
clone_result_comb <- assign_clone_id(
  scFvs = scFv_comb_annotation_ig_clean_4,
  clone_registry = NULL,
  info_columns = c(
    "scFv_type", "host", "chain_config", "linker_name", 
    "V1_V_call", "V1_D_call", "V1_J_call",
    "V2_V_call", "V2_D_call", "V2_J_call"),
  keep_key = TRUE
)

scFv_comb_annotation_ig_clean_5 <- clone_result_comb$annotated
clone_registry <- clone_result_comb$registry
write.csv(clone_registry, file.path(trainingData_out, "clone_registry.csv"), row.names = F)

### Step6: Classify barcode-level scFvs ----
scFv_comb_annotation_ig_clean_6 <- classify_scfv_barcode(scFv_comb_annotation_ig_clean_5) %>%
  select(c(-insert_seq_quals, clone_key))
scFv_comb_annotation_ig_clean_6 %>% # Count Single / Multiple
  dplyr::distinct(sample_id, barcode, clone_count, class) %>%
  dplyr::count(class)
scFv_comb_annotation_ig_clean_6_multiple <- scFv_comb_annotation_ig_clean_6 %>%
  dplyr::filter(class %in% "Multiple")
write.csv(scFv_comb_annotation_ig_clean_6, file.path(REP_out, "scFv_comb_annotation_ig_clean.csv"), row.names = F)

# Link scFv Data with Seurat obj ----
if (!exists("seu_sub")) seu_sub <- readRDS("seuratObjects/seu_sub.rds")
geneset_JAKSTAT <- read.csv("other/JAK_STAT_genes.csv") %>% pull(1)

write.csv(
  data.frame(geneset_JAKSTAT = geneset_JAKSTAT),
  file.path(trainingData_out, "geneset_usedGeneList.csv"),
  row.names = FALSE,
  na = ""
)

seu_meta <- seu_sub[[]] %>% rownames_to_column("barcode")
write.csv(seu_meta, file.path(trainingData_out, "metadata.csv"), row.names = F)

## Generate training data set ----
training_data <- collectTrainingData(
  scFv_annotation = scFv_comb_annotation_ig_clean_6,
  seu_obj = seu_sub,
  metadata = c(
    "orig.ident", "clusters_v1", "Gene set JAKSTAT"
  ),
  features = geneset_JAKSTAT,
  assay = "RNA",
  layer = "data",
  keep_cell_without_exp = F
)

# Test multiple
training_data %>%
  count(barcode, name = "barcode_row_counts") %>%
  count(barcode_row_counts, sort = TRUE)

write.csv(training_data, file.path(trainingData_out, "training_data.csv"), row.names = F)

## Append to Seurat metadata ----
scFv_annotation_seu <- training_data %>%
  select(barcode, scFv_type) %>%
  distinct() %>%
  group_by(barcode) %>%
  summarise(
    scFv_type = if (n_distinct(scFv_type) > 1) {
      "Multiple"
    } else {
      first(scFv_type)
    },
    .groups = "drop"
  )

seu_sub <- AppendSeuratMetadata(seu_sub, scFv_annotation_seu, fill_na_cols = "scFv_type")

seu_data <- LayerData(
  object = seu_sub,
  assay = "RNA",
  layer = "data",
  features = VariableFeatures(FindVariableFeatures(seu_sub, assay = "SCT", nfeatures = 1000)),
  cells = rownames(seu_sub[[]] %>% filter(!scFv_type %in% "Undefined"))
) %>% t() %>% as.data.frame() %>% rownames_to_column("barcode")
write.csv(seu_data, file.path(trainingData_out, "scFv_seurat_data.csv"), row.names = F)

p <- DimPlot(
  seu_sub, 
  group.by = "scFv_type", 
  pt.size = 1, 
  order = c("Novel", "1C4_LH", "STM073_LH", "STM418_HL", "STM418_LH", "STM004_HL", "STM004_LH", "EGFR_HL", "Multiple", "Undefined"),
  cols = c(
    "Novel" = "#F8766D",
    "1C4_LH" = "#B79F00",
    "EGFR_HL" = "#619CFF",
    "STM073_LH" = "#ffa600",
    "STM418_HL" = "#8B56ED",
    "STM418_LH" = "#8B56ED",
    "STM004_HL" = "#00BA38",
    "STM004_LH" = "#00BA38",
    "Multiple" = "#F564E3",
    "Undefined" = "gray90"
  )
) + 
  labs(title = "UMAP: scFv") + 
  ggPreset("DimPlot")
p <- improveUMAP(p)
saveSVG(p, file.path(out, "UMAP_scFv_type.svg"), width = 5, height = 5.8)

# Save RDS ----
saveRDS2(seu_sub, file.path(seu_obj_out, "seu_sub_anno.rds"), overwrite = F)

capture.output(sessionInfo(), file = "sessionInfo_REP.txt")