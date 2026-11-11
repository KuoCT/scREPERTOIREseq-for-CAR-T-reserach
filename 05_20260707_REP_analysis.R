# Set up environment ----
# Set working directory
working_dir <- file.path(this.path::here(),"..")
setwd(working_dir)

# Load libraries
# GEX
library(Seurat)
options(Seurat.object.assay.version = "v3")
# Plots
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
GEX_out <- file.path(getwd(), "out", exp_name, "GEX")
REP_out <- file.path(getwd(), "out", exp_name, "REP")
trainingData_out <- file.path(getwd(), "trainingData")

dir.create(seu_obj_out, showWarnings = F, recursive = T)
dir.create(GEX_out, showWarnings = F, recursive = T)
dir.create(REP_out, showWarnings = F, recursive = T)
dir.create(trainingData_out, showWarnings = F, recursive = T)

# Seed
set.seed(1234)

# Load functions
source("scripts/04_seuratFunctions.R")
source("scripts/04_toolFunctions.R")
source("scripts/05_scRepertoireTool.R")

# Import data ----
if (!exists("scFv_1")) scFv_1 <- readContig("cellrangerData/20260630_scREPERTOIREseq/scFv_PDL1MoCARJNG1_23H3TJLT4/")
if (!exists("scFv_2")) scFv_2 <- readContig("cellrangerData/20260630_scREPERTOIREseq/scFv_PDL1MoCARJNG2_23H3TJLT4/")

lapply(scFv_1[c("is_cell", "high_confidence")], table)
lapply(scFv_2[c("is_cell", "high_confidence")], table)

write.csv(scFv_1, file.path(REP_out, "scFv_1.csv"), row.names = F)
write.csv(scFv_2, file.path(REP_out, "scFv_2.csv"), row.names = F)

# UMI QC (filter 1) ----
scFv_1_filter1 <- scFv_1 %>% filter(umis >= 3)
scFv_2_filter1 <- scFv_2 %>% filter(umis >= 3)

write.csv(scFv_1_filter1, file.path(REP_out, "scFv_1_filter1_UMI.csv"), row.names = F)
write.csv(scFv_2_filter1, file.path(REP_out, "scFv_2_filter1_UMI.csv"), row.names = F)

# Franking site QC (filter 2) ----
franking_table <- data.frame(
  query_id = c("Myc_tag", "3_tail"),
  sequence = c(
    "GAGCAGAAACTTATTTCAGAGGAAGACCTG",
    "GGACCGACCACGACGCCAGCG"
  )
)

scFv_1_frankingSite <- findFrankingSite(
  df = scFv_1_filter1,
  franking_table = franking_table,
  sequence_col = "sequence",
  id_col = "contig_id",
  blast_program = "blastn",
  min_score = 90, 
  min_pident = 90,
  keep_reverse_strand = F,
  output_dir = file.path(REP_out, "BLAST_scFv_1_frankingSite")
)

scFv_2_frankingSite <- findFrankingSite(
  df = scFv_2_filter1,
  franking_table = franking_table,
  sequence_col = "sequence",
  id_col = "contig_id",
  blast_program = "blastn",
  min_score = 90, 
  min_pident = 90,
  keep_reverse_strand = F,
  output_dir = file.path(REP_out, "BLAST_scFv_2_frankingSite")
)

# franking_site_hits <- attr(scFv_1_frankingSite, "blast_hits")
# franking_site_hits <- attr(scFv_2_frankingSite, "blast_hits")

scFv_1_filter2 <- cutFrankingSite(scFv_1_frankingSite)
scFv_2_filter2 <- cutFrankingSite(scFv_2_frankingSite)

write.csv(scFv_1_filter2, file.path(REP_out, "scFv_1_filter2_frankingSite.csv"), row.names = F)
write.csv(scFv_2_filter2, file.path(REP_out, "scFv_2_filter2_frankingSite.csv"), row.names = F)

# scFv length QC (filter 3) ----
newPlotSet(TRUE)
p <- seqLengthPlot(scFv_1_filter2, colors = c("#F8766D", "black", "#D04132", "#7B1400")) + 
  coord_cartesian(xlim = c(0, 1000)) + 
  geom_vline(xintercept = c(600, 900), linetype = "dashed") +
  labs(title = "scFv Lengths", subtitle = "Batch 1") + 
  ggPreset("default")
p %>% addPlot("Batch 1")
saveSVG(p, file.path(REP_out, "scFv_1_lengths.svg"), width = 5, height = 4)

p <- seqLengthPlot(scFv_2_filter2, colors = c("#00BFC4", "black", "#038C90", "#044A4D")) + 
  coord_cartesian(xlim = c(0, 1000)) + 
  geom_vline(xintercept = c(600, 900), linetype = "dashed") +
  labs(title = "scFv Lengths", subtitle = "Batch 2") + 
  ggPreset("default")
p %>% addPlot("Batch 2")
saveSVG(p, file.path(REP_out, "scFv_2_lengths.svg"), width = 5, height = 4)

p <- seqLengthPlot(scFv_1_filter2 %>% filter(insert_seq_length > 100), colors = c("#F8766D", "black", "#D04132", "#7B1400")) + 
  coord_cartesian(xlim = c(0, 1000)) + 
  geom_vline(xintercept = c(600, 900), linetype = "dashed") +
  labs(title = "scFv Lengths", subtitle = "Batch 1 narrowed") + 
  ggPreset("default")
p %>% addPlot("Batch 1 narrowed")
saveSVG(p, file.path(REP_out, "scFv_1_lengths_v2.svg"), width = 5, height = 4)

p <- seqLengthPlot(scFv_2_filter2 %>% filter(insert_seq_length > 100), colors = c("#00BFC4", "black", "#038C90", "#044A4D")) + 
  coord_cartesian(xlim = c(0, 1000)) + 
  geom_vline(xintercept = c(600, 900), linetype = "dashed") +
  labs(title = "scFv Lengths", subtitle = "Batch 2 narrowed") + 
  ggPreset("default")
p %>% addPlot("Batch 2 narrowed")
saveSVG(p, file.path(REP_out, "scFv_2_lengths_v2.svg"), width = 5, height = 4)

listPlotSet()
p <- wrap_plots(PlotSet(), ncol = 2)
savePNG(p, file.path(REP_out, "scFv_lengthPlots.png"), width = 8, height = 6.5)
closePlotSet()

cols_to_remove <- c("Myc_tag_count", "Myc_tag_position", "Myc_tag_score", 
                    "3_tail_count", "3_tail_position", "3_tail_score", 
                    "left_hit_start", "left_hit_end", "left_hit_score", 
                    "right_hit_start", "right_hit_end", "right_hit_score", "sequence_length", "valid_franking_order")

scFv_1_filter3 <- scFv_1_filter2 %>% filter(insert_seq_length <= 900, insert_seq_length >= 600) %>% select(-all_of(cols_to_remove))
scFv_2_filter3 <- scFv_2_filter2 %>% filter(insert_seq_length <= 900, insert_seq_length >= 600) %>% select(-all_of(cols_to_remove))

write.csv(scFv_1_filter3, file.path(REP_out, "scFv_1_filter3_scFvLength.csv"), row.names = F)
write.csv(scFv_2_filter3, file.path(REP_out, "scFv_2_filter3_scFvLength.csv"), row.names = F)

# Find linker QC (filter 4)----
linkers <- read.csv("other/linker_list.csv")

scFv_1_linker <- findLinker(
  df = scFv_1_filter3,
  linker_table = linkers,
  sequence_col = "insert_seq",
  id_col = "contig_id",
  min_score = 90,
  min_pident = 90,
  keep_reverse_strand = FALSE,
  blast_program = "blastn",
  task = "blastn",
  dust = "no",
  output_dir = file.path(REP_out, "BLAST_scFv_1_linker")
)

scFv_2_linker <- findLinker(
  df = scFv_2_filter3,
  linker_table = linkers,
  sequence_col = "insert_seq",
  id_col = "contig_id",
  min_score = 80,
  min_pident = 80,
  keep_reverse_strand = FALSE,
  blast_program = "blastn",
  task = "blastn",
  dust = "no",
  output_dir = file.path(REP_out, "BLAST_scFv_2_linker")
)

# linker_hits <- attr(scFv_1_linker, "blast_hits")
# linker_hits <- attr(scFv_2_linker, "blast_hits")

scFv_1_filter4 <- scFv_1_linker %>% dplyr::filter(!is.na(linker_name)) %>% cutLinker()
scFv_2_filter4 <- scFv_2_linker %>% dplyr::filter(!is.na(linker_name)) %>% cutLinker()

write.csv(scFv_1_filter4, file.path(REP_out, "scFv_1_filter4_linker.csv"), row.names = F)
write.csv(scFv_2_filter4, file.path(REP_out, "scFv_2_filter4_linker.csv"), row.names = F)

scFv_1_nolinker <- scFv_1_linker %>% dplyr::filter(is.na(linker_name))
scFv_2_nolinker <- scFv_2_linker %>% dplyr::filter(is.na(linker_name))

write.csv(scFv_1_nolinker, file.path(REP_out, "scFv_1_filter4_nolinker.csv"), row.names = F)
write.csv(scFv_2_nolinker, file.path(REP_out, "scFv_2_filter4_nolinker.csv"), row.names = F)

# Find known sequences ----
known_sequences <- read.csv("other/known_CAR_sequence.csv") %>%
  dplyr::mutate(sequence = paste0(V_left, linker, V_right))

cols_to_remove <- c("linker_position", "linker_score", "linker_start", "linker_end")
scFv_1_annotation <- findKnownSequence(
  df = scFv_1_filter4 %>% select(-all_of(cols_to_remove)),
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
  df = scFv_2_filter4 %>% select(-all_of(cols_to_remove)),
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

# Summary ----
scFv_1_summary <- repSummarize(
  raw = scFv_1,
  UMI_filter = scFv_1_filter1,
  franking_site = scFv_1_filter2,
  length_filter = scFv_1_filter3,
  linker_filter = scFv_1_filter4,
  annotation = scFv_1_annotation
)

scFv_2_summary <- repSummarize(
  raw = scFv_2,
  UMI_filter = scFv_2_filter1,
  franking_site = scFv_2_filter2,
  length_filter = scFv_2_filter3,
  linker_filter = scFv_2_filter4,
  annotation = scFv_2_annotation
)

write.csv(scFv_1_summary, file.path(REP_out, "scFv_1_summary.csv"), row.names = F)
write.csv(scFv_2_summary, file.path(REP_out, "scFv_2_summary.csv"), row.names = F)

# Combine ----
scFv_batch_combine <- bind_rows(
  list(
    scFv_1_annotation,
    scFv_2_annotation
  )
) %>% select(
  "sample_id", 
  "contig_id", 
  "barcode",
  "is_cell",
  "high_confidence",
  "reads", 
  "umis", 
  "insert_seq_length", 
  "linker_name", 
  "sequence_name",
  "sequence",
  "insert_seq",
  "insert_V1_seq", 
  "insert_linker_seq", 
  "insert_V2_seq"
)
write.csv(scFv_batch_combine, file.path(REP_out, "scFv_batch_combine.csv"), row.names = F)
as_tibble(scFv_batch_combine)
str(scFv_batch_combine)

# IgBLAST ---- 
# scFv_annotation_ig <- scFv_batch_combine %>% 
#   select(
#     barcode,
#     sequence_name,
#     insert_V1_seq,
#     insert_linker_seq,
#     insert_V2_seq
#   ) %>%
#   distinct() %>%
#   rename(scFv_type = sequence_name)

scFv_annotation_ig <- annotate_scfv_igblast(
  scFvs = scFv_batch_combine %>% rename(scFv_type = sequence_name),
  work_dir = file.path(REP_out, "IgBLAST"),
  threads = 4
)

as_tibble(scFv_annotation_ig)
table(scFv_annotation_ig$is_cell)
table(scFv_annotation_ig$high_confidence)
unique(scFv_annotation_ig$V1_V_call)
unique(scFv_annotation_ig$V2_V_call)
table(nchar(scFv_annotation_ig$V1_CDR3))
table(nchar(scFv_annotation_ig$V2_CDR3))

## Post IgBLAST cleaning ----
### Step1: Internal QC ----
scfv_igblast_exception_table <- tribble(
  ~scFv_type,  ~slot, ~allow_missing_cdr3, ~allow_productive_false,
  ~allow_vj_in_frame_na, ~allow_vj_in_frame_false,
  ~allow_frameshift, ~allow_stop_codon, ~reason,
  
  "STM004_HL", "V1", TRUE, TRUE,
  TRUE, FALSE,
  FALSE, FALSE,
  "Known engineered construct; V1_CDR3 is not reliably recognized by IgBLAST."
)

scFv_annotation_ig_qc <- scFv_annotation_ig %>%
  apply_igblast_qc_with_exceptions(
    exception_table = scfv_igblast_exception_table
  ) %>%
  apply_insert_translation_qc()

scFv_annotation_ig_qc %>%
  count(
    pass_igblast_qc,
    pass_insert_translation_qc
  )

write.csv(scFv_annotation_ig_qc, file.path(REP_out, "scFv_annotation_IgBLAST.csv"), row.names = F)

scFv_annotation_ig_qc <- scFv_annotation_ig_qc %>%
  filter(
    pass_igblast_qc,
    pass_insert_translation_qc
  )

### Step2: Clean contig-level working table ----
scFv_contig_work <- scFv_annotation_ig_qc %>%
  mutate(
    contig_row_id = row_number(),
    
    is_cell_bool = to_logical_10x(is_cell),
    high_confidence_bool = to_logical_10x(high_confidence),
    
    scFv_class = case_when(
      scFv_type == "Novel" ~ "Novel",
      !is.na(scFv_type) & scFv_type != "Novel" ~ "Known",
      TRUE ~ "Unassigned"
    ),
    
    V1_V_gene = first_gene_call(V1_V_call),
    V1_D_gene = first_gene_call(V1_D_call),
    V1_J_gene = first_gene_call(V1_J_call),
    V1_C_gene = first_gene_call(V1_C_call),
    
    V2_V_gene = first_gene_call(V2_V_call),
    V2_D_gene = first_gene_call(V2_D_call),
    V2_J_gene = first_gene_call(V2_J_call),
    V2_C_gene = first_gene_call(V2_C_call),
    
    V1_CDR3_nt_length = nchar(V1_CDR3),
    V2_CDR3_nt_length = nchar(V2_CDR3),
    V1_CDR3_aa_length = nchar(V1_CDR3_aa),
    V2_CDR3_aa_length = nchar(V2_CDR3_aa),
    
    known_clone_id = case_when(
      scFv_class == "Known" ~ scFv_type,
      TRUE ~ NA_character_
    ),
    
    novel_clone_key = case_when(
      scFv_class == "Novel" ~ make_novel_clone_key(
        host = host,
        chain_config = chain_config,
        V1_chain_type = V1_chain_type,
        V1_V_gene = V1_V_gene,
        V1_J_gene = V1_J_gene,
        V1_CDR3_aa = V1_CDR3_aa,
        V1_CDR3 = V1_CDR3,
        V2_chain_type = V2_chain_type,
        V2_V_gene = V2_V_gene,
        V2_J_gene = V2_J_gene,
        V2_CDR3_aa = V2_CDR3_aa,
        V2_CDR3 = V2_CDR3
      ),
      TRUE ~ NA_character_
    ),
    
    clone_id = case_when(
      scFv_class == "Known" ~ known_clone_id,
      scFv_class == "Novel" ~ novel_clone_key,
      TRUE ~ NA_character_
    ),
    
    candidate_id = clone_id,
    
    quality_score =
      1000 * as.integer(pass_igblast_qc) +
      1000 * as.integer(pass_insert_translation_qc) +
      100 * as.integer(high_confidence_bool %in% TRUE) +
      50 * as.integer(is_cell_bool %in% TRUE) +
      10 * log10(umis + 1) +
      log10(reads + 1)
  ) %>%
  select(
    contig_row_id,
    
    sample_id,
    barcode,
    contig_id,
    
    scFv_type,
    scFv_class,
    clone_id,
    candidate_id,
    known_clone_id,
    novel_clone_key,
    
    host,
    chain_config,
    linker_name,
    
    is_cell,
    high_confidence,
    is_cell_bool,
    high_confidence_bool,
    reads,
    umis,
    quality_score,
    
    insert_seq_length,
    insert_seq,
    insert_seq_aa,
    insert_V1_seq,
    insert_linker_seq,
    insert_V2_seq,
    
    pass_igblast_qc,
    pass_insert_translation_qc,
    igblast_exception_applied,
    igblast_exception_reason,
    insert_translation_qc_fail_reason,
    
    V1_chain_type,
    V1_frameshift,
    V1_productive,
    V1_stop_codon,
    V1_vj_in_frame,
    V1_v_frameshift,
    V1_j_frameshift,
    V1_pass_igblast_qc,
    V1_CDR3_ok,
    V1_CDR3_aa_ok,
    V1_exception_applied = V1_igblast_exception_applied,
    
    V1_V_call,
    V1_D_call,
    V1_J_call,
    V1_C_call,
    V1_V_gene,
    V1_D_gene,
    V1_J_gene,
    V1_C_gene,
    V1_CDR3,
    V1_CDR3_aa,
    V1_CDR3_nt_length,
    V1_CDR3_aa_length,
    
    V2_chain_type,
    V2_frameshift,
    V2_productive,
    V2_stop_codon,
    V2_vj_in_frame,
    V2_v_frameshift,
    V2_j_frameshift,
    V2_pass_igblast_qc,
    V2_CDR3_ok,
    V2_CDR3_aa_ok,
    V2_exception_applied = V2_igblast_exception_applied,
    
    V2_V_call,
    V2_D_call,
    V2_J_call,
    V2_C_call,
    V2_V_gene,
    V2_D_gene,
    V2_J_gene,
    V2_C_gene,
    V2_CDR3,
    V2_CDR3_aa,
    V2_CDR3_nt_length,
    V2_CDR3_aa_length,
    
    V1_alignment_aa_in_insert,
    V2_alignment_aa_in_insert,
    V1_alignment_aa_found_in_other_frame,
    V2_alignment_aa_found_in_other_frame,
    
    any_of(c(
      "V1_v_identity",
      "V1_j_identity",
      "V2_v_identity",
      "V2_j_identity",
      "V1_hit_score",
      "V2_hit_score"
    ))
  )

# Check row counts per barcode
scFv_contig_work %>%
  count(barcode, name = "n_rows_per_barcode") %>%
  count(n_rows_per_barcode, sort = TRUE)

# Check candidate counts per barcode
scFv_contig_work %>%
  distinct(barcode, candidate_id, scFv_class) %>%
  count(barcode, name = "n_candidates_per_barcode") %>%
  count(n_candidates_per_barcode, sort = TRUE)

### Step3: Candidate-level table ----
scFv_candidate_work <- scFv_contig_work %>%
  group_by(barcode, candidate_id) %>%
  arrange(
    desc(quality_score),
    desc(umis),
    desc(reads),
    .by_group = TRUE
  ) %>%
  summarise(
    scFv_type = first(scFv_type),
    scFv_class = first(scFv_class),
    clone_id = first(clone_id),
    known_clone_id = first(known_clone_id),
    novel_clone_key = first(novel_clone_key),
    
    representative_contig_row_id = first(contig_row_id),
    representative_contig_id = first(contig_id),
    representative_sample_id = first(sample_id),
    
    n_rows_in_candidate = n(),
    n_samples_detected = n_distinct(sample_id),
    samples_detected = paste(sort(unique(sample_id)), collapse = ";"),
    
    total_reads = sum(reads, na.rm = TRUE),
    max_umis = max(umis, na.rm = TRUE),
    sum_umis = sum(umis, na.rm = TRUE),
    
    representative_reads = first(reads),
    representative_umis = first(umis),
    representative_quality_score = first(quality_score),
    
    insert_seq_length = first(insert_seq_length),
    insert_seq = first(insert_seq),
    insert_seq_aa = first(insert_seq_aa),
    insert_V1_seq = first(insert_V1_seq),
    insert_linker_seq = first(insert_linker_seq),
    insert_V2_seq = first(insert_V2_seq),
    
    host = first(host),
    chain_config = first(chain_config),
    linker_name = first(linker_name),
    
    any_is_cell = any(is_cell_bool %in% TRUE, na.rm = TRUE),
    any_high_confidence = any(high_confidence_bool %in% TRUE, na.rm = TRUE),
    all_high_confidence = all(high_confidence_bool %in% TRUE, na.rm = TRUE),
    
    V1_chain_type = first(V1_chain_type),
    V1_V_gene = first(V1_V_gene),
    V1_J_gene = first(V1_J_gene),
    V1_CDR3 = first(V1_CDR3),
    V1_CDR3_aa = first(V1_CDR3_aa),
    V1_CDR3_nt_length = first(V1_CDR3_nt_length),
    V1_CDR3_aa_length = first(V1_CDR3_aa_length),
    
    V2_chain_type = first(V2_chain_type),
    V2_V_gene = first(V2_V_gene),
    V2_J_gene = first(V2_J_gene),
    V2_CDR3 = first(V2_CDR3),
    V2_CDR3_aa = first(V2_CDR3_aa),
    V2_CDR3_nt_length = first(V2_CDR3_nt_length),
    V2_CDR3_aa_length = first(V2_CDR3_aa_length),
    
    .groups = "drop"
  )

### Step4: Barcode-level summary ----
scFv_barcode_summary <- scFv_candidate_work %>%
  group_by(barcode) %>%
  mutate(
    barcode_total_reads = sum(total_reads, na.rm = TRUE),
    barcode_total_max_umis = sum(max_umis, na.rm = TRUE),
    
    read_fraction = total_reads / barcode_total_reads,
    max_umi_fraction = max_umis / barcode_total_max_umis,
    
    rank_by_reads = dense_rank(desc(total_reads)),
    rank_by_max_umis = dense_rank(desc(max_umis))
  ) %>%
  summarise(
    n_candidates = n_distinct(candidate_id),
    n_known_candidates = n_distinct(candidate_id[scFv_class == "Known"]),
    n_novel_candidates = n_distinct(candidate_id[scFv_class == "Novel"]),
    
    n_total_rows = sum(n_rows_in_candidate, na.rm = TRUE),
    n_samples_detected = n_distinct(unlist(strsplit(samples_detected, ";"))),
    
    top_candidate_id = candidate_id[which.max(total_reads)],
    top_scFv_type = scFv_type[which.max(total_reads)],
    top_scFv_class = scFv_class[which.max(total_reads)],
    
    top_reads = max(total_reads, na.rm = TRUE),
    total_reads = sum(total_reads, na.rm = TRUE),
    top_read_fraction = top_reads / total_reads,
    
    top_max_umis = max(max_umis, na.rm = TRUE),
    total_max_umis = sum(max_umis, na.rm = TRUE),
    top_max_umi_fraction = top_max_umis / total_max_umis,
    
    candidate_ids = paste(candidate_id, collapse = ";"),
    scFv_types = paste(scFv_type, collapse = ";"),
    scFv_classes = paste(scFv_class, collapse = ";"),
    
    .groups = "drop"
  ) %>%
  mutate(
    barcode_scfv_status = case_when(
      n_candidates == 1 & top_scFv_class == "Known" ~
        "Single_known_scFv",
      
      n_candidates == 1 & top_scFv_class == "Novel" ~
        "Single_novel_scFv",
      
      n_candidates >= 2 &
        n_known_candidates >= 2 &
        n_novel_candidates == 0 ~
        "Multiple_known_scFv",
      
      n_candidates >= 2 &
        n_known_candidates == 0 &
        n_novel_candidates >= 2 ~
        "Multiple_novel_scFv",
      
      n_candidates >= 2 &
        n_known_candidates >= 1 &
        n_novel_candidates >= 1 ~
        "Multiple_known_plus_novel_scFv",
      
      TRUE ~ "Ambiguous_scFv"
    )
  )

scFv_barcode_summary %>%
  count(barcode_scfv_status, sort = TRUE)

scFv_barcode_summary %>%
  count(n_candidates, sort = TRUE)

scFv_barcode_summary %>%
  filter(n_candidates >= 2) %>%
  arrange(desc(n_candidates), desc(top_read_fraction))


### Step5: Final annotation table ----
# Build CDR3 set while preserving V1/V2 order.
make_cdr3_set <- function(V1_CDR3, V2_CDR3) {
  V1_CDR3 <- dplyr::coalesce(V1_CDR3, "NA")
  V2_CDR3 <- dplyr::coalesce(V2_CDR3, "NA")
  
  dplyr::if_else(
    V1_CDR3 == "NA" & V2_CDR3 == "NA",
    NA_character_,
    paste(V1_CDR3, V2_CDR3, sep = "__")
  )
}

# Translate a DNA segment in the same frame used for insert-level translation QC.
translate_segment_nt <- function(sequence) {
  vapply(
    sequence,
    translate_nt,
    character(1),
    frame_start = 1
  )
}

# Create stable novel clone IDs.
# V1/V2 order is preserved because novel_clone_key was built as V1 followed by V2.
# If the same AA-level clone has multiple nucleotide variants, they share one clone_ID.
novel_clone_id_map <- scFv_candidate_work %>%
  filter(scFv_class == "Novel") %>%
  mutate(
    clone_CDR3_set = make_cdr3_set(V1_CDR3, V2_CDR3),
    clone_CDR3_set_aa = make_cdr3_set(V1_CDR3_aa, V2_CDR3_aa)
  ) %>%
  group_by(novel_clone_key) %>%
  arrange(
    desc(total_reads),
    desc(max_umis),
    desc(n_samples_detected),
    .by_group = TRUE
  ) %>%
  summarise(
    representative_clone_CDR3_set = first(clone_CDR3_set),
    representative_clone_CDR3_set_aa = first(clone_CDR3_set_aa),
    
    n_candidates_in_novel_clone = n(),
    n_barcodes_in_novel_clone = n_distinct(barcode),
    n_V1_CDR3_nt_variants = n_distinct(V1_CDR3, na.rm = TRUE),
    n_V2_CDR3_nt_variants = n_distinct(V2_CDR3, na.rm = TRUE),
    
    host = first(host),
    chain_config = first(chain_config),
    V1_chain_type = first(V1_chain_type),
    V1_V_gene = first(V1_V_gene),
    V1_J_gene = first(V1_J_gene),
    V2_chain_type = first(V2_chain_type),
    V2_V_gene = first(V2_V_gene),
    V2_J_gene = first(V2_J_gene),
    
    .groups = "drop"
  ) %>%
  arrange(
    host,
    chain_config,
    V1_chain_type,
    V1_V_gene,
    V1_J_gene,
    V2_chain_type,
    V2_V_gene,
    V2_J_gene,
    representative_clone_CDR3_set_aa,
    representative_clone_CDR3_set,
    novel_clone_key
  ) %>%
  mutate(
    clone_ID_novel = paste0("novel_", row_number())
  ) %>%
  select(
    novel_clone_key,
    clone_ID_novel,
    representative_clone_CDR3_set,
    representative_clone_CDR3_set_aa,
    n_candidates_in_novel_clone,
    n_barcodes_in_novel_clone,
    n_V1_CDR3_nt_variants,
    n_V2_CDR3_nt_variants
  )

novel_clone_id_map %>%
  count(novel_clone_key, name = "n_map_rows") %>%
  filter(n_map_rows > 1) # Should be zero row

# Add candidate-level clone identity.
scFv_candidate_final_info <- scFv_candidate_work %>%
  left_join(
    novel_clone_id_map,
    by = "novel_clone_key",
    relationship = "many-to-one"
  ) %>%
  left_join(
    scFv_barcode_summary %>%
      select(
        barcode,
        n_candidates,
        n_known_candidates,
        n_novel_candidates,
        barcode_scfv_status
      ),
    by = "barcode",
    relationship = "many-to-one"
  ) %>%
  mutate(
    clone_ID = case_when(
      scFv_class == "Known" ~ scFv_type,
      scFv_class == "Novel" ~ clone_ID_novel,
      TRUE ~ NA_character_
    ),
    
    clone_CDR3_set = case_when(
      scFv_class == "Novel" &
        !is.na(representative_clone_CDR3_set) ~
        representative_clone_CDR3_set,
      TRUE ~ make_cdr3_set(V1_CDR3, V2_CDR3)
    ),
    
    clone_CDR3_set_aa = case_when(
      scFv_class == "Novel" &
        !is.na(representative_clone_CDR3_set_aa) ~
        representative_clone_CDR3_set_aa,
      TRUE ~ make_cdr3_set(V1_CDR3_aa, V2_CDR3_aa)
    ),
    
    final_scfv_type = case_when(
      n_candidates >= 2 ~ "multiple",
      scFv_class == "Known" ~ scFv_type,
      scFv_class == "Novel" ~ "novel",
      TRUE ~ "unassigned"
    ),
    
    clone_type = case_when(
      n_candidates >= 2 ~ "multiple",
      scFv_class == "Known" ~ "known",
      scFv_class == "Novel" ~ "unknown",
      TRUE ~ "unknown"
    )
  ) %>%
  select(
    barcode,
    candidate_id,
    
    final_scfv_type,
    clone_type,
    clone_ID,
    clone_CDR3_set,
    clone_CDR3_set_aa,
    
    n_candidates,
    n_known_candidates,
    n_novel_candidates,
    barcode_scfv_status,
    
    n_rows_in_candidate,
    n_samples_detected,
    samples_detected,
    
    total_reads,
    max_umis,
    sum_umis,
    
    n_candidates_in_novel_clone,
    n_barcodes_in_novel_clone,
    n_V1_CDR3_nt_variants,
    n_V2_CDR3_nt_variants
  )

nrow(scFv_candidate_work)
nrow(scFv_candidate_final_info)

# Recover full FR/CDR and VDJ-position columns from the full QC table.
# scFv_contig_work intentionally removed many of these columns.
scFv_full_region_cols <- scFv_annotation_ig_qc %>%
  mutate(contig_row_id = row_number()) %>%
  select(
    contig_row_id,
    
    V1_V,
    V1_D,
    V1_J,
    V1_C,
    V1_FR1,
    V1_FR1_aa,
    V1_CDR1,
    V1_CDR1_aa,
    V1_FR2,
    V1_FR2_aa,
    V1_CDR2,
    V1_CDR2_aa,
    V1_FR3,
    V1_FR3_aa,
    V1_CDR3,
    V1_CDR3_aa,
    V1_FR4,
    V1_FR4_aa,
    
    V2_V,
    V2_D,
    V2_J,
    V2_C,
    V2_FR1,
    V2_FR1_aa,
    V2_CDR1,
    V2_CDR1_aa,
    V2_FR2,
    V2_FR2_aa,
    V2_CDR2,
    V2_CDR2_aa,
    V2_FR3,
    V2_FR3_aa,
    V2_CDR3,
    V2_CDR3_aa,
    V2_FR4,
    V2_FR4_aa
  )

scFv_annotation_final <- scFv_contig_work %>%
  left_join(
    scFv_candidate_final_info,
    by = c("barcode", "candidate_id")
  ) %>%
  left_join(
    scFv_full_region_cols,
    by = "contig_row_id",
    suffix = c("", "_full")
  ) %>%
  mutate(
    insert_V1_seq_aa = translate_segment_nt(insert_V1_seq),
    insert_linker_seq_aa = translate_segment_nt(insert_linker_seq),
    insert_V2_seq_aa = translate_segment_nt(insert_V2_seq),
    
    # quality_score is used only to rank rows within the same barcode + clone_ID.
    # Formula:
    #   1000 if high_confidence == TRUE
    # +  500 if is_cell == TRUE
    # +  100 * number of NGS samples where this barcode-candidate was detected
    # +   50 * log10(umis + 1)
    # +   10 * log10(reads + 1)
    #
    # Since all rows here already passed IgBLAST QC and insert translation QC,
    # those pass flags are not included in the ranking score.
    quality_score =
      1000 * as.integer(high_confidence_bool %in% TRUE) +
      500 * as.integer(is_cell_bool %in% TRUE) +
      100 * n_samples_detected +
      50 * log10(umis + 1) +
      10 * log10(reads + 1)
  ) %>%
  transmute(
    barcode,
    
    scFv_type = final_scfv_type,
    host,
    chain_config,
    linker_name,
    
    insert_seq,
    insert_seq_aa,
    
    insert_V1_seq,
    insert_V1_seq_aa,
    insert_linker_seq,
    insert_linker_seq_aa,
    insert_V2_seq,
    insert_V2_seq_aa,
    
    V1_chain_type,
    V1_V,
    V1_D,
    V1_J,
    V1_C,
    V1_V_call,
    V1_D_call,
    V1_J_call,
    V1_C_call,
    V1_V_gene,
    V1_D_gene,
    V1_J_gene,
    V1_C_gene,
    V1_productive,
    V1_stop_codon,
    V1_vj_in_frame,
    V1_v_frameshift,
    V1_j_frameshift,
    V1_frameshift,
    
    V1_FR1,
    V1_FR1_aa,
    V1_CDR1,
    V1_CDR1_aa,
    V1_FR2,
    V1_FR2_aa,
    V1_CDR2,
    V1_CDR2_aa,
    V1_FR3,
    V1_FR3_aa,
    V1_CDR3,
    V1_CDR3_aa,
    V1_FR4,
    V1_FR4_aa,
    
    V2_chain_type,
    V2_V,
    V2_D,
    V2_J,
    V2_C,
    V2_V_call,
    V2_D_call,
    V2_J_call,
    V2_C_call,
    V2_V_gene,
    V2_D_gene,
    V2_J_gene,
    V2_C_gene,
    V2_productive,
    V2_stop_codon,
    V2_vj_in_frame,
    V2_v_frameshift,
    V2_j_frameshift,
    V2_frameshift,
    
    V2_FR1,
    V2_FR1_aa,
    V2_CDR1,
    V2_CDR1_aa,
    V2_FR2,
    V2_FR2_aa,
    V2_CDR2,
    V2_CDR2_aa,
    V2_FR3,
    V2_FR3_aa,
    V2_CDR3,
    V2_CDR3_aa,
    V2_FR4,
    V2_FR4_aa,
    
    clone_type,
    clone_ID,
    clone_CDR3_set,
    clone_CDR3_set_aa,
    
    quality_score,
    
    # Traceability columns.
    sample_id,
    contig_id,
    contig_row_id,
    candidate_id,
    reads,
    umis,
    insert_seq_length,
    n_candidates,
    n_known_candidates,
    n_novel_candidates,
    barcode_scfv_status,
    n_rows_in_candidate,
    n_samples_detected,
    samples_detected,
    total_reads,
    max_umis,
    sum_umis,
    is_cell,
    high_confidence,
    pass_igblast_qc,
    pass_insert_translation_qc,
    igblast_exception_applied,
    igblast_exception_reason
  ) %>%
  mutate(
    # Convert position string to sequence
    V1_FR1 = extract_seq_by_position_vec(insert_V1_seq, V1_FR1),
    V1_FR2 = extract_seq_by_position_vec(insert_V1_seq, V1_FR2),
    V1_FR3 = extract_seq_by_position_vec(insert_V1_seq, V1_FR3),
    V1_FR4 = extract_seq_by_position_vec(insert_V1_seq, V1_FR4),
    
    V2_FR1 = extract_seq_by_position_vec(insert_V2_seq, V2_FR1),
    V2_FR2 = extract_seq_by_position_vec(insert_V2_seq, V2_FR2),
    V2_FR3 = extract_seq_by_position_vec(insert_V2_seq, V2_FR3),
    V2_FR4 = extract_seq_by_position_vec(insert_V2_seq, V2_FR4)
  )

as_tibble(scFv_annotation_final)

scFv_annotation_final %>%
  summarise(
    n_rows = n(),
    n_barcodes = n_distinct(barcode),
    n_clone_IDs = n_distinct(clone_ID, na.rm = TRUE),
    n_known_rows = sum(clone_type == "known", na.rm = TRUE),
    n_unknown_rows = sum(clone_type == "unknown", na.rm = TRUE),
    n_multiple_rows = sum(clone_type == "multiple", na.rm = TRUE)
  )

scFv_annotation_final %>%
  count(scFv_type, clone_type, sort = TRUE)

scFv_annotation_final %>%
  count(n_candidates, barcode_scfv_status, sort = TRUE)

write.csv(
  scFv_annotation_final,
  file.path(REP_out, "scFv_annotation_final.csv"),
  row.names = FALSE
)


### Step6: Cleaned representative table ----
scFv_annotation_clean <- scFv_annotation_final %>%
  group_by(barcode, clone_ID) %>%
  arrange(
    desc(quality_score),
    desc(umis),
    desc(reads),
    .by_group = TRUE
  ) %>%
  slice(1) %>%
  ungroup()

scFv_annotation_clean %>%
  summarise(
    n_rows = n(),
    n_barcodes = n_distinct(barcode),
    n_clone_IDs = n_distinct(clone_ID, na.rm = TRUE)
  )

scFv_annotation_clean %>%
  count(scFv_type, clone_type, sort = TRUE)

scFv_annotation_clean %>%
  count(n_candidates, barcode_scfv_status, sort = TRUE)

scFv_annotation_clean %>%
  filter(scFv_type == "multiple") %>%
  arrange(barcode, clone_ID) %>%
  select(
    barcode,
    scFv_type,
    clone_type,
    clone_ID,
    clone_CDR3_set,
    clone_CDR3_set_aa,
    quality_score,
    reads,
    umis,
    n_candidates,
    barcode_scfv_status,
    sample_id,
    contig_id
  )

scFv_annotation_clean_2 <- scFv_annotation_clean %>%
  mutate(
    scFv_type = case_when(
      is.na(scFv_type) ~ NA_character_,
      stringr::str_to_lower(scFv_type) == "novel" ~ "Novel",
      stringr::str_to_lower(scFv_type) == "multiple" ~ "Multiple",
      stringr::str_to_lower(scFv_type) == "undefined" ~ "Undefined",
      TRUE ~ scFv_type
    )
  ) %>%
  select(
    # Cell / scFv identity
    barcode,
    scFv_type,
    host,
    chain_config,
    linker_name,
    
    # Full construct sequences
    insert_seq,
    insert_seq_aa,
    insert_V1_seq,
    insert_V1_seq_aa,
    insert_linker_seq,
    insert_linker_seq_aa,
    insert_V2_seq,
    insert_V2_seq_aa,
    
    # V1 VDJ annotation
    V1_chain_type,
    V1_V_call,
    V1_D_call,
    V1_J_call,
    V1_C_call,
    V1_V_gene,
    V1_D_gene,
    V1_J_gene,
    V1_C_gene,
    
    # V1 FR/CDR DNA and protein sequences
    V1_FR1,
    V1_FR1_aa,
    V1_CDR1,
    V1_CDR1_aa,
    V1_FR2,
    V1_FR2_aa,
    V1_CDR2,
    V1_CDR2_aa,
    V1_FR3,
    V1_FR3_aa,
    V1_CDR3,
    V1_CDR3_aa,
    V1_FR4,
    V1_FR4_aa,
    
    # V2 VDJ annotation
    V2_chain_type,
    V2_V_call,
    V2_D_call,
    V2_J_call,
    V2_C_call,
    V2_V_gene,
    V2_D_gene,
    V2_J_gene,
    V2_C_gene,
    
    # V2 FR/CDR DNA and protein sequences
    V2_FR1,
    V2_FR1_aa,
    V2_CDR1,
    V2_CDR1_aa,
    V2_FR2,
    V2_FR2_aa,
    V2_CDR2,
    V2_CDR2_aa,
    V2_FR3,
    V2_FR3_aa,
    V2_CDR3,
    V2_CDR3_aa,
    V2_FR4,
    V2_FR4_aa,
    
    # Clone annotation
    clone_type,
    clone_ID,
    clone_CDR3_set,
    clone_CDR3_set_aa,
    
    # Final quality ranking
    quality_score
  )


write.csv(
  scFv_annotation_clean_2,
  file.path(REP_out, "scFv_annotation_clean.csv"),
  row.names = FALSE
)


# Link scFv Data with Seurat obj ----
if (!exists("seu_sub"))  seu_sub <- readRDS("seuratObjects/20260707_seu_sub.rds")
# seu_sub <- readRDS("../20260630_primary_CART/seuratObjects/20260706_seu.rds")
scFv_annotation_clean_2 <- read.csv("../20260630_primary_CART/out/05_20260706_REP_makeExceptions2/REP/scFv_annotation_clean.csv")
geneset_list <- readRDS("other/20260622_geneset_list.rds")

write.csv(
  as.data.frame(lapply(geneset_list, `length<-`, max(lengths(geneset_list)))),
  file.path(trainingData_out, "20260707_geneset_usedGeneList.csv"),
  row.names = FALSE,
  na = ""
)

seu_meta <- seu_sub[[]] %>% rownames_to_column("barcode")
write.csv(seu_meta, file.path(trainingData_out, "20260707_metadata.csv"), row.names = F)

seu_data <- LayerData(
  object = seu_sub,
  assay = "RNA",
  layer = "data",
  features = geneset_list[["Set 6"]]
) %>% t() %>% as.data.frame() %>% select(EGFP, everything()) %>% rownames_to_column("barcode")
write.csv(seu_data, file.path(trainingData_out, "20260707_scFv_seurat_data.csv"), row.names = F)

## Generate training data set ----
training_data <- collectTrainingData(
  scFv_annotation = scFv_annotation_clean_2,
  seu_meatdata = seu_meta,
  seu_data = seu_data,
  keep_cols = c("clusters_v1")
)

# Test multiple
training_data %>%
  count(barcode, name = "n_rows") %>%
  count(n_rows, sort = TRUE)

write.csv(training_data, file.path(trainingData_out, "20260707_training_data.csv"), row.names = F)

## Append to Seurat metadata ----
scFv_annotation_seu <- scFv_annotation_clean_2 %>%
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
unique(seu_sub$scFv_type)
table(seu_sub$scFv_type)

p <- DimPlot(
  seu_sub, 
  group.by = "scFv_type", 
  pt.size = 1, 
  order = c(
    "Novel",
    "STM004_LH",
    "Multiple",
    "STM073_LH",
    "1C4_LH",
    "STM418_HL",
    "STM418_LH",
    "STM004_HL",
    "EGFR_HL",
    "Undefined"
  ),
  cols = c(
    "Undefined" = "gray90",
    "EGFR_HL" = "#25B8FA",
    "Multiple" = "#FF82C7",
    "STM073_LH" = "#00ffbc",
    "Novel" = "#EB4747",
    "1C4_LH" = "#25cc3c",
    "STM418_HL" = "#8B56ED",
    "STM418_LH" = "#8B56ED",
    "STM004_HL" = "#ffa600",
    "STM004_LH" = "#ffa600"
  )
) + 
  labs(title = "UMAP: scFv") + 
  ggPreset("DimPlot")
p <- improveUMAP(p)
saveSVG(p, file.path(REP_out, "UMAP_scFv_type.svg"), width = 4.5, height = 5.2)


# Heatmap ----
matrix <- training_data %>%
  select(barcode, any_of(geneset_list[[6]])) %>%
  mutate(barcode = make.unique(barcode)) %>%
  column_to_rownames("barcode") %>%
  t() %>%
  as.matrix()

matrix_cell_anno <- training_data %>%
  select(c(barcode, clusters_v1, scFv_type, `Gene set 6`)) %>%
  mutate(barcode = make.unique(barcode)) %>%
  column_to_rownames("barcode")

matrix_color_list <- list(
  clusters_v1 = c(
    "C1" = pal_hue()(8)[1],
    "C2" = pal_hue()(8)[2],
    "C3" = pal_hue()(8)[3],
    "C4" = pal_hue()(8)[4],
    "C5" = pal_hue()(8)[5],
    "C6" = pal_hue()(8)[6],
    "C7" = pal_hue()(8)[7],
    "C8" = pal_hue()(8)[8]
  ),
  scFv_type = c(
    "Undefined" = "gray90",
    "EGFR_HL" = "#25B8FA",
    "Multiple" = "#FF82C7",
    "STM073_LH" = "#00ffbc",
    "Novel" = "#EB4747",
    "1C4_LH" = "#25cc3c",
    "STM418_HL" = "#8B56ED",
    "STM418_LH" = "#8B56ED",
    "STM004_HL" = "#ffa600",
    "STM004_LH" = "#ffa600"
  ),
  `Gene set 6` = circlize::colorRamp2(c(-0.5, 0, 0.5), c("navy", "white", "firebrick3"))
)

matrix_col_anno <- HeatmapAnnotation(
  df = matrix_cell_anno, 
  col = matrix_color_list, 
  border = T,
  simple_anno_size = unit(3, "mm"),
  gap = unit(1, "mm"), 
  show_annotation_name = T,
  annotation_label = c(clusters_v1 = "Clusters", scFv_type = "scFv type", `Gene set 6` = "Gene set 6"),
  annotation_name_gp = gpar(fontsize = 8, fontface = "bold"),
  show_legend = F
)

## Colors
viridis_ht <- viridis_pal(option = "B")(100)
breaks <- seq(0, 2.5, length.out = length(viridis_ht))
col_fun <- circlize::colorRamp2(breaks, viridis_ht)

## Plot
ht <- Heatmap(
  matrix,
  show_heatmap_legend = F,
  name = "Expression",
  row_dend_width = unit(5, "mm"),
  column_dend_height = unit(5, "mm"),
  top_annotation = matrix_col_anno,
  show_row_names = T,
  row_names_gp = gpar(fontsize = 8),
  show_column_names = F,
  cluster_columns = T,
  col = col_fun,
  border = T,
  use_raster = T,
  clustering_distance_rows = "pearson",
  clustering_method_rows = "ward.D",
  clustering_distance_columns = "euclidean",
  clustering_method_columns = "ward.D"
)

## Legend list
lgd_list <- packLegend(
  Legend(
    title = "Clusters",
    labels = c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8"),
    legend_gp = gpar(fill = pal_hue()(8), col = "black"),
    labels_gp = gpar(fontsize = 8),
    border = "black"
  ),
  Legend(
    title = "scFv Type",
    labels = c("EGFR", "1C4", "STM073", "STM004", "STM418", "Novel", "Multiple"),
    legend_gp = gpar(fill = c("#25B8FA", "#25cc3c", "#00ffbc", "#ffa600", "#8B56ED", "#EB4747", "#FF82C7"), col = "black"),
    labels_gp = gpar(fontsize = 8),
    border = "black"
  ),
  Legend(
    title = "Gene set 6",
    col_fun = circlize::colorRamp2(c(-0.5, 0, 0.5), c("navy", "white", "firebrick3")),
    labels_gp = gpar(fontsize = 8),
    border = "black"
  ),
  Legend(
    title = "SCT Exp.",
    col_fun = col_fun,
    labels_gp = gpar(fontsize = 8),
    border = "black"
  ),
  direction = "vertical",
  gap = unit(4, "mm")
)

svglite(file.path(GEX_out, "Heatmap_gene_set_6.svg"), width = 7, height = 8)
draw(ht, annotation_legend_list = lgd_list, merge_legend = T)
dev.off()
