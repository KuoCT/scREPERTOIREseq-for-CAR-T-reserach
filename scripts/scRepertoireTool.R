# Reading Cellranger VDJ contig ----
# Split contig ID and update barcode index to match the Seurat object
# Example: AAACCAAAGAATGCAA-1_contig_2 -> barcode = AAACCAAAGAATGCAA-1, contig = 2
# With barcode_index = 3: AAACCAAAGAATGCAA-1_contig_2 -> AAACCAAAGAATGCAA-3_contig_2
parse_contig_id <- function(data, barcode_index = NULL) {
  data <- data %>%
    tidyr::separate(
      contig_id,
      into = c("barcode", "contig"),
      sep = "_contig_",
      remove = FALSE,
      convert = TRUE,
      fill = "right"
    )
  
  # Replace the barcode suffix when Seurat uses a different sample index
  if (!is.null(barcode_index)) {
    data <- data %>%
      dplyr::mutate(
        barcode = sub("-[0-9]+$", paste0("-", barcode_index), barcode),
        contig_id = paste0(barcode, "_contig_", contig)
      )
  }
  
  data
}


readContig <- function(folder, sample_id = NULL, barcode_index = NULL) {
  folder <- as.character(folder)
  if (length(folder) != 1L) stop("readContig() only accepts one folder at a time.")
  if (!dir.exists(folder)) stop("Folder does not exist: ", folder)
  
  # Normalize folder path (accept both outs folder or its above)
  folder <- normalizePath(folder, winslash = "/")
  sample_folder <- sub("/outs$", "", folder)
  outs_folder <- file.path(sample_folder, "outs")
  if (!dir.exists(outs_folder)) stop("Cannot find outs folder: ", outs_folder)
  
  # Add sample ID
  if (is.null(sample_id)) sample_id <- basename(sample_folder)
  
  # Search Cellranger VDJ contig annotation file
  json_file <- file.path(outs_folder, "all_contig_annotations.json")
  if (!file.exists(json_file)) stop("Cannot find: ", json_file)
  
  # Input Cellranger contig JSON data
  contig <- jsonlite::fromJSON(
    json_file,
    flatten = TRUE
  ) %>%
    tibble::as_tibble() %>%
    dplyr::transmute(
      sequence,
      sequence_quals = quals,
      sample_id = sample_id,
      contig_id = contig_name,
      is_cell,
      high_confidence,
      reads = as.integer(read_count),
      umis = as.integer(umi_count)
    ) %>%
    
    # Split contig ID and update barcode index
    parse_contig_id(barcode_index = barcode_index) %>%
    dplyr::select(
      sequence, sequence_quals, sample_id, contig_id, barcode, 
      contig, is_cell, high_confidence, reads, umis
    )
  
  contig
}


# BLAST functions ----
sequenceBlast <- function(
    df,
    query_table,
    sequence_col = "sequence",
    id_col = "contig_id",
    min_score = 70,
    min_pident = 70,
    keep_reverse_strand = TRUE,
    blast_program = c("blastn", "blastp"),
    task = NULL,
    dust = NULL,
    max_target_seqs = 1000000,
    output_dir = NULL
) {
  # Select BLAST engine
  blast_program <- match.arg(blast_program)
  db_type <- if (blast_program == "blastp") "prot" else "nucl"
  if (is.null(task)) task <- paste0(blast_program, "-short")
  
  # Validate inputs
  if (!sequence_col %in% colnames(df)) stop("Missing sequence column: ", sequence_col)
  if (!id_col %in% colnames(df)) stop("Missing id column: ", id_col)
  if (!all(c("query_id", "sequence") %in% colnames(query_table))) {
    stop("query_table must contain query_id + sequence")
  }
  
  # Create a temporary folder unless output_dir is provided
  if (is.null(output_dir)) {
    work_dir <- tempfile("sequenceBlast_")
    on.exit(unlink(work_dir, recursive = TRUE, force = TRUE), add = TRUE)
  } else {
    work_dir <- output_dir
  }
  dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Set working file paths
  subject_fasta <- file.path(work_dir, "subjects.fa")
  query_fasta <- file.path(work_dir, "queries.fa")
  blast_db <- file.path(work_dir, "blast_db")
  blast_out <- file.path(work_dir, "blast.tsv")
  
  # Write subject and query sequences in FASTA format
  writeLines(
    paste0(">", df[[id_col]], "\n", toupper(df[[sequence_col]])),
    subject_fasta
  )
  
  writeLines(
    paste0(">", query_table$query_id, "\n", toupper(query_table$sequence)),
    query_fasta
  )
  
  # Build BLAST database
  makeblastdb_status <- system2(
    "makeblastdb",
    c(
      "-in", shQuote(subject_fasta),
      "-dbtype", db_type,
      "-out", shQuote(blast_db)
    )
  )
  if (makeblastdb_status != 0) stop("makeblastdb failed.")
  
  # Set BLAST output fields and their corresponding R column names
  blast_fields <- c(
    query_id = "qseqid", subject_id = "sseqid", pident = "pident",
    alignment_length = "length", mismatch = "mismatch", gapopen = "gapopen",
    qstart = "qstart", qend = "qend", sstart = "sstart", send = "send",
    evalue = "evalue", bitscore = "bitscore", query_length = "qlen"
  )
  outfmt <- paste(c("6", unname(blast_fields)), collapse = " ")
  
  # Set common BLAST arguments
  blast_args <- c(
    "-task", task,
    "-query", shQuote(query_fasta),
    "-db", shQuote(blast_db)
  )
  
  # Add nucleotide-specific arguments
  if (blast_program == "blastn") {
    strand <- if (keep_reverse_strand) "both" else "plus"
    blast_args <- c(blast_args, "-strand", strand)
    
    if (!is.null(dust)) {
      blast_args <- c(blast_args, "-dust", dust)
    }
  }
  
  # Add output and result limit arguments
  blast_args <- c(
    blast_args,
    "-max_target_seqs", format(max_target_seqs, scientific = FALSE, trim = TRUE),
    "-out", shQuote(blast_out),
    "-outfmt", shQuote(outfmt)
  )
  
  # Run BLAST and capture command output
  blast_output <- suppressWarnings(
    system2(
      blast_program,
      blast_args,
      stdout = TRUE,
      stderr = TRUE
    )
  )
  blast_status <- attr(blast_output, "status")
  
  # Stop when BLAST returns a non-zero exit status
  if (!is.null(blast_status) && blast_status != 0) {
    stop(blast_program, " failed:\n", paste(blast_output, collapse = "\n"), call. = FALSE)
  }
  
  # Return an empty result when no alignment is found
  if (!file.exists(blast_out) || file.info(blast_out)$size == 0) {
    return(
      tibble::tibble(
        subject_id = character(), query_id = character(), position = character(),
        subject_start = integer(), subject_end = integer(), strand = character(),
        score = numeric(), pident = numeric(), alignment_length = integer(),
        mismatch = integer(), gapopen = integer(), evalue = numeric(),
        bitscore = numeric()
      )
    )
  }
  
  # Read BLAST results using the predefined column names
  hits <- read.delim(blast_out, header = FALSE, col.names = names(blast_fields), stringsAsFactors = FALSE)
  
  # Format alignment position and calculate normalized alignment score
  hits <- hits %>%
    dplyr::mutate(
      subject_start = pmin(sstart, send),
      subject_end = pmax(sstart, send),
      strand =  if (blast_program == "blastn") ifelse(sstart <= send, "+", "-") else NA_character_,
      score = pident * alignment_length / query_length,
      position = paste0(subject_start, "-", subject_end, "(", strand, ")")
    ) %>%
    
    # Keep hits that pass identity and score thresholds
    dplyr::filter(
      pident >= min_pident,
      score >= min_score
    ) %>%
    
    # Keep the final BLAST annotation columns
    dplyr::select(
      subject_id, query_id, position, subject_start, subject_end, strand, score,
      pident, alignment_length, mismatch, gapopen, evalue, bitscore
    )
  
  hits
}


findFrankingSite <- function(
    df,
    franking_table,
    sequence_col = "sequence",
    id_col = "contig_id",
    ...
) {
  # Search franking sequences
  hits <- sequenceBlast(
    df = df,
    query_table = franking_table,
    sequence_col = sequence_col,
    id_col = id_col,
    ...
  )
  
  # Keep the original table when no franking site is found
  if (nrow(hits) == 0) {
    attr(df, "blast_hits") <- hits
    return(df)
  }
  
  # Summarize all positions found for each query
  summary <- hits %>%
    dplyr::arrange(subject_id, query_id, subject_start) %>%
    dplyr::summarise(
      count = dplyr::n(),
      position = paste(position, collapse = "|"),
      score = paste(round(score, 2), collapse = "|"),
      .by = c(subject_id, query_id)
    ) %>%
    tidyr::pivot_wider(
      id_cols = subject_id,
      names_from = query_id,
      values_from = c(count, position, score),
      names_glue = "{query_id}_{.value}"
    )
  
  # Add franking site results to the original table
  result <- df %>%
    dplyr::left_join(
      summary,
      by = stats::setNames("subject_id", id_col)
    )
  
  attr(result, "blast_hits") <- hits
  
  result
}


findLinker <- function(
    df,
    linker_table,
    sequence_col = "insert_seq",
    id_col = "contig_id",
    ...
) {
  # Prepare linker sequences as BLAST queries
  query_table <- linker_table %>%
    dplyr::transmute(
      query_id = linker_name,
      sequence
    )
  
  # Search linker sequences
  hits <- sequenceBlast(
    df = df,
    query_table = query_table,
    sequence_col = sequence_col,
    id_col = id_col,
    ...
  )
  
  # Keep the best linker hit for each subject
  best_linker <- hits %>%
    dplyr::arrange(
      subject_id,
      dplyr::desc(score),
      dplyr::desc(bitscore),
      evalue,
      subject_start
    ) %>%
    dplyr::distinct(subject_id, .keep_all = TRUE) %>%
    dplyr::transmute(
      subject_id,
      linker_name = query_id,
      linker_position = position,
      linker_score = round(score, 2)
    )
  
  # Add linker annotation to the original table
  result <- df %>%
    dplyr::left_join(
      best_linker,
      by = stats::setNames("subject_id", id_col)
    ) %>%
    dplyr::select(
      dplyr::all_of(colnames(df)),
      linker_name,
      linker_position,
      linker_score
    )
  
  attr(result, "blast_hits") <- hits
  
  result
}


findKnownSequence <- function(
    df,
    known_table,
    sequence_col = "insert_seq",
    id_col = "contig_id",
    name_col = "CAR_name",
    ...
) {
  # Prepare known sequences as BLAST queries
  query_table <- known_table %>%
    dplyr::transmute(
      query_id = .data[[name_col]],
      sequence
    )
  
  # Search known sequences
  hits <- sequenceBlast(
    df = df,
    query_table = query_table,
    sequence_col = sequence_col,
    id_col = id_col,
    ...
  )
  
  # Keep the best known-sequence hit for each subject
  best_hit <- hits %>%
    dplyr::arrange(
      subject_id,
      dplyr::desc(score),
      dplyr::desc(bitscore),
      evalue,
      subject_start
    ) %>%
    dplyr::distinct(subject_id, .keep_all = TRUE) %>%
    dplyr::transmute(
      subject_id,
      scFv_type = query_id,
      known_sequence_position = position,
      known_sequence_score = round(score, 2)
    )
  
  # Add known-sequence annotation and mark unmatched subjects as Novel
  result <- df %>%
    dplyr::left_join(
      best_hit,
      by = stats::setNames("subject_id", id_col)
    ) %>%
    dplyr::mutate(
      scFv_type = dplyr::coalesce(scFv_type, "Novel")
    ) %>%
    dplyr::select(
      dplyr::all_of(colnames(df)),
      scFv_type,
      known_sequence_position,
      known_sequence_score
    )
  
  attr(result, "blast_hits") <- hits
  
  result
}


# Cut sequences ----
pick_hit <- function(position, score, criteria, min_score = NULL) {
  positions <- strsplit(position, "\\|")[[1]]
  scores <- as.numeric(strsplit(score, "\\|")[[1]])
  
  start <- as.integer(stringr::str_extract(positions, "^[0-9]+"))
  end <- as.integer(stringr::str_extract(positions, "(?<=-)[0-9]+"))
  
  if (!is.null(min_score)) {
    keep <- scores >= min_score
    positions <- positions[keep]
    scores <- scores[keep]
    start <- start[keep]
    end <- end[keep]
  }
  
  if (length(scores) == 0) return(c(start = NA_real_, end = NA_real_, score = NA_real_))
  
  order_list <- lapply(criteria, function(x) {
    if (x == "score") return(-scores)
    if (x == "leftmost") return(start)
    if (x == "rightmost") return(-end)
    stop("Unknown criteria: ", x)
  })
  
  index <- do.call(order, order_list)[1]
  c(start = start[index], end = end[index], score = scores[index])
}


cutFrankingSite <- function(
    df,
    sequence_col = "sequence",
    quality_col = "sequence_quals",
    left_target = "Myc_tag",
    right_target = "3_tail",
    left_criteria = c("score", "leftmost"),
    right_criteria = c("score", "rightmost"),
    flank_length = 5,
    min_score = NULL
) {
  left_count_col <- paste0(left_target, "_count")
  left_position_col <- paste0(left_target, "_position")
  left_score_col <- paste0(left_target, "_score")
  right_count_col <- paste0(right_target, "_count")
  right_position_col <- paste0(right_target, "_position")
  right_score_col <- paste0(right_target, "_score")
  
  # Check required columns
  required_columns <- c(
    sequence_col, quality_col,
    left_count_col, left_position_col, left_score_col,
    right_count_col, right_position_col, right_score_col
  )
  missing_columns <- setdiff(required_columns, colnames(df))
  if (length(missing_columns) > 0) stop("Missing columns: ", paste(missing_columns, collapse = ", "), call. = FALSE)
  
  # Check if sequence and quality strings have equal lengths
  invalid_quality <- df %>%
    dplyr::filter(
      is.na(.data[[sequence_col]]) |
        is.na(.data[[quality_col]]) |
        nchar(.data[[sequence_col]]) != nchar(.data[[quality_col]])
    ) %>%
    dplyr::transmute(
      contig_id,
      sequence_length = nchar(.data[[sequence_col]]),
      quality_length = nchar(.data[[quality_col]])
    )
  
  if (nrow(invalid_quality) > 0) {
    print(tibble::as_tibble(invalid_quality), n = Inf, width = Inf)
    stop("Sequence and quality lengths do not match.", call. = FALSE)
  }
  
  # Keep contigs containing both flanking sites
  result <- df %>%
    dplyr::filter(
      .data[[left_count_col]] > 0,
      .data[[right_count_col]] > 0
    )
  
  # Select the left and right flanking hits
  left_hit <- purrr::map2_dfr(
    result[[left_position_col]],
    result[[left_score_col]],
    ~ as.list(pick_hit(.x, .y, left_criteria, min_score))
  ) %>%
    dplyr::rename(
      left_hit_start = start,
      left_hit_end = end,
      left_hit_score = score
    )
  
  right_hit <- purrr::map2_dfr(
    result[[right_position_col]],
    result[[right_score_col]],
    ~ as.list(pick_hit(.x, .y, right_criteria, min_score))
  ) %>%
    dplyr::rename(
      right_hit_start = start,
      right_hit_end = end,
      right_hit_score = score
    )
  
  # Cut the insert sequence and its matching quality string
  result %>%
    dplyr::bind_cols(left_hit, right_hit) %>%
    dplyr::mutate(
      sequence_length = nchar(.data[[sequence_col]]),
      insert_start = left_hit_end + 1L,
      insert_end = right_hit_start - 1L,
      valid_franking_order = insert_start <= insert_end,
      
      left_franking_seq = dplyr::if_else(
        valid_franking_order,
        stringr::str_sub(
          .data[[sequence_col]],
          pmax(1L, left_hit_start - flank_length),
          left_hit_end
        ),
        NA_character_
      ),
      insert_seq = dplyr::if_else(
        valid_franking_order,
        stringr::str_sub(.data[[sequence_col]], insert_start, insert_end),
        NA_character_
      ),
      insert_seq_quals = dplyr::if_else(
        valid_franking_order,
        stringr::str_sub(.data[[quality_col]], insert_start, insert_end),
        NA_character_
      ),
      right_franking_seq = dplyr::if_else(
        valid_franking_order,
        stringr::str_sub(
          .data[[sequence_col]],
          right_hit_start,
          pmin(sequence_length, right_hit_end + flank_length)
        ),
        NA_character_
      ),
      insert_seq_length = nchar(insert_seq)
    ) %>%
    dplyr::select(
      dplyr::all_of(sequence_col),
      left_franking_seq,
      insert_seq,
      insert_seq_quals,
      right_franking_seq,
      dplyr::everything(),
      -dplyr::all_of(quality_col),
      -insert_start,
      -insert_end
    )
}


cutLinker <- function(
    df,
    sequence_col = "insert_seq",
    linker_position_col = "linker_position"
) {
  if (!sequence_col %in% colnames(df)) stop("Cannot find sequence_col: ", sequence_col)
  if (!linker_position_col %in% colnames(df)) stop("Cannot find linker_position_col: ", linker_position_col)
  
  df %>%
    dplyr::mutate(
      linker_start = as.integer(
        stringr::str_extract(.data[[linker_position_col]], "^[0-9]+")
      ),
      linker_end = as.integer(
        stringr::str_extract(.data[[linker_position_col]], "(?<=-)[0-9]+")
      ),
      
      V1_seq = stringr::str_sub(
        .data[[sequence_col]],
        1,
        linker_start - 1
      ),
      
      linker_seq = stringr::str_sub(
        .data[[sequence_col]],
        linker_start,
        linker_end
      ),
      
      V2_seq = stringr::str_sub(
        .data[[sequence_col]],
        linker_end + 1,
        nchar(.data[[sequence_col]])
      )
    ) %>%
    dplyr::select(
      sequence,
      left_franking_seq,
      insert_seq,
      right_franking_seq,
      V1_seq,
      linker_seq,
      V2_seq,
      dplyr::everything()
    )
}


# IgBLAST ----
# Prepare V1/V2 queries and write FASTA
prepare_scfv_igblast_queries <- function(scFvs, ids, work_dir) {
  # Assemble query map
  query_map <- scFvs %>% 
    dplyr::select(dplyr::all_of(ids), V1_seq, V2_seq) %>%
    
    # Split V1 and V2 sequences
    tidyr::pivot_longer(
      cols = c(V1_seq, V2_seq),
      names_to = "slot",
      values_to = "sequence"
    ) %>%
    dplyr::mutate(
      slot = dplyr::case_when(
        slot == "V1_seq" ~ "V1",
        slot == "V2_seq" ~ "V2",
        TRUE ~ slot
      )
    ) %>%
    dplyr::filter(!is.na(sequence), sequence != "") %>%
    
    # Assign one query ID to each unique sequence
    dplyr::mutate(
      query_id = sprintf(
        "q%06d",
        match(sequence, unique(sequence))
      )
    )
  
  # Write unique query sequences to FASTA
  query_fasta <- file.path(work_dir, "scfv_queries.fasta")
  query_map %>%
    dplyr::distinct(query_id, sequence) %>%
    dplyr::transmute(
      fasta_record = paste0(
        ">", query_id, "\n",
        sequence
      )
    ) %>%
    dplyr::pull(fasta_record) %>%
    writeLines(query_fasta)
  
  list(
    query_map = query_map,
    query_fasta = query_fasta
  )
}


# Run IgBLAST and read AIRR output
run_igblast_airr <- function(
    query_fasta,
    organism,
    output_file,
    igblast_share,
    germline_db_dir,
    threads = 4
) {
  # Assemble IgBLAST arguments
  igblast_args <- c(
    "-query", shQuote(query_fasta),
    "-organism", organism,
    "-ig_seqtype", "Ig",
    "-germline_db_V", shQuote(file.path(germline_db_dir, paste0("imgt_", organism, "_ig_v"))),
    "-germline_db_D", shQuote(file.path(germline_db_dir, paste0("imgt_", organism, "_ig_d"))),
    "-germline_db_J", shQuote(file.path(germline_db_dir, paste0("imgt_", organism, "_ig_j"))),
    "-auxiliary_data", shQuote(file.path(igblast_share, "optional_file", paste0(organism, "_gl.aux"))),
    "-domain_system", "imgt",
    "-outfmt", "19",
    "-out", shQuote(output_file),
    "-num_threads", as.character(threads)
  )
  
  # Run IgBLAST
  unlink(output_file)
  message("Running ", organism, " IgBLAST...")
  
  run_log <- system2("igblastn", args = igblast_args, stdout = TRUE, stderr = TRUE)
  run_status <- attr(run_log, "status")
  
  if (!is.null(run_status) && run_status != 0) {
    stop(
      stringr::str_to_title(organism),
      " IgBLAST failed:\n",
      paste(run_log, collapse = "\n")
    )
  }
  
  # Read AIRR result
  readr::read_tsv(
    output_file,
    show_col_types = FALSE,
    na = c("", "NA", "N/A")
  ) %>%
    dplyr::rename(query_id = sequence_id) %>%
    dplyr::mutate(host_candidate = stringr::str_to_title(organism))
}


# Select host, widen V1/V2 annotations and assemble final table
format_scfv_igblast_output <- function(
    scFvs,
    query_annotation,
    ids,
    work_dir,
    min_host_score_difference = 5
) {
  # Add missing optional AIRR columns as NA
  expected_airr_columns <- c(
    "locus", "v_call", "d_call", "j_call",
    "productive", "stop_codon", "vj_in_frame", "v_frameshift", "j_frameshift",
    "fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4",
    "fwr1_aa", "cdr1_aa", "fwr2_aa", "cdr2_aa",
    "fwr3_aa", "cdr3_aa", "fwr4_aa"
  )
  
  missing_airr_columns <- setdiff(expected_airr_columns, colnames(query_annotation))
  query_annotation[missing_airr_columns] <- NA
  
  # Select the best host for each contig
  host_selection <- query_annotation %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(ids)), host_candidate) %>%
    dplyr::summarise(host_score = sum(hit_score, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::across(dplyr::all_of(ids)), dplyr::desc(host_score), host_candidate) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(ids))) %>%
    dplyr::summarise(
      host_for_annotation = dplyr::first(host_candidate),
      best_host_score = dplyr::first(host_score),
      second_host_score = dplyr::nth(host_score, 2, default = NA_real_),
      host = dplyr::case_when(
        is.na(second_host_score) ~ host_for_annotation,
        best_host_score - second_host_score < min_host_score_difference ~ "Ambiguous_Human_Mouse",
        TRUE ~ host_for_annotation
      ),
      .groups = "drop"
    )
  
  # Keep selected AIRR annotations and widen V1/V2
  slot_annotation <- query_annotation %>%
    dplyr::inner_join(
      host_selection %>%
        dplyr::select(dplyr::all_of(ids), host_for_annotation, host),
      by = ids,
      relationship = "many-to-one"
    ) %>%
    dplyr::filter(host_candidate == host_for_annotation) %>%
    dplyr::transmute(
      dplyr::across(dplyr::all_of(ids)), host, slot, locus,
      V_call = v_call, D_call = d_call, J_call = j_call,
      productive, stop_codon, vj_in_frame, v_frameshift, j_frameshift,
      FR1 = fwr1, CDR1 = cdr1, FR2 = fwr2, CDR2 = cdr2,
      FR3 = fwr3, CDR3 = cdr3, FR4 = fwr4,
      FR1_aa = fwr1_aa, CDR1_aa = cdr1_aa, FR2_aa = fwr2_aa, CDR2_aa = cdr2_aa,
      FR3_aa = fwr3_aa, CDR3_aa = cdr3_aa, FR4_aa = fwr4_aa
    ) %>%
    tidyr::pivot_wider(
      id_cols = dplyr::all_of(c(ids, "host")),
      names_from = slot,
      values_from = c(
        locus, V_call, D_call, J_call,
        productive, stop_codon, vj_in_frame, v_frameshift, j_frameshift,
        FR1, CDR1, FR2, CDR2, FR3, CDR3, FR4,
        FR1_aa, CDR1_aa, FR2_aa, CDR2_aa, FR3_aa, CDR3_aa, FR4_aa
      ),
      names_glue = "{slot}_{.value}"
    )
  
  # Join IgBLAST annotations back to source rows
  final_table <- scFvs %>%
    dplyr::left_join(
      slot_annotation,
      by = ids,
      relationship = "one-to-one"
    ) %>%
    dplyr::mutate(
      chain_config = dplyr::case_when(
        V1_locus == "IGH" & V2_locus == "IGK" ~ "HL_kappa",
        V1_locus == "IGK" & V2_locus == "IGH" ~ "LH_kappa",
        V1_locus == "IGH" & V2_locus == "IGL" ~ "HL_lambda",
        V1_locus == "IGL" & V2_locus == "IGH" ~ "LH_lambda",
        TRUE ~ "other"
      )
    ) %>%
    dplyr::select(
      dplyr::all_of(ids), barcode,
      insert_seq_length,
      scFv_type, host, chain_config, linker_name,
      insert_seq, insert_seq_quals, V1_seq, linker_seq, V2_seq,
      is_cell, high_confidence, reads, umis,
      
      V1_V_call, V1_D_call, V1_J_call,
      V1_productive, V1_stop_codon, V1_vj_in_frame, V1_v_frameshift, V1_j_frameshift,
      V1_FR1, V1_CDR1, V1_FR2, V1_CDR2, V1_FR3, V1_CDR3, V1_FR4,
      V1_FR1_aa, V1_CDR1_aa, V1_FR2_aa, V1_CDR2_aa,
      V1_FR3_aa, V1_CDR3_aa, V1_FR4_aa,
      
      V2_V_call, V2_D_call, V2_J_call,
      V2_productive, V2_stop_codon, V2_vj_in_frame, V2_v_frameshift, V2_j_frameshift,
      V2_FR1, V2_CDR1, V2_FR2, V2_CDR2, V2_FR3, V2_CDR3, V2_FR4,
      V2_FR1_aa, V2_CDR1_aa, V2_FR2_aa, V2_CDR2_aa,
      V2_FR3_aa, V2_CDR3_aa, V2_FR4_aa
    ) %>%
    readr::write_tsv(
      file.path(work_dir, "scfv_annotation_ig.tsv"),
      na = ""
    )
  
  final_table
}


annotate_scfv_igblast <- function(
    scFvs,
    ids = c("sample_id", "contig_id"),
    work_dir = NULL, # Save working data into a folder
    threads = 4,
    min_host_score_difference = 5
) {
  # Check required columns
  missing_columns <- setdiff(c(ids, "V1_seq", "V2_seq"), colnames(scFvs))
  if (length(missing_columns) > 0) stop("Missing columns: ", paste(missing_columns, collapse = ", "))
  
  # Check if IDs uniquely identify each source row
  duplicated_ids <- scFvs %>%
    dplyr::count(dplyr::across(dplyr::all_of(ids)), name = "n") %>%
    dplyr::filter(n > 1)
  
  if (nrow(duplicated_ids) > 0) {
    print(tibble::as_tibble(duplicated_ids), n = Inf, width = Inf)
    stop("These IDs do not uniquely identify each row.", call. = FALSE)
  }
  
  # Use tempfile if work_dir is not provided
  if (is.null(work_dir)) {
    work_dir <- tempfile("scfv_igblast_")
    on.exit(unlink(work_dir, recursive = TRUE, force = TRUE), add = TRUE)
  }
  dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Read paths configured by the Conda environment
  igblast_share <- Sys.getenv("IGBLAST_SHARE")
  germline_db_dir <- Sys.getenv("GERMLINE_DB_DIR")
  
  # Prepare query data
  query_data <- prepare_scfv_igblast_queries(
    scFvs = scFvs,
    ids = ids,
    work_dir = work_dir
  )
  
  # Run mouse and human IgBLAST
  mouse_airr <- run_igblast_airr(
    query_fasta = query_data$query_fasta,
    organism = "mouse",
    output_file = file.path(work_dir, "scfv_mouse.airr.tsv"),
    igblast_share = igblast_share,
    germline_db_dir = germline_db_dir,
    threads = threads
  )
  
  human_airr <- run_igblast_airr(
    query_fasta = query_data$query_fasta,
    organism = "human",
    output_file = file.path(work_dir, "scfv_human.airr.tsv"),
    igblast_share = igblast_share,
    germline_db_dir = germline_db_dir,
    threads = threads
  )
  
  # Map AIRR results back to source contigs and V1/V2 slots
  query_annotation <- query_data$query_map %>%
    dplyr::left_join(
      dplyr::bind_rows(mouse_airr, human_airr),
      by = "query_id",
      relationship = "many-to-many"
    ) %>%
    dplyr::mutate(
      hit_score = rowSums(
        dplyr::across(
          dplyr::any_of(c("v_score", "d_score", "j_score")),
          ~ suppressWarnings(as.numeric(.x))
        ),
        na.rm = TRUE
      )
    ) %>%
    readr::write_tsv(
      file.path(work_dir, "scfv_query_annotations.airr.tsv"),
      na = ""
    )
  
  # Assemble final output table
  format_scfv_igblast_output(
    scFvs = scFvs,
    query_annotation = query_annotation,
    ids = ids,
    work_dir = work_dir,
    min_host_score_difference = min_host_score_difference
  )
}


# Post IgBLAST cleaning ----
translate_nt <- function(sequence) {
  result <- rep(NA_character_, length(sequence))
  valid <- !is.na(sequence) & sequence != ""
  
  if (any(valid)) {
    result[valid] <- suppressWarnings(
      as.character(
        Biostrings::translate(
          Biostrings::DNAStringSet(toupper(sequence[valid])),
          if.fuzzy.codon = "X"
        )
      )
    )
  }
  
  result
}

igblast_qc <- function(
    scFvs,
    keep_qc = FALSE,
    keep_fail = FALSE
) {
  # Set IgBLAST amino-acid annotation columns
  igblast_aa_columns <- c(
    "V1_FR1_aa", "V1_CDR1_aa", "V1_FR2_aa", "V1_CDR2_aa",
    "V1_FR3_aa", "V1_CDR3_aa", "V1_FR4_aa",
    "V2_FR1_aa", "V2_CDR1_aa", "V2_FR2_aa", "V2_CDR2_aa",
    "V2_FR3_aa", "V2_CDR3_aa", "V2_FR4_aa"
  )
  
  # Check required columns
  required_columns <- c(igblast_aa_columns, "insert_seq", "insert_seq_length","V1_productive","V2_productive")
  missing_columns <- setdiff(required_columns, colnames(scFvs))
  if (length(missing_columns) > 0) stop("Missing: ", paste(missing_columns, collapse = ", "), call. = FALSE)
  
  # Apply internal QC
  result <- scFvs %>%
    dplyr::mutate(
      # Reject only sequences explicitly classified as nonproductive
      productive_fail = (V1_productive %in% FALSE | V2_productive %in% FALSE),
      
      # Complete insert must contain full codons
      insert_frame_fail = (is.na(insert_seq_length) | insert_seq_length %% 3 != 0),
      
      # Translate all available insert sequences
      insert_seq_aa = translate_nt(insert_seq),
      V1_seq_aa = translate_nt(V1_seq),
      linker_seq_aa = translate_nt(linker_seq),
      V2_seq_aa = translate_nt(V2_seq),
      
      # Compare IgBLAST AA annotations with the translated insert
      igblast_aa_fail = dplyr::if_any(
        dplyr::all_of(igblast_aa_columns),
        ~ !is.na(.x) & .x != "" & !is.na(insert_seq_aa) &
          !stringr::str_detect(insert_seq_aa, stringr::fixed(toupper(.x)))
      )
    ) %>%
    dplyr::relocate(insert_seq_aa, .after = insert_seq_quals) %>%
    dplyr::relocate(V1_seq_aa, .after = V2_seq) %>%
    dplyr::relocate(linker_seq_aa, .after = V1_seq_aa) %>%
    dplyr::relocate(V2_seq_aa, .after = linker_seq_aa)
  
  # Remove failed sequences unless requested
  if (!keep_fail) {
    result <- result %>%
      dplyr::filter(!productive_fail, !insert_frame_fail, !igblast_aa_fail) %>%
      dplyr::select(
        -c(productive_fail, insert_frame_fail, igblast_aa_fail)
      )
  }
  
  # Remove QC unless requested
  if (!keep_qc) {
    result <- result %>%
      dplyr::select(
        -c(V1_productive, V1_stop_codon, V1_vj_in_frame, V1_v_frameshift, V1_j_frameshift,
           V2_productive, V2_stop_codon, V2_vj_in_frame, V2_v_frameshift, V2_j_frameshift)
      )
  }
  
  result
}


deduplicate_exact_scfv <- function(
    scFvs,
    min_reads = 100
) {
  # Check required columns
  required_columns <- c(
    "sample_id", "barcode", "contig_id", "insert_seq",
    "insert_seq_quals", "reads", "umis"
  )
  missing_columns <- setdiff(required_columns, colnames(scFvs))
  if (length(missing_columns) > 0) {
    stop("Missing columns: ", paste(missing_columns, collapse = ", "), call. = FALSE)
  }
  
  if (length(min_reads) != 1L || is.na(min_reads) || min_reads < 0) {
    stop("min_reads must be a single non-negative number.", call. = FALSE)
  }
  
  # Calculate quality summaries and apply the read-support gate
  result <- scFvs %>%
    dplyr::mutate(
      insert_q_scores = purrr::map(
        insert_seq_quals,
        ~ { utf8ToInt(.x) - 33L }
      ),
      insert_min_q = purrr::map_dbl(
        insert_q_scores,
        ~ { if (length(.x) > 0) min(.x) else NA_real_ }
      ),
      insert_mean_q = purrr::map_dbl(
        insert_q_scores,
        ~ { if (length(.x) > 0) mean(.x) else NA_real_ }
      ),
      reads_pass = !is.na(reads) & reads >= min_reads
    ) %>%
    
    # Keep the best exact sequence within each cell
    dplyr::arrange(
      sample_id,
      barcode,
      insert_seq,
      dplyr::desc(umis),
      dplyr::desc(reads_pass),
      dplyr::desc(insert_min_q),
      dplyr::desc(insert_mean_q),
      contig_id
    ) %>%
    dplyr::distinct(
      sample_id,
      barcode,
      insert_seq,
      .keep_all = TRUE
    ) %>%
    dplyr::select(-insert_q_scores, -reads_pass)
  
  result
}


scfv_error_correction <- function(
    scFvs,
    max_nt_mismatches = 1,
    min_reads = 100,
    min_major_mismatch_q = 30,
    max_minor_mismatch_q = 20
) {
  # Check required columns
  required_columns <- c(
    "sample_id", "barcode", "contig_id", "insert_seq",
    "insert_seq_quals", "reads", "umis"
  )
  missing_columns <- setdiff(required_columns, colnames(scFvs))
  if (length(missing_columns) > 0) {
    stop("Missing columns: ", paste(missing_columns, collapse = ", "), call. = FALSE)
  }
  
  # Check parameters
  if (length(max_nt_mismatches) != 1L || is.na(max_nt_mismatches) ||
      max_nt_mismatches < 1 || max_nt_mismatches %% 1 != 0) {
    stop("max_nt_mismatches must be a positive integer.", call. = FALSE)
  }
  if (length(min_reads) != 1L || is.na(min_reads) || min_reads < 0) {
    stop("min_reads must be a non-negative number.", call. = FALSE)
  }
  if (length(min_major_mismatch_q) != 1L || is.na(min_major_mismatch_q)) {
    stop("min_major_mismatch_q must be a single number.", call. = FALSE)
  }
  if (length(max_minor_mismatch_q) != 1L || is.na(max_minor_mismatch_q)) {
    stop("max_minor_mismatch_q must be a single number.", call. = FALSE)
  }
  
  # Check sequence and quality lengths
  invalid_quality <- scFvs %>%
    dplyr::filter(
      is.na(insert_seq) |
        is.na(insert_seq_quals) |
        nchar(insert_seq) != nchar(insert_seq_quals)
    ) %>%
    dplyr::transmute(
      sample_id, barcode, contig_id,
      sequence_length = nchar(insert_seq),
      quality_length = nchar(insert_seq_quals)
    )
  
  if (nrow(invalid_quality) > 0) {
    print(tibble::as_tibble(invalid_quality), n = Inf, width = Inf)
    stop("Insert sequence and quality lengths do not match.", call. = FALSE)
  }
  
  # Prepare sequence quality and support metrics
  working_table <- scFvs %>%
    dplyr::mutate(
      .row_id = dplyr::row_number(),
      insert_q_scores = purrr::map(insert_seq_quals, ~ utf8ToInt(.x) - 33L),
      insert_min_q = purrr::map_dbl(
        insert_q_scores,
        ~ { if (length(.x) > 0) min(.x) else NA_real_ }
      ),
      insert_mean_q = purrr::map_dbl(
        insert_q_scores,
        ~ { if (length(.x) > 0) mean(.x) else NA_real_ }
      ),
      reads_pass = !is.na(reads) & reads >= min_reads
    ) %>%
    dplyr::arrange(
      sample_id,
      barcode,
      dplyr::desc(umis),
      dplyr::desc(reads_pass),
      dplyr::desc(insert_min_q),
      dplyr::desc(insert_mean_q),
      contig_id
    ) %>%
    dplyr::group_by(sample_id, barcode) %>%
    dplyr::mutate(.support_rank = dplyr::row_number()) %>%
    dplyr::ungroup()
  
  # Compare all sequence pairs within one cell
  compare_group <- function(group) {
    if (nrow(group) < 2) return(tibble::tibble())
    
    pair_indices <- utils::combn(seq_len(nrow(group)), 2, simplify = FALSE)
    
    purrr::map_dfr(pair_indices, function(index) {
      row_1 <- group[index[1], ]
      row_2 <- group[index[2], ]
      
      # Only compare equal-length sequences
      if (nchar(row_1$insert_seq) != nchar(row_2$insert_seq)) {
        return(tibble::tibble())
      }
      
      bases_1 <- strsplit(row_1$insert_seq, "", fixed = TRUE)[[1]]
      bases_2 <- strsplit(row_2$insert_seq, "", fixed = TRUE)[[1]]
      mismatch_positions <- which(bases_1 != bases_2)
      nt_mismatches <- length(mismatch_positions)
      
      if (nt_mismatches == 0 || nt_mismatches > max_nt_mismatches) {
        return(tibble::tibble())
      }
      
      # Select provisional major by UMI, read gate and overall quality
      if (row_1$.support_rank < row_2$.support_rank) {
        major <- row_1
        minor <- row_2
      } else {
        major <- row_2
        minor <- row_1
      }
      
      major_umis <- dplyr::coalesce(major$umis, -Inf)
      minor_umis <- dplyr::coalesce(minor$umis, -Inf)
      
      # Require equal or higher UMI support and a passed read gate
      support_pass <- major_umis >= minor_umis && major$reads_pass
      
      major_mismatch_q <- major$insert_q_scores[[1]][mismatch_positions]
      minor_mismatch_q <- minor$insert_q_scores[[1]][mismatch_positions]
      
      # Require high-quality major bases and low-quality minor bases
      quality_pass <- all(major_mismatch_q >= min_major_mismatch_q) &&
        all(minor_mismatch_q <= max_minor_mismatch_q)
      
      tibble::tibble(
        sample_id = major$sample_id,
        barcode = major$barcode,
        major_row_id = major$.row_id,
        minor_row_id = minor$.row_id,
        major_contig_id = major$contig_id,
        minor_contig_id = minor$contig_id,
        major_umis = major$umis,
        minor_umis = minor$umis,
        major_reads = major$reads,
        minor_reads = minor$reads,
        major_reads_pass = major$reads_pass,
        minor_reads_pass = minor$reads_pass,
        nt_mismatches = nt_mismatches,
        mismatch_positions = paste(mismatch_positions, collapse = ","),
        major_mismatch_q = paste(major_mismatch_q, collapse = ","),
        minor_mismatch_q = paste(minor_mismatch_q, collapse = ","),
        support_pass = support_pass,
        quality_pass = quality_pass,
        correction_decision = dplyr::if_else(
          support_pass & quality_pass,
          "Correct",
          "Retain"
        ),
        major_support_rank = major$.support_rank
      )
    })
  }
  
  # Build candidate correction report
  correction_report <- working_table %>%
    split(interaction(.$sample_id, .$barcode, drop = TRUE, lex.order = TRUE)) %>%
    purrr::map_dfr(compare_group)
  
  # Keep one best target when a minor matches multiple major sequences
  if (nrow(correction_report) == 0L) {
    applied_corrections <- tibble::tibble(
      minor_row_id = integer()
    )
  } else {
    applied_corrections <- correction_report %>%
      dplyr::filter(correction_decision == "Correct") %>%
      dplyr::arrange(
        minor_row_id,
        nt_mismatches,
        major_support_rank,
        major_contig_id
      ) %>%
      dplyr::distinct(minor_row_id, .keep_all = TRUE)
  }
  
  # Remove likely sequencing or assembly errors
  result <- working_table %>%
    dplyr::filter(!.row_id %in% applied_corrections$minor_row_id) %>%
    dplyr::arrange(.row_id) %>%
    dplyr::select(
      -insert_q_scores,
      -reads_pass,
      -.row_id,
      -.support_rank
    )
  
  attr(result, "correction_report") <- correction_report
  attr(result, "applied_corrections") <- applied_corrections
  
  result
}


collapse_gene_calls <- function(values) {
  values <- values[!is.na(values) & trimws(as.character(values)) != ""]
  
  if (length(values) == 0) {
    return(NA_character_)
  }
  
  calls <- unlist(strsplit(as.character(values), ",", fixed = TRUE))
  calls <- sort(unique(trimws(calls)))
  
  paste(calls, collapse = ",")
}


assign_clone_id <- function(
    scFvs,
    clone_registry = NULL,
    info_columns = "scFv_type",
    prefix = "Clone",
    keep_key = FALSE
) {
  # Set FR/CDR amino-acid columns used to define exact clones
  clone_columns <- c(
    "V1_FR1_aa", "V1_CDR1_aa", "V1_FR2_aa", "V1_CDR2_aa",
    "V1_FR3_aa", "V1_CDR3_aa", "V1_FR4_aa",
    "V2_FR1_aa", "V2_CDR1_aa", "V2_FR2_aa", "V2_CDR2_aa",
    "V2_FR3_aa", "V2_CDR3_aa", "V2_FR4_aa"
  )
  
  # Check source columns
  required_columns <- c(clone_columns, info_columns)
  missing_columns <- setdiff(required_columns, colnames(scFvs))
  if (length(missing_columns) > 0) {
    stop(
      "Missing columns in scFvs: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  
  # Build clone keys from exact FR/CDR amino-acid combinations
  build_clone_key <- function(data) {
    normalized_sequences <- data %>%
      dplyr::transmute(
        dplyr::across(
          dplyr::all_of(clone_columns),
          ~ dplyr::if_else(
            is.na(.x) | trimws(.x) == "",
            "<MISSING>",
            toupper(trimws(.x))
          )
        )
      )
    
    do.call(
      paste,
      c(normalized_sequences, sep = "|")
    )
  }
  
  # Collapse one information column for each clone
  collapse_clone_info <- function(values, column_name, clone_key) {
    values <- unique(values[!is.na(values) & trimws(as.character(values)) != ""])
    
    if (length(values) > 1) {
      stop(
        "Conflicting ", column_name, " values for clone key: ", clone_key,
        "\nValues: ", paste(values, collapse = ", "),
        call. = FALSE
      )
    }
    
    if (length(values) == 0) {
      return(NA_character_)
    }
    
    as.character(values)
  }
  
  # Prepare current clone candidates
  current_table <- scFvs %>%
    dplyr::mutate(
      .row_id = dplyr::row_number(),
      clone_key = build_clone_key(scFvs)
    )
  
  # IgBLAST gene call columns can contain ambiguous assignments
  gene_call_columns <- c(
    "V1_V_call", "V1_D_call", "V1_J_call",
    "V2_V_call", "V2_D_call", "V2_J_call"
  )
  
  strict_info_columns <- setdiff(info_columns, gene_call_columns)
  gene_call_columns <- intersect(gene_call_columns, info_columns)
  
  # Merge ambiguous IgBLAST calls within the same clone
  collapse_gene_calls <- function(values) {
    values <- values[
      !is.na(values) &
        trimws(as.character(values)) != ""
    ]
    
    if (length(values) == 0) {
      return(NA_character_)
    }
    
    calls <- unlist(
      strsplit(as.character(values), ",", fixed = TRUE)
    )
    
    calls <- sort(unique(trimws(calls)))
    
    paste(calls, collapse = ",")
  }
  
  current_registry <- current_table %>%
    dplyr::group_by(clone_key) %>%
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(clone_columns),
        dplyr::first
      ),
      dplyr::across(
        dplyr::all_of(strict_info_columns),
        ~ collapse_clone_info(
          .x,
          dplyr::cur_column(),
          dplyr::first(clone_key)
        )
      ),
      dplyr::across(
        dplyr::all_of(gene_call_columns),
        collapse_gene_calls
      ),
      .groups = "drop"
    )
  
  # Prepare an empty registry when no known table is provided
  if (is.null(clone_registry)) {
    clone_registry <- current_registry[0, ] %>%
      dplyr::mutate(
        clone_id = character(),
        .before = 1
      )
  } else {
    registry_required <- c("clone_id", clone_columns, info_columns)
    missing_registry_columns <- setdiff(
      registry_required,
      colnames(clone_registry)
    )
    
    if (length(missing_registry_columns) > 0) {
      stop(
        "Missing columns in clone_registry: ",
        paste(missing_registry_columns, collapse = ", "),
        call. = FALSE
      )
    }
    
    clone_registry <- clone_registry %>%
      dplyr::mutate(
        clone_key = build_clone_key(clone_registry)
      ) %>%
      dplyr::select(
        clone_id,
        clone_key,
        dplyr::all_of(clone_columns),
        dplyr::all_of(info_columns)
      )
    
    # Check whether one clone key maps to multiple clone IDs
    duplicated_keys <- clone_registry %>%
      dplyr::distinct(clone_key, clone_id) %>%
      dplyr::count(clone_key, name = "n_clone_ids") %>%
      dplyr::filter(n_clone_ids > 1)
    
    if (nrow(duplicated_keys) > 0) {
      print(duplicated_keys, n = Inf, width = Inf)
      stop(
        "Some clone keys map to multiple clone IDs.",
        call. = FALSE
      )
    }
    
    # Check whether one clone ID maps to multiple clone keys
    duplicated_ids <- clone_registry %>%
      dplyr::distinct(clone_id, clone_key) %>%
      dplyr::count(clone_id, name = "n_clone_keys") %>%
      dplyr::filter(n_clone_keys > 1)
    
    if (nrow(duplicated_ids) > 0) {
      print(duplicated_ids, n = Inf, width = Inf)
      stop(
        "Some clone IDs map to multiple clone keys.",
        call. = FALSE
      )
    }
    
    clone_registry <- clone_registry %>%
      dplyr::distinct(clone_key, .keep_all = TRUE)
  }
  
  # Find clone keys not yet present in the registry
  new_clones <- current_registry %>%
    dplyr::anti_join(
      clone_registry %>% dplyr::select(clone_key),
      by = "clone_key"
    ) %>%
    dplyr::arrange(clone_key)
  
  # Continue numbering from the largest existing clone ID
  existing_numbers <- stringr::str_match(
    clone_registry$clone_id,
    paste0("^", stringr::str_escape(prefix), "_(\\d+)$")
  )[, 2]
  
  last_clone_number <- suppressWarnings(
    max(as.integer(existing_numbers), na.rm = TRUE)
  )
  
  if (!is.finite(last_clone_number)) {
    last_clone_number <- 0L
  }
  
  if (nrow(new_clones) > 0) {
    new_clones <- new_clones %>%
      dplyr::mutate(
        clone_id = paste0(
          prefix,
          "_",
          last_clone_number + dplyr::row_number()
        )
      )
  } else {
    new_clones <- clone_registry[0, ]
  }
  
  # Update known information without changing existing clone IDs
  registry_updated <- clone_registry %>%
    dplyr::full_join(
      current_registry,
      by = c("clone_key", clone_columns),
      suffix = c(".registry", ".current")
    )
  
  for (column_name in info_columns) {
    registry_column <- paste0(column_name, ".registry")
    current_column <- paste0(column_name, ".current")
    
    registry_updated[[column_name]] <- dplyr::coalesce(
      registry_updated[[registry_column]],
      registry_updated[[current_column]]
    )
  }
  
  registry_updated <- registry_updated %>%
    dplyr::select(
      clone_id,
      clone_key,
      dplyr::all_of(clone_columns),
      dplyr::all_of(info_columns)
    ) %>%
    dplyr::filter(!is.na(clone_id))
  
  # Append newly identified clones
  registry_updated <- dplyr::bind_rows(
    registry_updated,
    new_clones %>%
      dplyr::select(
        clone_id,
        clone_key,
        dplyr::all_of(clone_columns),
        dplyr::all_of(info_columns)
      )
  ) %>%
    dplyr::distinct(clone_key, .keep_all = TRUE) %>%
    dplyr::arrange(as.numeric(stringr::str_extract(clone_id, "\\d+$")))
  
  # Add stable clone IDs to the current dataset
  annotated <- current_table %>%
    dplyr::left_join(
      registry_updated %>%
        dplyr::select(clone_key, clone_id),
      by = "clone_key",
      relationship = "many-to-one"
    ) %>%
    dplyr::arrange(.row_id) %>%
    dplyr::select(-.row_id) %>%
    dplyr::relocate(clone_id, .after = scFv_type)
  
  if (!keep_key) {
    annotated <- annotated %>%
      dplyr::select(-clone_key)
  } else {
    annotated <- annotated %>%
      dplyr::relocate(clone_key, .after = clone_id)
  }
  
  list(
    annotated = annotated,
    registry = registry_updated
  )
}


classify_scfv_barcode <- function(scFvs) {
  required_columns <- c("sample_id", "barcode", "clone_id")
  missing_columns <- setdiff(required_columns, colnames(scFvs))
  
  if (length(missing_columns) > 0) {
    stop("Missing columns: ", paste(missing_columns, collapse = ", "), call. = FALSE)
  }
  
  scFvs %>%
    dplyr::group_by(sample_id, barcode) %>%
    dplyr::mutate(
      clone_count = dplyr::n_distinct(clone_id, na.rm = TRUE),
      class = dplyr::case_when(
        clone_count == 0 ~ "Undefined",
        clone_count == 1 ~ "Single",
        clone_count > 1 ~ "Multiple"
      )
    ) %>%
    dplyr::ungroup()
}


# Plots ----
seqLengthPlot <- function(
    df,
    length_col = "insert_seq_length",
    binwidth = 10,
    limits = NULL,
    colors = c("#F8766D", "black", "#de6158", "#7a342f"),
    density_alpha = 0.5,
    xlab = "Sequence Length",
    ylab = "Frequency"
) {
  if (!length_col %in% colnames(df)) stop("Cannot find length_col: ", length_col)
  if (length(colors) != 4) stop("colors must be length 4:\n", "c(hist_fill, hist_border, density_fill, density_line)")
  
  histogram_fill <- colors[1]
  histogram_color <- colors[2]
  density_fill <- colors[3]
  density_color <- colors[4]
  
  plot_df <- df %>%
    dplyr::filter(!is.na(.data[[length_col]]))
  
  if (is.null(limits)) {
    limits <- c(
      0,
      ceiling(max(plot_df[[length_col]], na.rm = TRUE) / binwidth) *
        binwidth + binwidth
    )
  }
  
  ggplot2::ggplot(plot_df, ggplot2::aes(x = .data[[length_col]])) +
    ggplot2::geom_histogram(
      binwidth = binwidth,
      fill = histogram_fill,
      color = histogram_color
    ) +
    ggplot2::geom_density(
      ggplot2::aes(
        y = ggplot2::after_stat(density) * nrow(plot_df) * binwidth
      ),
      fill = density_fill,
      color = density_color,
      linewidth = 0.2,
      alpha = density_alpha
    ) +
    ggplot2::scale_x_continuous() +
    ggplot2::coord_cartesian(xlim = limits) +
    ggplot2::labs(
      x = xlab,
      y = ylab
    )
}


# Summary ----
repSummarize <- function(
    ...,
    annotation_col = "scFv_type"
) {
  data_list <- list(...)
  
  if (is.null(names(data_list)) || any(names(data_list) == "")) stop("Please provide named tables.")
  if (!"raw" %in% names(data_list)) stop("Please provide raw table.")
  if (!"annotation" %in% names(data_list)) stop("Please provide annotation table.")
  
  annotation_df <- data_list$annotation
  filter_list <- data_list[names(data_list) != "annotation"]
  
  if (!annotation_col %in% colnames(annotation_df)) stop("Cannot find annotation_col: ", annotation_col)
  
  filter_summary <- tibble::tibble(
    stage = names(filter_list),
    group = NA_character_,
    n = purrr::map_int(filter_list, nrow)
  ) %>%
    dplyr::mutate(
      parent_n = dplyr::lag(n),
      total_n = dplyr::first(n),
      percent_parent = n / parent_n,
      percent_total = n / total_n
    )
  
  annotation_parent_n <- dplyr::last(filter_summary$n)
  total_n <- dplyr::first(filter_summary$n)
  
  annotation_summary <- annotation_df %>%
    dplyr::count(.data[[annotation_col]], name = "n") %>%
    dplyr::rename_with(
      ~ "group",
      dplyr::all_of(annotation_col)
    ) %>%
    dplyr::mutate(
      stage = "annotation",
      parent_n = annotation_parent_n,
      total_n = total_n,
      percent_parent = n / parent_n,
      percent_total = n / total_n
    ) %>%
    dplyr::select(stage, group, n, parent_n, total_n, percent_parent, percent_total)
  
  filter_summary %>%
    dplyr::select(stage, group, n, percent_parent, percent_total) %>%
    dplyr::bind_rows(
      annotation_summary %>% dplyr::select(stage, group, n, percent_parent, percent_total)
    ) %>%
    dplyr::mutate(
      stage = dplyr::case_when(
        stage == "raw" ~ "Raw",
        stage == "annotation" ~ "Annotation",
        TRUE ~ paste0("Filter: ", stringr::str_remove(stage, "_filter$"))
      )
    ) %>%
    dplyr::rename(
      Stage = stage,
      Group = group,
      `%parent` = percent_parent,
      `%total` = percent_total
    ) %>%
    dplyr::mutate(
      `%parent` = round(`%parent` * 100, 4),
      `%total` = round(`%total` * 100, 4)
    )
}


# Build data sets for AI training ----
collectTrainingData <- function(
    scFv_annotation,
    seu_obj,
    metadata = NULL,
    features = NULL,
    assay = "RNA",
    layer = "data",
    barcode_col = "barcode",
    keep_cell_without_exp = FALSE
) {
  # Check barcode column
  if (!(barcode_col %in% colnames(scFv_annotation))) {
    stop("scFv_annotation is missing required column: ", barcode_col, call. = FALSE)
  }
  
  # Extract selected Seurat metadata
  seu_metadata <- seu_obj[[]] %>%
    tibble::rownames_to_column(barcode_col)
  
  missing_metadata <- setdiff(metadata, colnames(seu_metadata))
  if (length(missing_metadata) > 0) {
    stop("Missing metadata columns: ", paste(missing_metadata, collapse = ", "), call. = FALSE)
  }
  
  metadata_tbl <- seu_metadata %>%
    dplyr::select(dplyr::all_of(unique(c(barcode_col, metadata)))) %>%
    dplyr::distinct(.data[[barcode_col]], .keep_all = TRUE)
  
  # Extract selected feature expression
  if (is.null(features) || length(features) == 0) {
    expression_tbl <- tibble::tibble(
      !!barcode_col := colnames(seu_obj),
      .has_expression = TRUE
    )
  } else {
    available_features <- rownames(seu_obj[[assay]])
    missing_features <- setdiff(features, available_features)
    
    if (length(missing_features) > 0) {
      warning(
        "These features were not found and will be ignored: ",
        paste(missing_features, collapse = ", ")
      )
    }
    
    features <- intersect(features, available_features)
    
    expression_tbl <- SeuratObject::LayerData(
      object = seu_obj,
      assay = assay,
      layer = layer,
      features = features
    ) %>%
      t() %>%
      as.data.frame() %>%
      tibble::rownames_to_column(barcode_col) %>%
      dplyr::mutate(.has_expression = TRUE)
  }
  
  # Remove duplicated columns before joining
  metadata_cols <- setdiff(
    colnames(metadata_tbl),
    c(barcode_col, colnames(scFv_annotation))
  )
  
  feature_cols <- setdiff(
    colnames(expression_tbl),
    c(barcode_col, ".has_expression", colnames(scFv_annotation), metadata_cols)
  )
  
  metadata_tbl <- metadata_tbl %>%
    dplyr::select(dplyr::all_of(c(barcode_col, metadata_cols)))
  
  expression_tbl <- expression_tbl %>%
    dplyr::select(dplyr::all_of(c(barcode_col, feature_cols)), .has_expression)
  
  # Combine scFv annotation, metadata and expression
  training_data <- scFv_annotation %>%
    dplyr::left_join(metadata_tbl, by = barcode_col) %>%
    dplyr::left_join(expression_tbl, by = barcode_col)
  
  if (!keep_cell_without_exp) {
    training_data <- dplyr::filter(training_data, .has_expression %in% TRUE)
  }
  
  training_data %>%
    dplyr::select(
      dplyr::all_of(colnames(scFv_annotation)),
      dplyr::all_of(metadata_cols),
      dplyr::all_of(feature_cols)
    )
}
