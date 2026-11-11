# Reading Cellranger VDJ contig ----
readContig <- function(folder, barcode_index = NULL) {
  folder <- as.character(folder)
  
  if (length(folder) != 1) {
    stop("read_contig() only accepts one folder at a time.")
  }
  
  if (!dir.exists(folder)) {
    stop("Folder does not exist: ", folder)
  }
  
  # Check whether the input folder is a Cellranger sample folder or outs folder
  if (basename(normalizePath(folder)) == "outs") {
    outs_folder <- folder
    sample_id <- basename(dirname(normalizePath(folder)))
  } else {
    outs_folder <- file.path(folder, "outs")
    sample_id <- basename(normalizePath(folder))
  }
  
  if (!dir.exists(outs_folder)) {
    stop("Cannot find outs folder: ", outs_folder)
  }
  
  # Search key Cellranger VDJ output files
  fasta_file <- file.path(outs_folder, "all_contig.fasta")
  anno_file <- file.path(outs_folder, "all_contig_annotations.csv")
  
  if (!file.exists(fasta_file)) {
    stop("Cannot find all_contig.fasta: ", fasta_file)
  }
  
  if (!file.exists(anno_file)) {
    stop("Cannot find all_contig_annotations.csv: ", anno_file)
  }
  
  # Helper function: normalize contig ID and split barcode / contig number
  parse_contig_id <- function(contig_id, barcode_index = NULL) {
    parsed <- tibble::tibble(
      contig_id = contig_id
    ) %>%
      tidyr::separate(
        contig_id,
        into = c("barcode", "contig"),
        sep = "_contig_",
        remove = FALSE,
        convert = TRUE,
        fill = "right"
      )
    
    if (!is.null(barcode_index)) {
      barcode_index <- as.character(barcode_index)
      
      parsed <- parsed %>%
        dplyr::mutate(
          barcode = stringr::str_replace(
            barcode,
            "-[0-9]+$",
            paste0("-", barcode_index)
          ),
          contig_id = paste0(barcode, "_contig_", contig)
        )
    }
    
    parsed
  }
  
  # Input de-novo contig FASTA data
  fasta <- data.frame(
    V1 = readLines(fasta_file),
    stringsAsFactors = FALSE
  )
  
  # Create a new column to group FASTA records
  fasta <- fasta %>%
    dplyr::mutate(id = cumsum(grepl("^>", V1)))
  
  # Convert FASTA records into contig table
  contig <- fasta %>%
    dplyr::group_by(id) %>%
    dplyr::summarise(
      contig_id_raw = dplyr::first(V1),
      sequence = paste(V1[!grepl("^>", V1)], collapse = ""),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      contig_id_raw = gsub("^>", "", contig_id_raw),
      sample_id = sample_id
    )
  
  contig_parsed <- parse_contig_id(
    contig_id = contig$contig_id_raw,
    barcode_index = barcode_index
  )
  
  contig <- contig %>%
    dplyr::bind_cols(
      contig_parsed %>%
        dplyr::select(contig_id, barcode, contig)
    ) %>%
    dplyr::select(
      sequence,
      sample_id,
      contig_id,
      barcode,
      contig
    )
  
  # Input Cellranger contig annotation data
  contig_anno_raw <- read.csv(
    anno_file,
    stringsAsFactors = FALSE
  )
  
  required_cols <- c(
    "contig_id",
    "barcode",
    "is_cell",
    "high_confidence",
    "reads",
    "umis"
  )
  
  missing_cols <- setdiff(required_cols, colnames(contig_anno_raw))
  
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in all_contig_annotations.csv: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  contig_anno_parsed <- parse_contig_id(
    contig_id = contig_anno_raw$contig_id,
    barcode_index = barcode_index
  )
  
  contig_anno <- contig_anno_raw %>%
    dplyr::select(
      is_cell,
      high_confidence,
      reads,
      umis
    ) %>%
    dplyr::bind_cols(
      contig_anno_parsed %>%
        dplyr::select(contig_id, barcode, contig)
    ) %>%
    dplyr::mutate(sample_id = sample_id) %>%
    dplyr::select(
      sample_id,
      contig_id,
      barcode,
      contig,
      is_cell,
      high_confidence,
      reads,
      umis
    )
  
  # Merge contig sequence and annotation data
  contig <- contig %>%
    dplyr::left_join(
      contig_anno,
      by = c("sample_id", "contig_id", "barcode", "contig")
    ) %>%
    dplyr::select(
      sequence,
      sample_id,
      contig_id,
      barcode,
      contig,
      is_cell,
      high_confidence,
      reads,
      umis
    )
  
  contig
}


# BLAST function ----
sequenceBlast <- function(
    df,
    query_table,
    sequence_col = "sequence",
    id_col = "contig_id",
    min_score = 70,
    min_pident = 70,
    keep_reverse_strand = TRUE,
    blast_program = c("blastn", "blastp"),
    db_type = NULL,
    task = NULL,
    dust = NULL,
    max_target_seqs = 1000000,
    output_dir = NULL,
    remove_temp = is.null(output_dir)
) {
  
  # Select BLAST engine
  blast_program <- match.arg(blast_program)
  
  if (is.null(db_type)) {
    db_type <- ifelse(blast_program == "blastp", "prot", "nucl")
  }
  
  if (is.null(task)) {
    task <- ifelse(blast_program == "blastp", "blastp-short", "blastn-short")
  }
  
  # Check BLAST availability
  for (tool in c(blast_program, "makeblastdb")) {
    
    if (system2(tool, "-version", stdout = NULL, stderr = NULL) != 0) {
      stop(tool, " not found in PATH.")
    }
  }
  
  # Validate inputs
  if (!sequence_col %in% colnames(df)) {
    stop("Missing sequence column: ", sequence_col)
  }
  
  if (!id_col %in% colnames(df)) {
    stop("Missing id column: ", id_col)
  }
  
  if (!all(c("query_id", "sequence") %in% colnames(query_table))) {
    stop("query_table must contain query_id + sequence")
  }
  
  # Clean query IDs
  query_table <- query_table %>%
    dplyr::mutate(
      query_id = stringr::str_replace_all(query_id, "[^A-Za-z0-9_]", "_"),
      query_length = nchar(sequence)
    )
  
  if (anyDuplicated(query_table$query_id)) {
    stop("Duplicated query_id after cleaning")
  }
  
  # Internal subject ID mapping
  df_internal <- df %>%
    dplyr::mutate(
      .subject_id = .data[[id_col]],
      .blast_subject_id = paste0("subject_", dplyr::row_number())
    )
  
  subject_map <- df_internal %>%
    dplyr::select(.blast_subject_id, .subject_id)
  
  # Working directory
  work_dir <- if (is.null(output_dir)) {
    tempfile("sequenceBlast_")
  } else {
    output_dir
  }
  
  dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
  
  if (is.null(output_dir) && remove_temp) {
    on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)
  }
  
  subject_fasta <- file.path(work_dir, "subjects.fa")
  query_fasta <- file.path(work_dir, "queries.fa")
  blast_db <- file.path(work_dir, "blast_db")
  blast_out <- file.path(work_dir, "blast.tsv")
  
  # Write FASTA
  writeLines(
    purrr::map2_chr(
      df_internal$.blast_subject_id,
      df_internal[[sequence_col]],
      ~ paste0(">", .x, "\n", toupper(.y))
    ),
    subject_fasta
  )
  
  writeLines(
    purrr::map2_chr(
      query_table$query_id,
      query_table$sequence,
      ~ paste0(">", .x, "\n", toupper(.y))
    ),
    query_fasta
  )
  
  # Build BLAST DB
  makeblastdb_status <- system2(
    "makeblastdb",
    c("-in", subject_fasta,
      "-dbtype", db_type,
      "-out", blast_db)
  )
  
  if (makeblastdb_status != 0) {
    stop("makeblastdb failed")
  }
  
  # Output format
  outfmt <- paste(
    "6",
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "qlen"
  )
  
  blast_args <- c()
  
  if (!is.null(task)) {
    blast_args <- c(blast_args, "-task", task)
  }
  
  blast_args <- c(
    blast_args,
    "-query", query_fasta,
    "-db", blast_db
  )
  
  if (blast_program == "blastn") {
    strand_option <- ifelse(
      keep_reverse_strand,
      "both",
      "plus"
    )
    
    blast_args <- c(blast_args, "-strand", strand_option)
    
    if (!is.null(dust)) {
      blast_args <- c(blast_args, "-dust", dust)
    }
  }
  
  blast_args <- c(
    blast_args,
    "-max_target_seqs", format(max_target_seqs, scientific = FALSE),
    "-out", blast_out,
    "-outfmt", shQuote(outfmt)
  )
  
  blast_output <- system2(
    blast_program,
    blast_args,
    stdout = TRUE,
    stderr = TRUE
  )
  
  blast_status <- attr(blast_output, "status")
  
  if (!is.null(blast_status) && blast_status != 0) {
    stop(
      blast_program,
      " failed.\nCommand:\n",
      paste(blast_program, paste(blast_args, collapse = " ")),
      "\n\nBLAST output:\n",
      paste(blast_output, collapse = "\n"),
      call. = FALSE
    )
  }
  
  if (file.info(blast_out)$size == 0) {
    return(tibble::tibble())
  }
  
  hits <- read.delim(blast_out, header = FALSE)
  
  colnames(hits) <- c(
    "query_id",
    ".blast_subject_id",
    "pident",
    "alignment_length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "query_length"
  )
  
  hits %>%
    dplyr::left_join(subject_map, by = ".blast_subject_id") %>%
    dplyr::mutate(
      subject_start = pmin(sstart, send),
      subject_end = pmax(sstart, send),
      strand = ifelse(sstart <= send, "+", "-"),
      score = pident * alignment_length / query_length,
      position = paste0(subject_start, "-", subject_end, "(", strand, ")")
    ) %>%
    dplyr::filter(
      pident >= min_pident,
      score >= min_score
    ) %>%
    dplyr::select(
      subject_id = .subject_id,
      query_id,
      position,
      subject_start,
      subject_end,
      strand,
      score,
      pident,
      alignment_length,
      mismatch,
      gapopen,
      evalue,
      bitscore
    )
}


findFrankingSite <- function(
    df,
    franking_table,
    sequence_col = "sequence",
    id_col = "contig_id",
    ...
) {
  if (!all(c("query_id", "sequence") %in% colnames(franking_table))) {
    stop("franking_table must contain query_id + sequence")
  }
  
  hits <- sequenceBlast(
    df = df,
    query_table = franking_table,
    sequence_col = sequence_col,
    id_col = id_col,
    ...
  )
  
  if (nrow(hits) == 0) {
    attr(df, "blast_hits") <- hits
    return(df)
  }
  
  summary <- hits %>%
    dplyr::arrange(subject_id, query_id, subject_start) %>%
    dplyr::group_by(subject_id, query_id) %>%
    dplyr::summarise(
      count = dplyr::n(),
      position = paste(position, collapse = "|"),
      score = paste(round(score, 2), collapse = "|"),
      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(
      id_cols = subject_id,
      names_from = query_id,
      values_from = c(count, position, score),
      names_glue = "{query_id}_{.value}"
    )
  
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
  
  query_table <- linker_table %>%
    dplyr::transmute(
      query_id = linker_name,
      sequence = sequence
    )
  
  hits <- sequenceBlast(
    df = df,
    query_table = query_table,
    sequence_col = sequence_col,
    id_col = id_col,
    ...
  )
  
  if (nrow(hits) == 0) {
    return(df)
  }
  
  best_linker <- hits %>%
    dplyr::arrange(
      subject_id,
      dplyr::desc(score),
      dplyr::desc(bitscore),
      evalue,
      subject_start
    ) %>%
    dplyr::group_by(subject_id) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(
      subject_id,
      linker_name = query_id,
      linker_position = position,
      linker_score = round(score, 2)
    )
  
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
  required_cols <- c(name_col, "sequence")
  
  if (!all(required_cols %in% colnames(known_table))) {
    stop(
      "known_table must contain columns: ",
      paste(required_cols, collapse = ", ")
    )
  }
  
  query_table <- known_table %>%
    dplyr::transmute(
      query_id = .data[[name_col]],
      sequence = sequence
    )
  
  hits <- sequenceBlast(
    df = df,
    query_table = query_table,
    sequence_col = sequence_col,
    id_col = id_col,
    ...
  )
  
  if (nrow(hits) == 0) {
    result <- df %>%
      dplyr::mutate(
        sequence_name = "novel",
        known_sequence_position = NA_character_,
        known_sequence_score = NA_real_
      )
    
    attr(result, "blast_hits") <- hits
    return(result)
  }
  
  best_hit <- hits %>%
    dplyr::arrange(
      subject_id,
      dplyr::desc(score),
      dplyr::desc(bitscore),
      evalue,
      subject_start
    ) %>%
    dplyr::group_by(subject_id) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(
      subject_id,
      sequence_name = query_id,
      known_sequence_position = position,
      known_sequence_score = round(score, 2)
    )
  
  result <- df %>%
    dplyr::left_join(
      best_hit,
      by = stats::setNames("subject_id", id_col)
    ) %>%
    dplyr::mutate(
      sequence_name = dplyr::coalesce(
        sequence_name,
        "Novel"
      )
    ) %>%
    dplyr::select(
      dplyr::all_of(colnames(df)),
      sequence_name,
      known_sequence_position,
      known_sequence_score
    )
  
  attr(result, "blast_hits") <- hits
  
  result
}


# IgBLAST ----
normalize_nt <- function(sequence) {
  sequence <- toupper(sequence)
  sequence <- gsub("[^ACGTN]", "", sequence)
  sequence
}

wrap_fasta_sequence <- function(sequence, width = 80) {
  starts <- seq(1, nchar(sequence), by = width)
  pieces <- substring(sequence, starts, pmin(starts + width - 1, nchar(sequence)))
  paste(pieces, collapse = "\n")
}

write_fasta_table <- function(query_table, fasta_file) {
  fasta_lines <- purrr::map2_chr(
    query_table$query_id,
    query_table$sequence,
    function(query_id, sequence) {
      paste0(">", query_id, "\n", wrap_fasta_sequence(sequence))
    }
  )
  
  writeLines(fasta_lines, fasta_file)
}

is_blast_db_prefix <- function(prefix) {
  any(file.exists(paste0(prefix, c(".nhr", ".nin", ".nsq", ".ndb"))))
}

find_db_prefix <- function(igdata, candidates) {
  search_dirs <- unique(c(
    file.path(igdata, "database"),
    file.path(dirname(igdata), "database"),
    igdata
  ))
  
  for (search_dir in search_dirs) {
    for (candidate in candidates) {
      prefix <- file.path(search_dir, candidate)
      
      if (is_blast_db_prefix(prefix)) {
        return(prefix)
      }
    }
  }
  
  all_files <- list.files(
    igdata,
    recursive = TRUE,
    full.names = TRUE,
    pattern = "\\.(nhr|nin|nsq|ndb)$"
  )
  
  prefixes <- unique(sub("\\.(nhr|nin|nsq|ndb)$", "", all_files))
  
  for (candidate in candidates) {
    matched <- prefixes[basename(prefixes) == candidate]
    
    if (length(matched) > 0) {
      return(matched[[1]])
    }
  }
  
  NA_character_
}

find_optional_file <- function(igdata, filename) {
  candidates <- c(
    file.path(igdata, "optional_file", filename),
    file.path(dirname(igdata), "optional_file", filename)
  )
  
  existing <- candidates[file.exists(candidates)]
  
  if (length(existing) == 0) {
    return(NA_character_)
  }
  
  existing[[1]]
}

get_igblast_config <- function(organism) {
  igblastn_path <- Sys.which("igblastn")
  
  if (igblastn_path == "") {
    stop("Cannot find igblastn in PATH.")
  }
  
  igblast_env <- dirname(dirname(igblastn_path))
  igdata <- file.path(igblast_env, "share", "igblast")
  
  if (!dir.exists(igdata)) {
    stop("Cannot find IGDATA directory: ", igdata)
  }
  
  if (organism == "human") {
    v_candidates <- c("imgt_human_ig_v", "human_gl_V")
    d_candidates <- c("imgt_human_ig_d", "human_gl_D")
    j_candidates <- c("imgt_human_ig_j", "human_gl_J")
    c_candidates <- c("imgt_human_ig_c", "ncbi_human_c_genes", "human_gl_C")
    aux_file <- "human_gl.aux"
  } else if (organism == "mouse") {
    v_candidates <- c("imgt_mouse_ig_v", "mouse_gl_V")
    d_candidates <- c("imgt_mouse_ig_d", "mouse_gl_D")
    j_candidates <- c("imgt_mouse_ig_j", "mouse_gl_J")
    c_candidates <- c("imgt_mouse_ig_c", "ncbi_mouse_c_genes", "mouse_gl_C")
    aux_file <- "mouse_gl.aux"
  } else {
    stop("Unsupported organism: ", organism)
  }
  
  config <- list(
    igblastn = igblastn_path,
    igdata = igdata,
    organism = organism,
    db_v = find_db_prefix(igdata, v_candidates),
    db_d = find_db_prefix(igdata, d_candidates),
    db_j = find_db_prefix(igdata, j_candidates),
    db_c = find_db_prefix(igdata, c_candidates),
    auxiliary_data = find_optional_file(igdata, aux_file)
  )
  
  required_names <- c("db_v", "db_d", "db_j", "auxiliary_data")
  missing_required <- required_names[is.na(unlist(config[required_names]))]
  
  if (length(missing_required) > 0) {
    stop(
      "Missing IgBLAST resources for ",
      organism,
      ": ",
      paste(missing_required, collapse = ", ")
    )
  }
  
  config
}

make_query_tables <- function(scFvs) {
  long_table <- scFvs %>%
    mutate(row_id = dplyr::row_number()) %>%
    select(
      row_id,
      barcode,
      scFv_type,
      insert_V1_seq,
      insert_linker_seq,
      insert_V2_seq
    ) %>%
    pivot_longer(
      cols = c(insert_V1_seq, insert_V2_seq),
      names_to = "slot",
      values_to = "sequence"
    ) %>%
    mutate(
      slot = dplyr::case_when(
        slot == "insert_V1_seq" ~ "V1",
        slot == "insert_V2_seq" ~ "V2",
        TRUE ~ slot
      ),
      sequence = normalize_nt(sequence)
    ) %>%
    filter(!is.na(sequence), sequence != "")
  
  unique_table <- long_table %>%
    distinct(sequence) %>%
    mutate(query_id = sprintf("q%06d", dplyr::row_number()))
  
  long_table <- long_table %>%
    left_join(unique_table, by = "sequence")
  
  list(
    long = long_table,
    unique = unique_table
  )
}

run_igblast_airr <- function(
    fasta_file,
    output_file,
    config,
    threads = 4
) {
  Sys.setenv(IGDATA = config$igdata)
  
  args <- c(
    "-query", fasta_file,
    "-out", output_file,
    "-outfmt", "19",
    "-organism", config$organism,
    "-ig_seqtype", "Ig",
    "-domain_system", "imgt",
    "-germline_db_V", config$db_v,
    "-germline_db_D", config$db_d,
    "-germline_db_J", config$db_j,
    "-auxiliary_data", config$auxiliary_data,
    "-num_threads", as.character(threads)
  )
  
  if (!is.na(config$db_c)) {
    args <- c(args, "-c_region_db", config$db_c)
  }
  
  message("Running IgBLAST for ", config$organism)
  
  run_log <- system2(
    config$igblastn,
    args = args,
    stdout = TRUE,
    stderr = TRUE
  )
  
  status <- attr(run_log, "status")
  
  if (!is.null(status) && status != 0) {
    stop(
      "IgBLAST failed for ",
      config$organism,
      ":\n",
      paste(run_log, collapse = "\n")
    )
  }
  
  if (!file.exists(output_file) || file.size(output_file) == 0) {
    stop("IgBLAST output file is empty: ", output_file)
  }
  
  readr::read_tsv(
    output_file,
    col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE,
    na = c("", "NA", "N/A")
  ) %>%
    mutate(host_candidate = stringr::str_to_title(config$organism))
}

first_existing_column <- function(data, candidates) {
  existing <- candidates[candidates %in% colnames(data)]
  
  if (length(existing) == 0) {
    return(rep(NA, nrow(data)))
  }
  
  data[[existing[[1]]]]
}

to_numeric_safe <- function(x) {
  suppressWarnings(as.numeric(x))
}

format_position <- function(start, end) {
  start <- to_numeric_safe(start)
  end <- to_numeric_safe(end)
  
  ifelse(
    is.na(start) | is.na(end),
    NA_character_,
    paste0(start, "-", end)
  )
}

get_position_start <- function(position) {
  start <- stringr::str_extract(position, "^[0-9]+")
  suppressWarnings(as.numeric(start))
}

normalized_nt_width <- function(sequence, missing_as_zero = FALSE) {
  sequence <- normalize_nt(sequence)
  width <- nchar(sequence)
  
  if (missing_as_zero) {
    width[is.na(width)] <- 0L
  }
  
  width
}

format_frameshift <- function(start_position) {
  frameshift <- (start_position - 1) %% 3
  
  dplyr::case_when(
    is.na(frameshift) ~ NA_character_,
    frameshift == 0 ~ "0",
    frameshift == 1 ~ "+1(-2)",
    frameshift == 2 ~ "+2(-1)",
    TRUE ~ NA_character_
  )
}

parse_locus <- function(v_call, j_call = NA_character_, c_call = NA_character_) {
  combined <- dplyr::coalesce(v_call, j_call, c_call, "")
  
  dplyr::case_when(
    grepl("^IGH", combined) ~ "IGH",
    grepl("^IGK", combined) ~ "IGK",
    grepl("^IGL", combined) ~ "IGL",
    TRUE ~ NA_character_
  )
}

locus_to_chain_type <- function(locus) {
  dplyr::case_when(
    locus == "IGH" ~ "Heavy",
    locus == "IGK" ~ "Light_kappa",
    locus == "IGL" ~ "Light_lambda",
    TRUE ~ NA_character_
  )
}

score_igblast_hit <- function(data) {
  score_columns <- c(
    "v_score",
    "d_score",
    "j_score",
    "c_score",
    "v_identity",
    "j_identity"
  )
  
  existing <- score_columns[score_columns %in% colnames(data)]
  
  if (length(existing) == 0) {
    call_columns <- intersect(
      c("v_call", "d_call", "j_call", "c_call"),
      colnames(data)
    )
    
    if (length(call_columns) == 0) {
      return(rep(0, nrow(data)))
    }
    
    return(rowSums(!is.na(data[call_columns]), na.rm = TRUE))
  }
  
  score_matrix <- data %>%
    select(all_of(existing)) %>%
    mutate(across(everything(), to_numeric_safe)) %>%
    as.data.frame()
  
  rowSums(score_matrix, na.rm = TRUE)
}

get_segment_position <- function(data, segment) {
  segment <- tolower(segment)
  
  start <- first_existing_column(
    data,
    c(
      paste0(segment, "_sequence_start"),
      paste0(segment, "_alignment_start"),
      paste0(segment, "_seq_start"),
      paste0(segment, "_start")
    )
  )
  
  end <- first_existing_column(
    data,
    c(
      paste0(segment, "_sequence_end"),
      paste0(segment, "_alignment_end"),
      paste0(segment, "_seq_end"),
      paste0(segment, "_end")
    )
  )
  
  format_position(start, end)
}

clean_region_sequence <- function(sequence) {
  sequence <- toupper(sequence)
  sequence <- gsub("[-.[:space:]]", "", sequence)
  sequence[!is.na(sequence) & sequence == ""] <- NA_character_
  sequence
}

to_airr_logical <- function(x) {
  if (is.logical(x)) {
    return(x)
  }
  
  x <- stringr::str_to_lower(trimws(as.character(x)))
  
  dplyr::case_when(
    x %in% c("true", "t", "1", "yes", "y") ~ TRUE,
    x %in% c("false", "f", "0", "no", "n") ~ FALSE,
    TRUE ~ NA
  )
}

clean_aa_sequence <- function(sequence) {
  sequence <- toupper(sequence)
  sequence <- gsub("[-.[:space:]]", "", sequence)
  sequence[!is.na(sequence) & sequence == ""] <- NA_character_
  sequence
}

format_airr_frame_status <- function(v_frameshift, j_frameshift, vj_in_frame) {
  v_has_frameshift <- dplyr::coalesce(v_frameshift, FALSE)
  j_has_frameshift <- dplyr::coalesce(j_frameshift, FALSE)
  
  dplyr::case_when(
    v_has_frameshift | j_has_frameshift ~ "frameshift",
    !is.na(vj_in_frame) & !vj_in_frame ~ "out_of_frame",
    !is.na(vj_in_frame) & vj_in_frame ~ "0",
    !is.na(v_frameshift) & !v_frameshift &
      (is.na(j_frameshift) | !j_frameshift) ~ "0",
    TRUE ~ NA_character_
  )
}

pass_airr_domain_qc <- function(
    stop_codon,
    vj_in_frame,
    v_frameshift,
    j_frameshift
) {
  stop_ok <- !is.na(stop_codon) & !stop_codon
  frame_ok <- !is.na(vj_in_frame) & vj_in_frame
  v_shift_ok <- is.na(v_frameshift) | !v_frameshift
  j_shift_ok <- is.na(j_frameshift) | !j_frameshift
  
  stop_ok & frame_ok & v_shift_ok & j_shift_ok
}

find_ordered_region_positions <- function(full_sequence, regions) {
  full_sequence <- normalize_nt(full_sequence)
  current_start <- 1
  output <- rep(NA_character_, length(regions))
  names(output) <- names(regions)
  
  for (region_name in names(regions)) {
    region_sequence <- clean_region_sequence(regions[[region_name]])
    
    if (is.na(region_sequence) || nchar(region_sequence) == 0) {
      next
    }
    
    search_area <- substring(full_sequence, current_start)
    local_start <- regexpr(region_sequence, search_area, fixed = TRUE)[[1]]
    
    if (local_start < 0) {
      next
    }
    
    absolute_start <- current_start + local_start - 1
    absolute_end <- absolute_start + nchar(region_sequence) - 1
    
    output[[region_name]] <- paste0(absolute_start, "-", absolute_end)
    current_start <- absolute_end + 1
  }
  
  output
}

extract_airr_annotation <- function(airr_data) {
  data <- airr_data %>%
    rename(query_id = sequence_id) %>%
    mutate(
      hit_score = score_igblast_hit(pick(everything())),
      
      locus = dplyr::coalesce(
        first_existing_column(pick(everything()), "locus"),
        parse_locus(v_call, j_call, c_call)
      ),
      
      chain_type = locus_to_chain_type(locus),
      
      V = get_segment_position(pick(everything()), "v"),
      D = get_segment_position(pick(everything()), "d"),
      J = get_segment_position(pick(everything()), "j"),
      C = get_segment_position(pick(everything()), "c"),
      
      productive = to_airr_logical(
        first_existing_column(pick(everything()), "productive")
      ),
      stop_codon = to_airr_logical(
        first_existing_column(pick(everything()), "stop_codon")
      ),
      vj_in_frame = to_airr_logical(
        first_existing_column(pick(everything()), "vj_in_frame")
      ),
      v_frameshift = to_airr_logical(
        first_existing_column(pick(everything()), "v_frameshift")
      ),
      j_frameshift = to_airr_logical(
        first_existing_column(pick(everything()), "j_frameshift")
      ),
      
      sequence_alignment_aa = clean_aa_sequence(
        first_existing_column(
          pick(everything()),
          "sequence_alignment_aa"
        )
      ),
      v_sequence_alignment_aa = clean_aa_sequence(
        first_existing_column(
          pick(everything()),
          "v_sequence_alignment_aa"
        )
      ),
      j_sequence_alignment_aa = clean_aa_sequence(
        first_existing_column(
          pick(everything()),
          "j_sequence_alignment_aa"
        )
      ),
      
      FR1_seq = clean_region_sequence(
        first_existing_column(pick(everything()), c("fwr1", "fr1"))
      ),
      CDR1 = clean_region_sequence(
        first_existing_column(pick(everything()), "cdr1")
      ),
      FR2_seq = clean_region_sequence(
        first_existing_column(pick(everything()), c("fwr2", "fr2"))
      ),
      CDR2 = clean_region_sequence(
        first_existing_column(pick(everything()), "cdr2")
      ),
      FR3_seq = clean_region_sequence(
        first_existing_column(pick(everything()), c("fwr3", "fr3"))
      ),
      CDR3 = clean_region_sequence(
        first_existing_column(pick(everything()), c("cdr3", "junction"))
      ),
      FR4_seq = clean_region_sequence(
        first_existing_column(pick(everything()), c("fwr4", "fr4"))
      ),
      
      FR1_aa = clean_aa_sequence(
        first_existing_column(pick(everything()), c("fwr1_aa", "fr1_aa"))
      ),
      CDR1_aa = clean_aa_sequence(
        first_existing_column(pick(everything()), "cdr1_aa")
      ),
      FR2_aa = clean_aa_sequence(
        first_existing_column(pick(everything()), c("fwr2_aa", "fr2_aa"))
      ),
      CDR2_aa = clean_aa_sequence(
        first_existing_column(pick(everything()), "cdr2_aa")
      ),
      FR3_aa = clean_aa_sequence(
        first_existing_column(pick(everything()), c("fwr3_aa", "fr3_aa"))
      ),
      CDR3_aa = clean_aa_sequence(
        first_existing_column(pick(everything()), c("cdr3_aa", "junction_aa"))
      ),
      FR4_aa = clean_aa_sequence(
        first_existing_column(pick(everything()), c("fwr4_aa", "fr4_aa"))
      )
    )
  
  region_positions <- purrr::pmap_dfr(
    list(
      data$sequence,
      data$FR1_seq,
      data$CDR1,
      data$FR2_seq,
      data$CDR2,
      data$FR3_seq,
      data$CDR3,
      data$FR4_seq
    ),
    function(
    sequence,
    fr1,
    cdr1,
    fr2,
    cdr2,
    fr3,
    cdr3,
    fr4
    ) {
      positions <- find_ordered_region_positions(
        sequence,
        c(
          FR1 = fr1,
          CDR1_pos = cdr1,
          FR2 = fr2,
          CDR2_pos = cdr2,
          FR3 = fr3,
          CDR3_pos = cdr3,
          FR4 = fr4
        )
      )
      
      tibble(
        FR1 = positions[["FR1"]],
        FR2 = positions[["FR2"]],
        FR3 = positions[["FR3"]],
        FR4 = positions[["FR4"]]
      )
    }
  )
  
  bind_cols(data, region_positions) %>%
    select(
      query_id,
      host_candidate,
      hit_score,
      chain_type,
      locus,
      V,
      D,
      J,
      C,
      FR1,
      CDR1,
      FR2,
      CDR2,
      FR3,
      CDR3,
      FR4,
      FR1_aa,
      CDR1_aa,
      FR2_aa,
      CDR2_aa,
      FR3_aa,
      CDR3_aa,
      FR4_aa,
      v_call,
      d_call,
      j_call,
      c_call,
      productive,
      stop_codon,
      vj_in_frame,
      v_frameshift,
      j_frameshift,
      sequence_alignment_aa,
      v_sequence_alignment_aa,
      j_sequence_alignment_aa,
      any_of(c("v_identity", "j_identity"))
    )
}

pick_best_pair_host <- function(
    long_table,
    annotation_table,
    min_score_margin = 5
) {
  query_host_scores <- annotation_table %>%
    select(query_id, host_candidate, hit_score) %>%
    distinct(query_id, host_candidate, .keep_all = TRUE)
  
  pair_scores <- long_table %>%
    select(row_id, slot, query_id) %>%
    left_join(
      query_host_scores,
      by = "query_id",
      relationship = "many-to-many"
    ) %>%
    group_by(row_id, host_candidate) %>%
    summarise(
      pair_score = sum(hit_score, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(row_id, desc(pair_score))
  
  pair_scores %>%
    group_by(row_id) %>%
    summarise(
      host_for_annotation = first(host_candidate),
      best_score = first(pair_score),
      second_score = dplyr::nth(pair_score, 2, default = NA_real_),
      host = dplyr::case_when(
        is.na(second_score) ~ host_for_annotation,
        best_score - second_score < min_score_margin ~
          "Ambiguous_Human_Mouse",
        TRUE ~ host_for_annotation
      ),
      .groups = "drop"
    )
}

make_slot_wide <- function(slot_annotation, slot_name) {
  prefix <- paste0(slot_name, "_")
  
  slot_annotation %>%
    filter(slot == slot_name) %>%
    transmute(
      row_id,
      !!paste0(prefix, "chain_type") := chain_type,
      !!paste0(prefix, "locus") := locus,
      !!paste0(prefix, "hit_score") := hit_score,
      
      !!paste0(prefix, "V") := V,
      !!paste0(prefix, "D") := D,
      !!paste0(prefix, "J") := J,
      !!paste0(prefix, "C") := C,
      
      !!paste0(prefix, "FR1") := FR1,
      !!paste0(prefix, "CDR1") := CDR1,
      !!paste0(prefix, "FR2") := FR2,
      !!paste0(prefix, "CDR2") := CDR2,
      !!paste0(prefix, "FR3") := FR3,
      !!paste0(prefix, "CDR3") := CDR3,
      !!paste0(prefix, "FR4") := FR4,
      
      !!paste0(prefix, "FR1_aa") := FR1_aa,
      !!paste0(prefix, "CDR1_aa") := CDR1_aa,
      !!paste0(prefix, "FR2_aa") := FR2_aa,
      !!paste0(prefix, "CDR2_aa") := CDR2_aa,
      !!paste0(prefix, "FR3_aa") := FR3_aa,
      !!paste0(prefix, "CDR3_aa") := CDR3_aa,
      !!paste0(prefix, "FR4_aa") := FR4_aa,
      
      !!paste0(prefix, "V_call") := v_call,
      !!paste0(prefix, "D_call") := d_call,
      !!paste0(prefix, "J_call") := j_call,
      !!paste0(prefix, "C_call") := c_call,
      
      !!paste0(prefix, "productive") := productive,
      !!paste0(prefix, "stop_codon") := stop_codon,
      !!paste0(prefix, "vj_in_frame") := vj_in_frame,
      !!paste0(prefix, "v_frameshift") := v_frameshift,
      !!paste0(prefix, "j_frameshift") := j_frameshift,
      
      !!paste0(prefix, "sequence_alignment_aa") := sequence_alignment_aa,
      !!paste0(prefix, "v_sequence_alignment_aa") := v_sequence_alignment_aa,
      !!paste0(prefix, "j_sequence_alignment_aa") := j_sequence_alignment_aa,
      
      !!paste0(prefix, "v_identity") := first_existing_column(
        pick(everything()),
        "v_identity"
      ),
      !!paste0(prefix, "j_identity") := first_existing_column(
        pick(everything()),
        "j_identity"
      )
    )
}

annotate_scfv_igblast <- function(
    scFvs,
    work_dir = "scfv_igblast_work",
    threads = 4,
    min_host_score_margin = 5
) {
  dir.create(work_dir, showWarnings = FALSE, recursive = TRUE)
  
  slot_prefixes <- c("V1", "V2")
  slot_columns_to_refresh <- c(
    "_chain_type",
    "_locus",
    "_hit_score",
    "_frameshift",
    "_V",
    "_D",
    "_J",
    "_C",
    "_FR1",
    "_CDR1",
    "_FR2",
    "_CDR2",
    "_FR3",
    "_CDR3",
    "_FR4",
    "_FR1_aa",
    "_CDR1_aa",
    "_FR2_aa",
    "_CDR2_aa",
    "_FR3_aa",
    "_CDR3_aa",
    "_FR4_aa",
    "_V_call",
    "_D_call",
    "_J_call",
    "_C_call",
    "_productive",
    "_stop_codon",
    "_vj_in_frame",
    "_v_frameshift",
    "_j_frameshift",
    "_sequence_alignment_aa",
    "_v_sequence_alignment_aa",
    "_j_sequence_alignment_aa",
    "_v_identity",
    "_j_identity",
    "_pass_igblast_basic_qc",
    "_pass_igblast_productive_qc"
  )
  
  old_annotation_columns <- c(
    "host",
    "chain_config",
    "pass_igblast_basic_qc",
    "pass_igblast_productive_qc",
    paste0(
      rep(slot_prefixes, each = length(slot_columns_to_refresh)),
      rep(slot_columns_to_refresh, times = length(slot_prefixes))
    )
  )
  
  scFvs_base <- scFvs %>%
    select(-any_of(old_annotation_columns))
  
  query_tables <- make_query_tables(scFvs_base)
  
  fasta_file <- file.path(work_dir, "unique_scfv_v_sequences.fasta")
  write_fasta_table(query_tables$unique, fasta_file)
  
  human_config <- get_igblast_config("human")
  mouse_config <- get_igblast_config("mouse")
  
  human_airr <- run_igblast_airr(
    fasta_file = fasta_file,
    output_file = file.path(work_dir, "human_airr.tsv"),
    config = human_config,
    threads = threads
  )
  
  mouse_airr <- run_igblast_airr(
    fasta_file = fasta_file,
    output_file = file.path(work_dir, "mouse_airr.tsv"),
    config = mouse_config,
    threads = threads
  )
  
  annotation_table <- bind_rows(human_airr, mouse_airr) %>%
    extract_airr_annotation() %>%
    group_by(query_id, host_candidate) %>%
    arrange(desc(hit_score), .by_group = TRUE) %>%
    slice(1) %>%
    ungroup()
  
  pair_host <- pick_best_pair_host(
    long_table = query_tables$long,
    annotation_table = annotation_table,
    min_score_margin = min_host_score_margin
  )
  
  slot_annotation <- query_tables$long %>%
    left_join(
      pair_host,
      by = "row_id",
      relationship = "many-to-one"
    ) %>%
    left_join(
      annotation_table,
      by = c(
        "query_id" = "query_id",
        "host_for_annotation" = "host_candidate"
      ),
      relationship = "many-to-one"
    )
  
  v1_wide <- make_slot_wide(slot_annotation, "V1")
  v2_wide <- make_slot_wide(slot_annotation, "V2")
  
  final_table <- scFvs_base %>%
    mutate(row_id = dplyr::row_number()) %>%
    left_join(pair_host %>% select(row_id, host), by = "row_id") %>%
    left_join(v1_wide, by = "row_id") %>%
    left_join(v2_wide, by = "row_id") %>%
    mutate(
      V1_frameshift = format_airr_frame_status(
        V1_v_frameshift,
        V1_j_frameshift,
        V1_vj_in_frame
      ),
      V2_frameshift = format_airr_frame_status(
        V2_v_frameshift,
        V2_j_frameshift,
        V2_vj_in_frame
      ),
      
      V1_pass_igblast_basic_qc = pass_airr_domain_qc(
        V1_stop_codon,
        V1_vj_in_frame,
        V1_v_frameshift,
        V1_j_frameshift
      ),
      V2_pass_igblast_basic_qc = pass_airr_domain_qc(
        V2_stop_codon,
        V2_vj_in_frame,
        V2_v_frameshift,
        V2_j_frameshift
      ),
      
      V1_pass_igblast_productive_qc =
        V1_productive %in% TRUE &
        V1_pass_igblast_basic_qc,
      
      V2_pass_igblast_productive_qc =
        V2_productive %in% TRUE &
        V2_pass_igblast_basic_qc,
      
      pass_igblast_basic_qc =
        V1_pass_igblast_basic_qc &
        V2_pass_igblast_basic_qc,
      
      pass_igblast_productive_qc =
        V1_pass_igblast_productive_qc &
        V2_pass_igblast_productive_qc,
      
      chain_config = dplyr::case_when(
        V1_chain_type == "Heavy" &
          V2_chain_type %in% c("Light_kappa", "Light_lambda") ~ "HL",
        V1_chain_type %in% c("Light_kappa", "Light_lambda") &
          V2_chain_type == "Heavy" ~ "LH",
        TRUE ~ "Ambiguous"
      )
    ) %>%
    select(
      barcode,
      scFv_type,
      insert_V1_seq,
      insert_linker_seq,
      insert_V2_seq,
      host,
      chain_config,
      
      V1_chain_type,
      V1_frameshift,
      V1_productive,
      V1_stop_codon,
      V1_vj_in_frame,
      V1_v_frameshift,
      V1_j_frameshift,
      V1_pass_igblast_basic_qc,
      V1_pass_igblast_productive_qc,
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
      V1_V_call,
      V1_D_call,
      V1_J_call,
      V1_C_call,
      
      V2_chain_type,
      V2_frameshift,
      V2_productive,
      V2_stop_codon,
      V2_vj_in_frame,
      V2_v_frameshift,
      V2_j_frameshift,
      V2_pass_igblast_basic_qc,
      V2_pass_igblast_productive_qc,
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
      V2_FR4_aa,
      V2_V_call,
      V2_D_call,
      V2_J_call,
      V2_C_call,
      
      pass_igblast_basic_qc,
      pass_igblast_productive_qc,
      everything(),
      -row_id
    )
  
  final_table
}


# Post IgBLAST cleaning ----
is_present <- function(x) {
  !is.na(x) & x != ""
}

make_slot_exception_table <- function(exception_table, slot_name) {
  exception_table %>%
    filter(slot == slot_name) %>%
    select(-slot) %>%
    rename_with(
      ~ paste0(slot_name, "_exception_", .x),
      -scFv_type
    )
}

apply_igblast_qc_with_exceptions <- function(
    data,
    exception_table = scfv_igblast_exception_table
) {
  v1_exception <- make_slot_exception_table(exception_table, "V1")
  v2_exception <- make_slot_exception_table(exception_table, "V2")
  
  data %>%
    left_join(v1_exception, by = "scFv_type") %>%
    left_join(v2_exception, by = "scFv_type") %>%
    mutate(
      across(
        matches("^V[12]_exception_allow_"),
        ~ coalesce(.x, FALSE)
      ),
      
      V1_has_CDR3 = is_present(V1_CDR3),
      V2_has_CDR3 = is_present(V2_CDR3),
      V1_has_CDR3_aa = is_present(V1_CDR3_aa),
      V2_has_CDR3_aa = is_present(V2_CDR3_aa),
      
      V1_CDR3_ok =
        V1_has_CDR3 |
        V1_exception_allow_missing_cdr3,
      
      V2_CDR3_ok =
        V2_has_CDR3 |
        V2_exception_allow_missing_cdr3,
      
      V1_CDR3_aa_ok =
        V1_has_CDR3_aa |
        V1_exception_allow_missing_cdr3,
      
      V2_CDR3_aa_ok =
        V2_has_CDR3_aa |
        V2_exception_allow_missing_cdr3,
      
      V1_stop_ok =
        !(V1_stop_codon %in% TRUE) |
        V1_exception_allow_stop_codon,
      
      V2_stop_ok =
        !(V2_stop_codon %in% TRUE) |
        V2_exception_allow_stop_codon,
      
      V1_vj_in_frame_ok =
        V1_vj_in_frame %in% TRUE |
        (is.na(V1_vj_in_frame) & V1_exception_allow_vj_in_frame_na) |
        (V1_vj_in_frame %in% FALSE & V1_exception_allow_vj_in_frame_false),
      
      V2_vj_in_frame_ok =
        V2_vj_in_frame %in% TRUE |
        (is.na(V2_vj_in_frame) & V2_exception_allow_vj_in_frame_na) |
        (V2_vj_in_frame %in% FALSE & V2_exception_allow_vj_in_frame_false),
      
      V1_frameshift_ok =
        (
          !(V1_v_frameshift %in% TRUE) &
            !(V1_j_frameshift %in% TRUE)
        ) |
        V1_exception_allow_frameshift,
      
      V2_frameshift_ok =
        (
          !(V2_v_frameshift %in% TRUE) &
            !(V2_j_frameshift %in% TRUE)
        ) |
        V2_exception_allow_frameshift,
      
      V1_productive_ok =
        V1_productive %in% TRUE |
        V1_exception_allow_productive_false,
      
      V2_productive_ok =
        V2_productive %in% TRUE |
        V2_exception_allow_productive_false,
      
      V1_pass_igblast_qc =
        V1_CDR3_ok &
        V1_CDR3_aa_ok &
        V1_stop_ok &
        V1_vj_in_frame_ok &
        V1_frameshift_ok &
        V1_productive_ok,
      
      V2_pass_igblast_qc =
        V2_CDR3_ok &
        V2_CDR3_aa_ok &
        V2_stop_ok &
        V2_vj_in_frame_ok &
        V2_frameshift_ok &
        V2_productive_ok,
      
      pass_igblast_qc =
        V1_pass_igblast_qc &
        V2_pass_igblast_qc,
      
      V1_igblast_exception_applied =
        !is.na(V1_exception_reason) &
        (
          (!V1_has_CDR3 & V1_exception_allow_missing_cdr3) |
            (!V1_has_CDR3_aa & V1_exception_allow_missing_cdr3) |
            (!(V1_productive %in% TRUE) & V1_exception_allow_productive_false) |
            (is.na(V1_vj_in_frame) & V1_exception_allow_vj_in_frame_na) |
            (V1_vj_in_frame %in% FALSE & V1_exception_allow_vj_in_frame_false) |
            (
              (V1_v_frameshift %in% TRUE | V1_j_frameshift %in% TRUE) &
                V1_exception_allow_frameshift
            ) |
            (V1_stop_codon %in% TRUE & V1_exception_allow_stop_codon)
        ),
      
      V2_igblast_exception_applied =
        !is.na(V2_exception_reason) &
        (
          (!V2_has_CDR3 & V2_exception_allow_missing_cdr3) |
            (!V2_has_CDR3_aa & V2_exception_allow_missing_cdr3) |
            (!(V2_productive %in% TRUE) & V2_exception_allow_productive_false) |
            (is.na(V2_vj_in_frame) & V2_exception_allow_vj_in_frame_na) |
            (V2_vj_in_frame %in% FALSE & V2_exception_allow_vj_in_frame_false) |
            (
              (V2_v_frameshift %in% TRUE | V2_j_frameshift %in% TRUE) &
                V2_exception_allow_frameshift
            ) |
            (V2_stop_codon %in% TRUE & V2_exception_allow_stop_codon)
        ),
      
      igblast_exception_applied =
        V1_igblast_exception_applied |
        V2_igblast_exception_applied,
      
      igblast_exception_reason = case_when(
        V1_igblast_exception_applied & V2_igblast_exception_applied ~
          paste(V1_exception_reason, V2_exception_reason, sep = "; "),
        V1_igblast_exception_applied ~ V1_exception_reason,
        V2_igblast_exception_applied ~ V2_exception_reason,
        TRUE ~ NA_character_
      )
    )
}


translate_nt <- function(sequence, frame_start = 1) {
  codon_table <- c(
    TTT = "F", TTC = "F", TTA = "L", TTG = "L",
    TCT = "S", TCC = "S", TCA = "S", TCG = "S",
    TAT = "Y", TAC = "Y", TAA = "*", TAG = "*",
    TGT = "C", TGC = "C", TGA = "*", TGG = "W",
    CTT = "L", CTC = "L", CTA = "L", CTG = "L",
    CCT = "P", CCC = "P", CCA = "P", CCG = "P",
    CAT = "H", CAC = "H", CAA = "Q", CAG = "Q",
    CGT = "R", CGC = "R", CGA = "R", CGG = "R",
    ATT = "I", ATC = "I", ATA = "I", ATG = "M",
    ACT = "T", ACC = "T", ACA = "T", ACG = "T",
    AAT = "N", AAC = "N", AAA = "K", AAG = "K",
    AGT = "S", AGC = "S", AGA = "R", AGG = "R",
    GTT = "V", GTC = "V", GTA = "V", GTG = "V",
    GCT = "A", GCC = "A", GCA = "A", GCG = "A",
    GAT = "D", GAC = "D", GAA = "E", GAG = "E",
    GGT = "G", GGC = "G", GGA = "G", GGG = "G"
  )
  
  if (is.na(sequence) || sequence == "") {
    return(NA_character_)
  }
  
  sequence <- normalize_nt(sequence)
  
  if (nchar(sequence) < frame_start + 2) {
    return(NA_character_)
  }
  
  sequence <- substring(sequence, frame_start)
  full_codon_width <- floor(nchar(sequence) / 3) * 3
  
  if (full_codon_width < 3) {
    return(NA_character_)
  }
  
  sequence <- substring(sequence, 1, full_codon_width)
  
  codon_starts <- seq(1, nchar(sequence), by = 3)
  codons <- substring(sequence, codon_starts, codon_starts + 2)
  
  amino_acids <- unname(codon_table[codons])
  amino_acids[is.na(amino_acids)] <- "X"
  
  paste0(amino_acids, collapse = "")
}

clean_aa_sequence <- function(sequence) {
  sequence <- toupper(sequence)
  sequence <- gsub("[-.[:space:]]", "", sequence)
  sequence[!is.na(sequence) & sequence == ""] <- NA_character_
  sequence
}

has_internal_stop <- function(amino_acid_sequence, ignore_terminal_stop = TRUE) {
  if (is.na(amino_acid_sequence) || amino_acid_sequence == "") {
    return(NA)
  }
  
  if (ignore_terminal_stop) {
    amino_acid_sequence <- sub("\\*$", "", amino_acid_sequence)
  }
  
  grepl("\\*", amino_acid_sequence)
}

aa_contains <- function(parent_aa, child_aa) {
  if (
    is.na(parent_aa) || parent_aa == "" ||
    is.na(child_aa) || child_aa == ""
  ) {
    return(NA)
  }
  
  parent_aa <- clean_aa_sequence(parent_aa)
  child_aa <- clean_aa_sequence(child_aa)
  
  grepl(child_aa, parent_aa, fixed = TRUE)
}

aa_contains_any <- function(parent_aa, child_candidates) {
  child_candidates <- child_candidates[!is.na(child_candidates)]
  child_candidates <- child_candidates[child_candidates != ""]
  
  if (is.na(parent_aa) || parent_aa == "" || length(child_candidates) == 0) {
    return(NA)
  }
  
  any(vapply(
    child_candidates,
    function(child_aa) {
      aa_contains(parent_aa, child_aa) %in% TRUE
    },
    logical(1)
  ))
}

apply_insert_translation_qc <- function(
    data,
    ndei_site = "CATATG"
) {
  data %>%
    mutate(
      insert_V1_seq_norm = normalize_nt(insert_V1_seq),
      insert_seq_norm = normalize_nt(insert_seq),
      
      V1_starts_with_ndei =
        !is.na(insert_V1_seq_norm) &
        str_starts(insert_V1_seq_norm, ndei_site),
      
      insert_seq_aa = vapply(
        insert_seq_norm,
        translate_nt,
        character(1),
        frame_start = 1
      ),
      
      insert_seq_aa_frame2 = vapply(
        insert_seq_norm,
        translate_nt,
        character(1),
        frame_start = 2
      ),
      
      insert_seq_aa_frame3 = vapply(
        insert_seq_norm,
        translate_nt,
        character(1),
        frame_start = 3
      ),
      
      insert_seq_aa_has_internal_stop = vapply(
        insert_seq_aa,
        has_internal_stop,
        logical(1),
        ignore_terminal_stop = TRUE
      ),
      
      V1_alignment_aa_in_insert = pmap_lgl(
        list(
          insert_seq_aa,
          V1_sequence_alignment_aa,
          V1_v_sequence_alignment_aa
        ),
        function(parent_aa, sequence_alignment_aa, v_alignment_aa) {
          aa_contains_any(
            parent_aa,
            c(sequence_alignment_aa, v_alignment_aa)
          )
        }
      ),
      
      V2_alignment_aa_in_insert = pmap_lgl(
        list(
          insert_seq_aa,
          V2_sequence_alignment_aa,
          V2_v_sequence_alignment_aa
        ),
        function(parent_aa, sequence_alignment_aa, v_alignment_aa) {
          aa_contains_any(
            parent_aa,
            c(sequence_alignment_aa, v_alignment_aa)
          )
        }
      ),
      
      V1_alignment_aa_in_frame2 = pmap_lgl(
        list(
          insert_seq_aa_frame2,
          V1_sequence_alignment_aa,
          V1_v_sequence_alignment_aa
        ),
        function(parent_aa, sequence_alignment_aa, v_alignment_aa) {
          aa_contains_any(
            parent_aa,
            c(sequence_alignment_aa, v_alignment_aa)
          )
        }
      ),
      
      V2_alignment_aa_in_frame2 = pmap_lgl(
        list(
          insert_seq_aa_frame2,
          V2_sequence_alignment_aa,
          V2_v_sequence_alignment_aa
        ),
        function(parent_aa, sequence_alignment_aa, v_alignment_aa) {
          aa_contains_any(
            parent_aa,
            c(sequence_alignment_aa, v_alignment_aa)
          )
        }
      ),
      
      V1_alignment_aa_in_frame3 = pmap_lgl(
        list(
          insert_seq_aa_frame3,
          V1_sequence_alignment_aa,
          V1_v_sequence_alignment_aa
        ),
        function(parent_aa, sequence_alignment_aa, v_alignment_aa) {
          aa_contains_any(
            parent_aa,
            c(sequence_alignment_aa, v_alignment_aa)
          )
        }
      ),
      
      V2_alignment_aa_in_frame3 = pmap_lgl(
        list(
          insert_seq_aa_frame3,
          V2_sequence_alignment_aa,
          V2_v_sequence_alignment_aa
        ),
        function(parent_aa, sequence_alignment_aa, v_alignment_aa) {
          aa_contains_any(
            parent_aa,
            c(sequence_alignment_aa, v_alignment_aa)
          )
        }
      ),
      
      V1_alignment_aa_found_in_other_frame =
        V1_alignment_aa_in_frame2 |
        V1_alignment_aa_in_frame3,
      
      V2_alignment_aa_found_in_other_frame =
        V2_alignment_aa_in_frame2 |
        V2_alignment_aa_in_frame3,
      
      pass_insert_orf_qc =
        !is.na(insert_seq_aa) &
        !(insert_seq_aa_has_internal_stop %in% TRUE),
      
      pass_insert_alignment_qc =
        V1_alignment_aa_in_insert %in% TRUE &
        V2_alignment_aa_in_insert %in% TRUE,
      
      pass_insert_translation_qc =
        V1_starts_with_ndei &
        pass_insert_orf_qc &
        pass_insert_alignment_qc,
      
      insert_translation_qc_fail_reason = case_when(
        pass_insert_translation_qc ~ NA_character_,
        !V1_starts_with_ndei ~ "V1 does not start with CATATG",
        is.na(insert_seq_aa) ~ "insert_seq cannot be translated in frame 1",
        insert_seq_aa_has_internal_stop %in% TRUE ~
          "insert_seq_aa contains internal stop codon",
        !(V1_alignment_aa_in_insert %in% TRUE) &
          V1_alignment_aa_found_in_other_frame ~
          "V1 alignment AA is found in another frame",
        !(V2_alignment_aa_in_insert %in% TRUE) &
          V2_alignment_aa_found_in_other_frame ~
          "V2 alignment AA is found in another frame",
        !(V1_alignment_aa_in_insert %in% TRUE) ~
          "V1 alignment AA is not found in insert_seq_aa",
        !(V2_alignment_aa_in_insert %in% TRUE) ~
          "V2 alignment AA is not found in insert_seq_aa",
        TRUE ~ "Unknown insert translation QC failure"
      )
    ) %>%
    select(
      -insert_V1_seq_norm,
      -insert_seq_norm
    )
}


first_gene_call <- function(x) {
  x <- as.character(x)
  
  x %>%
    str_split(",", simplify = TRUE) %>%
    .[, 1] %>%
    str_remove("\\*.*$") %>%
    na_if("")
}

coalesce_text <- function(...) {
  dplyr::coalesce(..., "")
}

make_novel_clone_key <- function(
    host,
    chain_config,
    V1_chain_type,
    V1_V_gene,
    V1_J_gene,
    V1_CDR3_aa,
    V1_CDR3,
    V2_chain_type,
    V2_V_gene,
    V2_J_gene,
    V2_CDR3_aa,
    V2_CDR3
) {
  V1_CDR3_key <- dplyr::coalesce(V1_CDR3_aa, V1_CDR3, "NA")
  V2_CDR3_key <- dplyr::coalesce(V2_CDR3_aa, V2_CDR3, "NA")
  
  paste(
    "Novel",
    host,
    chain_config,
    V1_chain_type,
    V1_V_gene,
    V1_J_gene,
    V1_CDR3_key,
    V2_chain_type,
    V2_V_gene,
    V2_J_gene,
    V2_CDR3_key,
    sep = "|"
  )
}

to_logical_10x <- function(x) {
  if (is.logical(x)) {
    return(x)
  }
  
  x <- str_to_lower(as.character(x))
  
  case_when(
    x == "true" ~ TRUE,
    x == "false" ~ FALSE,
    TRUE ~ NA
  )
}

# Convert position string to sequence
extract_seq_by_position <- function(parent_seq, position) {
  if (is.na(parent_seq) || is.na(position) || position == "") {
    return(NA_character_)
  }
  
  position <- as.character(position)
  
  if (!grepl("^[0-9]+-[0-9]+$", position)) {
    return(position)
  }
  
  start <- as.integer(sub("-.*$", "", position))
  end <- as.integer(sub("^.*-", "", position))
  
  if (is.na(start) || is.na(end) || start < 1 || end < start) {
    return(NA_character_)
  }
  
  if (end > nchar(parent_seq)) {
    return(NA_character_)
  }
  
  substr(parent_seq, start, end)
}

extract_seq_by_position_vec <- function(parent_seq, position) {
  mapply(
    extract_seq_by_position,
    parent_seq,
    position,
    USE.NAMES = FALSE
  )
}

# Cut sequences ----
cutFrankingSite <- function(
    df,
    sequence_col = "sequence",
    left_target = "Myc_tag",
    right_target = "3_tail",
    left_criteria = c("score", "leftmost"),
    right_criteria = c("score", "rightmost"),
    flank_length = 5,
    min_score = NULL
) {
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
    
    if (length(scores) == 0) {
      return(c(start = NA_real_, end = NA_real_, score = NA_real_))
    }
    
    order_list <- lapply(criteria, function(x) {
      if (x == "score") {
        return(-scores)
      }
      
      if (x == "leftmost") {
        return(start)
      }
      
      if (x == "rightmost") {
        return(-end)
      }
      
      stop("Unknown criteria: ", x)
    })
    
    index <- do.call(order, order_list)[1]
    
    c(
      start = start[index],
      end = end[index],
      score = scores[index]
    )
  }
  
  left_count_col <- paste0(left_target, "_count")
  left_position_col <- paste0(left_target, "_position")
  left_score_col <- paste0(left_target, "_score")
  
  right_count_col <- paste0(right_target, "_count")
  right_position_col <- paste0(right_target, "_position")
  right_score_col <- paste0(right_target, "_score")
  
  result <- df %>%
    dplyr::filter(
      .data[[left_count_col]] > 0,
      .data[[right_count_col]] > 0
    )
  
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
  
  result %>%
    dplyr::bind_cols(left_hit, right_hit) %>%
    dplyr::mutate(
      sequence_length = nchar(.data[[sequence_col]]),
      valid_franking_order = left_hit_end < right_hit_start,
      
      left_franking_seq = dplyr::if_else(
        valid_franking_order,
        stringr::str_sub(
          .data[[sequence_col]],
          pmax(1, left_hit_start - flank_length),
          left_hit_end
        ),
        NA_character_
      ),
      
      insert_seq = dplyr::if_else(
        valid_franking_order,
        stringr::str_sub(
          .data[[sequence_col]],
          left_hit_end + 1,
          right_hit_start - 1
        ),
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
      insert_seq_length = nchar(insert_seq),
    ) %>%
    dplyr::select(
      dplyr::all_of(sequence_col),
      left_franking_seq,
      insert_seq,
      right_franking_seq,
      dplyr::everything()
    )
}


cutLinker <- function(
    df,
    sequence_col = "insert_seq",
    linker_position_col = "linker_position"
) {
  if (!sequence_col %in% colnames(df)) {
    stop("Cannot find sequence_col: ", sequence_col)
  }
  
  if (!linker_position_col %in% colnames(df)) {
    stop("Cannot find linker_position_col: ", linker_position_col)
  }
  
  df %>%
    dplyr::mutate(
      linker_start = as.integer(
        stringr::str_extract(.data[[linker_position_col]], "^[0-9]+")
      ),
      linker_end = as.integer(
        stringr::str_extract(.data[[linker_position_col]], "(?<=-)[0-9]+")
      ),
      
      insert_V1_seq = stringr::str_sub(
        .data[[sequence_col]],
        1,
        linker_start - 1
      ),
      
      insert_linker_seq = stringr::str_sub(
        .data[[sequence_col]],
        linker_start,
        linker_end
      ),
      
      insert_V2_seq = stringr::str_sub(
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
      insert_V1_seq,
      insert_linker_seq,
      insert_V2_seq,
      dplyr::everything()
    )
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
  if (!length_col %in% colnames(df)) {
    stop("Cannot find length_col: ", length_col)
  }
  
  if (length(colors) != 4) {
    stop("colors must be length 4:\n",
         "c(hist_fill, hist_border, density_fill, density_line)")
  }
  
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
    annotation_col = "sequence_name"
) {
  data_list <- list(...)
  
  if (is.null(names(data_list)) || any(names(data_list) == "")) {
    stop("Please provide named tables.")
  }
  
  if (!"raw" %in% names(data_list)) {
    stop("Please provide raw table.")
  }
  
  if (!"annotation" %in% names(data_list)) {
    stop("Please provide annotation table.")
  }
  
  annotation_df <- data_list$annotation
  filter_list <- data_list[names(data_list) != "annotation"]
  
  if (!annotation_col %in% colnames(annotation_df)) {
    stop("Cannot find annotation_col: ", annotation_col)
  }
  
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
    dplyr::select(
      stage,
      group,
      n,
      parent_n,
      total_n,
      percent_parent,
      percent_total
    )
  
  filter_summary %>%
    dplyr::select(
      stage,
      group,
      n,
      percent_parent,
      percent_total
    ) %>%
    dplyr::bind_rows(
      annotation_summary %>%
        dplyr::select(
          stage,
          group,
          n,
          percent_parent,
          percent_total
        )
    ) %>%
    dplyr::mutate(
      stage = dplyr::case_when(
        stage == "raw" ~ "Raw",
        stage == "annotation" ~ "Annotation",
        TRUE ~ paste0(
          "Filter: ",
          stringr::str_remove(stage, "_filter$")
        )
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
    seu_meatdata,
    seu_data,
    keep_cols = NULL,
    gene_set_pattern = "^Gene set [0-9]+$",
    keep_cell_without_exp = FALSE,
    barcode_col = "barcode"
) {
  stop_missing_cols <- function(data, required_cols, data_name) {
    missing_cols <- setdiff(required_cols, colnames(data))
    
    if (length(missing_cols) > 0) {
      stop(
        data_name,
        " is missing required columns: ",
        paste(missing_cols, collapse = ", "),
        call. = FALSE
      )
    }
  }
  
  resolve_keep_cols <- function(default_cols, keep_cols, available_cols) {
    default_cols <- intersect(default_cols, available_cols)
    
    if (is.null(keep_cols) || length(keep_cols) == 0) {
      return(default_cols)
    }
    
    keep_cols <- as.character(keep_cols)
    
    remove_cols <- keep_cols[startsWith(keep_cols, "-")]
    remove_cols <- sub("^-", "", remove_cols)
    
    add_cols <- keep_cols[!startsWith(keep_cols, "-")]
    
    unknown_add_cols <- setdiff(add_cols, available_cols)
    unknown_remove_cols <- setdiff(remove_cols, available_cols)
    
    if (length(unknown_add_cols) > 0) {
      warning(
        "These keep_cols were not found and will be ignored: ",
        paste(unknown_add_cols, collapse = ", ")
      )
    }
    
    if (length(unknown_remove_cols) > 0) {
      warning(
        "These columns requested for removal were not found and will be ignored: ",
        paste(unknown_remove_cols, collapse = ", ")
      )
    }
    
    final_cols <- unique(c(default_cols, intersect(add_cols, available_cols)))
    final_cols <- setdiff(final_cols, remove_cols)
    
    final_cols
  }
  
  stop_missing_cols(scFv_annotation, barcode_col, "scFv_annotation")
  stop_missing_cols(seu_meatdata, barcode_col, "seu_meatdata")
  stop_missing_cols(seu_data, barcode_col, "seu_data")
  
  gene_set_cols <- grep(
    gene_set_pattern,
    colnames(seu_meatdata),
    value = TRUE
  )
  
  expression_cols <- setdiff(colnames(seu_data), barcode_col)
  
  default_scfv_cols <- c(
    barcode_col,
    "scFv_type",
    "clone_type",
    "clone_ID",
    "clone_CDR3_set",
    "clone_CDR3_set_aa",
    "host",
    "chain_config",
    "linker_name",
    "quality_score",
    
    "insert_seq",
    "insert_seq_aa",
    "insert_V1_seq",
    "insert_V1_seq_aa",
    "insert_linker_seq",
    "insert_linker_seq_aa",
    "insert_V2_seq",
    "insert_V2_seq_aa",
    
    "V1_chain_type",
    "V1_V_call",
    "V1_D_call",
    "V1_J_call",
    "V1_C_call",
    "V1_V_gene",
    "V1_D_gene",
    "V1_J_gene",
    "V1_C_gene",
    "V1_FR1",
    "V1_FR1_aa",
    "V1_CDR1",
    "V1_CDR1_aa",
    "V1_FR2",
    "V1_FR2_aa",
    "V1_CDR2",
    "V1_CDR2_aa",
    "V1_FR3",
    "V1_FR3_aa",
    "V1_CDR3",
    "V1_CDR3_aa",
    "V1_FR4",
    "V1_FR4_aa",
    
    "V2_chain_type",
    "V2_V_call",
    "V2_D_call",
    "V2_J_call",
    "V2_C_call",
    "V2_V_gene",
    "V2_D_gene",
    "V2_J_gene",
    "V2_C_gene",
    "V2_FR1",
    "V2_FR1_aa",
    "V2_CDR1",
    "V2_CDR1_aa",
    "V2_FR2",
    "V2_FR2_aa",
    "V2_CDR2",
    "V2_CDR2_aa",
    "V2_FR3",
    "V2_FR3_aa",
    "V2_CDR3",
    "V2_CDR3_aa",
    "V2_FR4",
    "V2_FR4_aa"
  )
  
  default_metadata_cols <- c(
    barcode_col,
    "orig.ident",
    gene_set_cols
  )
  
  default_expression_cols <- c(
    barcode_col,
    expression_cols
  )
  
  available_cols <- unique(c(
    colnames(scFv_annotation),
    colnames(seu_meatdata),
    colnames(seu_data)
  ))
  
  default_cols <- unique(c(
    default_scfv_cols,
    default_metadata_cols,
    default_expression_cols
  ))
  
  final_keep_cols <- resolve_keep_cols(
    default_cols = default_cols,
    keep_cols = keep_cols,
    available_cols = available_cols
  )
  
  scfv_cols <- intersect(final_keep_cols, colnames(scFv_annotation))
  metadata_cols <- intersect(final_keep_cols, colnames(seu_meatdata))
  expression_keep_cols <- intersect(final_keep_cols, colnames(seu_data))
  
  # 避免 join 後產生 .x / .y。左側資料優先。
  metadata_cols <- setdiff(metadata_cols, setdiff(scfv_cols, barcode_col))
  expression_keep_cols <- setdiff(
    expression_keep_cols,
    setdiff(c(scfv_cols, metadata_cols), barcode_col)
  )
  
  if (!(barcode_col %in% scfv_cols)) {
    scfv_cols <- c(barcode_col, scfv_cols)
  }
  
  if (!(barcode_col %in% metadata_cols)) {
    metadata_cols <- c(barcode_col, metadata_cols)
  }
  
  if (!(barcode_col %in% expression_keep_cols)) {
    expression_keep_cols <- c(barcode_col, expression_keep_cols)
  }
  
  metadata_tbl <- seu_meatdata %>%
    select(all_of(metadata_cols)) %>%
    distinct(.data[[barcode_col]], .keep_all = TRUE)
  
  expression_tbl <- seu_data %>%
    select(all_of(expression_keep_cols)) %>%
    distinct(.data[[barcode_col]], .keep_all = TRUE) %>%
    mutate(.has_expression = TRUE)
  
  training_data <- scFv_annotation %>%
    select(all_of(scfv_cols)) %>%
    left_join(metadata_tbl, by = barcode_col) %>%
    left_join(expression_tbl, by = barcode_col)
  
  if (!keep_cell_without_exp) {
    training_data <- training_data %>%
      filter(.has_expression %in% TRUE)
  }
  
  training_data %>%
    select(-any_of(".has_expression"))
}