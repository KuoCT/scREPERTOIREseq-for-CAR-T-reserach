# Standard Seurat QC, filtering, normalization, and feature selection ----
seurat_qc_flow <- function(
    obj,
    out_path = "",
    min_nFeature_RNA = 200,
    max_nFeature_RNA = 7500,
    max_percent_mt = 10,
    plot_limit_nFeature_RNA = c(NA_real_, NA_real_),  # e.g. c(0, 10000), c(NA_real_, 10000), c(0, NA_real_)
    plot_limit_nCount_RNA = c(NA_real_, NA_real_),    # e.g. c(0, 45000), c(NA_real_, 45000), c(0, NA_real_)
    plot_limit_percent_mt = c(NA_real_, NA_real_),    # e.g. c(0, 35), c(NA_real_, 35), c(0, NA_real_)
    variable_nfeatures = 3000,
    mt_pattern = "^MT-",
    violin_adjust = 2,
    violin_width = 8,
    violin_height = 5,
    scatter_width = 9,
    scatter_height = 4
) {
  # Set RNA as the default assay
  DefaultAssay(obj) <- "RNA"
  
  # Calculate the percentage of mitochondrial gene expression
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = mt_pattern)
  
  # Extract project name for output file naming
  obj_name <- Project(obj)
  
  # Create output directory only when out_path is provided
  save_plot <- !is.null(out_path) && nzchar(out_path)
  if (save_plot && !dir.exists(out_path)) {
    dir.create(out_path, recursive = TRUE)
  }
  
  # Apply y-axis limits with support for partial limits
  add_y_limits <- function(p, limits) {
    lower <- limits[1]
    upper <- limits[2]
    
    if (all(is.na(c(lower, upper)))) {
      return(p)
    }
    
    p + coord_cartesian(ylim = c(lower, upper))
  }
  
  # Apply x/y-axis limits with support for partial limits
  add_xy_limits <- function(p, xlim_vals = c(NA_real_, NA_real_),
                            ylim_vals = c(NA_real_, NA_real_)) {
    has_x <- !all(is.na(xlim_vals))
    has_y <- !all(is.na(ylim_vals))
    
    if (!has_x && !has_y) {
      return(p)
    }
    
    if (has_x && has_y) {
      return(p + coord_cartesian(xlim = xlim_vals, ylim = ylim_vals))
    }
    
    if (has_x) {
      return(p + coord_cartesian(xlim = xlim_vals))
    }
    
    p + coord_cartesian(ylim = ylim_vals)
  }
  
  # QC violin plots
  qc_features <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
  qc_plot_limits <- list(
    plot_limit_nFeature_RNA,
    plot_limit_nCount_RNA,
    plot_limit_percent_mt
  )
  
  # Avoid generating an extra invalid violin polygon in SVG
  update_geom_defaults("violin", aes(linewidth = NA, alpha = 0))
  
  violin_plots <- lapply(seq_along(qc_features), function(i) {
    p <- VlnPlot(
      object = obj,
      features = qc_features[i],
      pt.size = 0,
      adjust = violin_adjust
    ) +
      geom_jitter(width = 0.45, height = 0, alpha = 0.03) +
      geom_violin(
        adjust = violin_adjust,
        linewidth = 1,
        alpha = 0.8
      ) +
      theme_prism(base_size = 12) +
      NoLegend()
    
    add_y_limits(p, qc_plot_limits[[i]])
  })
  
  p_violin <- wrap_plots(
    violin_plots,
    ncol = 3,
    widths = c(1, 1.05, 1.08)
  )
  
  # Always print plots
  print(p_violin)
  
  # Save violin plot only when out_path is provided
  if (save_plot) {
    violin_file <- file.path(out_path, paste0(obj_name, "_QC_vlnplot.svg"))
    svglite(violin_file, width = violin_width, height = violin_height)
    print(p_violin)
    dev.off()
    
    # Remove duplicated violin polygon objects in the SVG output
    # This is a workaround for malformed SVG elements caused by geom_violin
    svg_lines <- readLines(violin_file)
    pattern <- paste0(
      "polygon\\s+points=.*style=['\"]?[^'\"]*",
      "stroke-width:\\s*nan;\\s*stroke:\\s*none;",
      "\\s*stroke-linecap:\\s*butt[^'\"]*['\"]?"
    )
    svg_lines <- svg_lines[!grepl(pattern, svg_lines)]
    writeLines(svg_lines, violin_file)
  }
  
  # QC scatter plots
  p_scatter_mt <- FeatureScatter(
    obj,
    feature1 = "nCount_RNA",
    feature2 = "percent.mt",
    shuffle = TRUE,
    pt.size = 1.5
  ) +
    theme_prism(base_size = 10)
  
  p_scatter_mt <- add_xy_limits(
    p_scatter_mt,
    xlim_vals = plot_limit_nCount_RNA,
    ylim_vals = plot_limit_percent_mt
  )
  
  p_scatter_feature <- FeatureScatter(
    obj,
    feature1 = "nCount_RNA",
    feature2 = "nFeature_RNA",
    shuffle = TRUE,
    pt.size = 1.5
  ) +
    theme_prism(base_size = 10)
  
  p_scatter_feature <- add_xy_limits(
    p_scatter_feature,
    xlim_vals = plot_limit_nCount_RNA,
    ylim_vals = plot_limit_nFeature_RNA
  )
  
  p_scatter <- p_scatter_mt + p_scatter_feature
  
  # Always print plots
  print(p_scatter)
  
  # Save scatter plot only when out_path is provided
  if (save_plot) {
    scatter_file <- file.path(
      out_path,
      paste0(obj_name, "_QC_Scatterplot.svg")
    )
    svglite(scatter_file, width = scatter_width, height = scatter_height)
    print(p_scatter)
    dev.off()
  }
  
  # Filter cells based on QC thresholds
  obj <- subset(
    obj,
    subset =
      nFeature_RNA > min_nFeature_RNA &
      nFeature_RNA < max_nFeature_RNA &
      percent.mt < max_percent_mt
  )
  
  # Normalize the data
  obj <- NormalizeData(obj)
  
  # Identify highly variable features
  obj <- FindVariableFeatures(obj, nfeatures = variable_nfeatures)
  
  # Scale all features
  obj <- ScaleData(obj)
  
  return(obj)
}


# Rename and reorder Seurat clusters ----
RenameSeuratClusters <- function(
    seu_object, 
    cluster_names, 
    order_levels = NULL, 
    new_column_name = "renamed_clusters"
) {
  # Ensure the number of provided cluster names matches the number of clusters
  if (length(cluster_names) != length(unique(Idents(seu_object)))) {
    stop("The number of cluster names must match the number of clusters in the object.")
  }
  
  # Create a named list for renaming clusters
  cluster_mapping <- setNames(cluster_names, levels(Idents(seu_object)))
  
  # Rename the clusters
  seu_object <- RenameIdents(seu_object, cluster_mapping)
  
  # Reorder the clusters if `order_levels` is provided
  if (!is.null(order_levels)) {
    if (!all(order_levels %in% cluster_names)) {
      stop("All order levels must be present in the cluster names.")
    }
    Idents(seu_object) <- factor(Idents(seu_object), levels = order_levels)
  }
  
  # Save the renamed clusters in a new metadata column
  seu_object[[new_column_name]] <- Idents(seu_object)
  
  return(seu_object)
}


# Modified FeaturePlot ----
FeaturePlot2 <- function(
    obj, 
    features, 
    min.cutoff = NA,
    max.cutoff = NA, 
    colours = replace(viridis_pal(option = "B")(10), 1, "gray90"),
    colour_limits = NULL,
    dim_x_range = NULL,
    dim_y_range = NULL,
    x_limits = NULL,
    x_breaks = waiver(),
    y_limits = NULL,
    y_breaks = waiver()
) {
  plots <- list()
  
  for (feature in features) {
    p <- FeaturePlot(
      obj,
      features = feature,
      order = TRUE,
      pt.size = 1,
      min.cutoff = min.cutoff,
      max.cutoff = max.cutoff
    ) +
      scale_color_gradientn(
        colours = colours, 
        limits = colour_limits, 
        name = "Scaled Exp."
      ) +
      coord_cartesian(xlim = dim_x_range, ylim = dim_y_range) +
      scale_x_continuous(breaks = x_breaks, limits = x_limits) +
      scale_y_continuous(breaks = y_breaks, limits = y_limits) + 
      theme_prism(base_size = 10) +
      theme(
        plot.title = element_text(hjust = 0),
        plot.subtitle = element_text(hjust = 0),
        legend.title = element_text(vjust = 0.75),
        legend.ticks = element_line(colour = "black"),
        legend.frame = element_rect(colour = "black"),
        legend.position = "bottom"
      )
    
    plots[[feature]] <- p
  }
  
  if (length(features) == 1) {
    return(plots[[1]])  # ggplot object
  } else {
    return(plots)       # list of ggplot objects
  }
}

# Modified AddModuleScore ----
AddModuleScore2 <- function(obj, markers, ScoreName, assay = "RNA") {
  existing_scores <- colnames(obj@meta.data)
  if (ScoreName %in% existing_scores) {
    obj@meta.data <- obj@meta.data[, !(colnames(obj@meta.data) %in% ScoreName)]
  }
  
  obj <- AddModuleScore(
    obj,
    features = list(markers),
    name = ScoreName,
    assay = assay
  )
  
  colnames(obj@meta.data) <- gsub(
    pattern = paste0(ScoreName, 1),
    replacement = ScoreName,
    colnames(obj@meta.data)
  )
  
  return(obj)
}


# Append metadata to Seurat object ----
AppendSeuratMetadata <- function(
    seu_object,
    df_to_Append,
    barcode_col = "barcode",
    conflict_action = c("overwrite", "make_unique"),
    fill_na_cols = NULL,
    fill_na_value = "Undefined"
) {
  conflict_action <- match.arg(conflict_action)
  
  metadata <- seu_object[[]] %>%
    tibble::rownames_to_column(barcode_col)
  
  overlap_cols <- intersect(
    colnames(metadata),
    colnames(df_to_Append)
  )
  overlap_cols <- setdiff(overlap_cols, barcode_col)
  
  if (length(overlap_cols) > 0) {
    if (conflict_action == "overwrite") {
      metadata <- metadata %>%
        dplyr::select(-dplyr::all_of(overlap_cols))
    } else if (conflict_action == "make_unique") {
      df_to_Append <- df_to_Append %>%
        dplyr::rename_with(
          ~ paste0(.x, "_new"),
          dplyr::all_of(overlap_cols)
        )
    }
  }
  
  metadata <- metadata %>%
    dplyr::left_join(df_to_Append, by = barcode_col)
  
  if (!is.null(fill_na_cols)) {
    metadata <- metadata %>%
      dplyr::mutate(
        dplyr::across(
          dplyr::all_of(fill_na_cols),
          ~ dplyr::if_else(is.na(.x), fill_na_value, as.character(.x))
        )
      )
  }
  
  metadata <- metadata %>%
    tibble::column_to_rownames(barcode_col)
  
  seu_object <- Seurat::AddMetaData(
    seu_object,
    metadata = metadata
  )
  
  return(seu_object)
}

# Color GSEA Enrichment Score Plot ----
improveGSEAplot <- function(
    svg_path,
    score_path,
    output_folder,
    y_pos = NULL,
    height = NULL,
    x_start = NULL,
    x_end = NULL,
    palette = c("navy", "white", "firebrick3"),
    hit_height_scale = 0.7,
    border_color = "black",
    border_width = 0.75
) {
  # Read SVG source code
  svg_lines <- readLines(svg_path, warn = FALSE)
  
  # Auto-detect ranked list metric range from SVG
  if (is.null(x_start) || is.null(x_end)) {
    polygon_text <- grep("<polygon\\s+points=", svg_lines, value = TRUE)
    
    polygon_values <- lapply(polygon_text, function(x) {
      points_text <- sub(".*<polygon[^>]*points=['\"]([^'\"]+)['\"].*", "\\1", x)
      points <- strsplit(trimws(points_text), "\\s+")[[1]]
      xy <- do.call(rbind, strsplit(points, ","))
      
      data.frame(
        x = suppressWarnings(as.numeric(xy[, 1])),
        y = suppressWarnings(as.numeric(xy[, 2]))
      )
    })
    
    polygon_values <- polygon_values[vapply(polygon_values, function(x) {
      sum(!is.na(x$x) & !is.na(x$y)) > 20 && diff(range(x$x, na.rm = TRUE)) > 100
    }, logical(1))]
    
    if (length(polygon_values) == 0) stop("Cannot auto-detect ranked list metric polygon from SVG: ", svg_path)
    
    point_counts <- vapply(polygon_values, nrow, integer(1))
    rank_polygon <- polygon_values[[which.max(point_counts)]]
    
    if (is.null(x_start)) x_start <- min(rank_polygon$x, na.rm = TRUE)
    if (is.null(x_end)) x_end <- max(rank_polygon$x, na.rm = TRUE)
  }
  
  # Detect hitting gene vertical lines
  line_idx <- grep("<line\\s+", svg_lines)
  line_text <- svg_lines[line_idx]
  
  if (length(line_idx) > 0) {
    line_values <- lapply(line_text, function(x) {
      vals <- sub(".*x1=['\"]([0-9.]+)['\"].*y1=['\"]([0-9.]+)['\"].*x2=['\"]([0-9.]+)['\"].*y2=['\"]([0-9.]+)['\"].*", "\\1 \\2 \\3 \\4", x)
      as.numeric(strsplit(vals, "\\s+")[[1]])
    })
    
    line_values <- do.call(rbind, line_values)
    colnames(line_values) <- c("x1", "y1", "x2", "y2")
    
    line_df <- data.frame(
      line_id = line_idx,
      line_text = line_text,
      x1 = line_values[, "x1"],
      y1 = line_values[, "y1"],
      x2 = line_values[, "x2"],
      y2 = line_values[, "y2"],
      stringsAsFactors = FALSE
    )
    
    # Pick the repeated solid vertical lines used for hitting genes
    hit_candidates <- line_df[
      abs(line_df$x1 - line_df$x2) < 1e-6 &
        abs(line_df$y1 - line_df$y2) > 20 &
        line_df$y1 > 200 & line_df$y2 > 200 &
        grepl("stroke-width:\\s*1\\.50", line_df$line_text) &
        !grepl("stroke-dasharray", line_df$line_text),
      ,
      drop = FALSE
    ]
    
    if (nrow(hit_candidates) > 0) {
      pair_key <- paste(
        round(pmax(hit_candidates$y1, hit_candidates$y2), 2),
        round(pmin(hit_candidates$y1, hit_candidates$y2), 2),
        sep = "_"
      )
      
      hit_key <- names(sort(table(pair_key), decreasing = TRUE))[1]
      hit_lines <- hit_candidates[pair_key == hit_key, , drop = FALSE]
      
      hit_bottom_old <- unique(round(pmax(hit_lines$y1, hit_lines$y2), 4))
      hit_top_old <- unique(round(pmin(hit_lines$y1, hit_lines$y2), 4))
      
      hit_bottom_old <- hit_bottom_old[1]
      hit_top_old <- hit_top_old[1]
      
      hit_center <- (hit_top_old + hit_bottom_old) / 2
      hit_height_old <- hit_bottom_old - hit_top_old
      hit_height_new <- hit_height_old * hit_height_scale
      hit_top_new <- hit_center - hit_height_new / 2
      hit_bottom_new <- hit_center + hit_height_new / 2
      
      # Update hitting gene line height
      for (k in seq_len(nrow(hit_lines))) {
        idx <- hit_lines$line_id[k]
        old_y1 <- hit_lines$y1[k]
        old_y2 <- hit_lines$y2[k]
        
        new_y1 <- if (old_y1 > old_y2) hit_bottom_new else hit_top_new
        new_y2 <- if (old_y1 > old_y2) hit_top_new else hit_bottom_new
        
        svg_lines[idx] <- sub(
          "(y1=['\"])[0-9.]+(['\"])",
          paste0("\\1", format(round(new_y1, 2), nsmall = 2), "\\2"),
          svg_lines[idx]
        )
        svg_lines[idx] <- sub(
          "(y2=['\"])[0-9.]+(['\"])",
          paste0("\\1", format(round(new_y2, 2), nsmall = 2), "\\2"),
          svg_lines[idx]
        )
      }
      
      # Make color bar match hitting gene line range unless manually set
      if (is.null(y_pos)) y_pos <- hit_top_new
      if (is.null(height)) height <- hit_height_new
    }
  }
  
  # Fallback values if hitting gene lines cannot be detected
  if (is.null(y_pos)) y_pos <- 300
  if (is.null(height)) height <- 25
  
  # Read score table
  df <- read.csv(score_path, row.names = NULL, sep = "\t")
  if (!"score" %in% colnames(df)) stop("Score table does not contain column: score")
  
  # Order score
  df <- df[order(-df$score), ]
  scores <- df$score
  n <- length(scores)
  if (n == 0) stop("No score values found in: ", score_path)
  
  # Color mapping
  lower_q <- quantile(scores, 0.10, na.rm = TRUE)
  upper_q <- quantile(scores, 0.90, na.rm = TRUE)
  mid <- 0
  scores_clipped <- scales::squish(scores, range = c(lower_q, upper_q))
  
  neg_fn <- scales::col_numeric(palette = palette[1:2], domain = c(lower_q, mid))
  pos_fn <- scales::col_numeric(palette = palette[2:3], domain = c(mid, upper_q))
  
  mapped_colors <- vapply(scores_clipped, function(s) {
    if (s == mid) {
      palette[2]
    } else if (s < mid) {
      neg_fn(s)
    } else {
      pos_fn(s)
    }
  }, character(1))
  
  # Coordinate mapping
  rect_width <- (x_end - x_start) / n
  x_coords <- x_start + seq(0, n - 1) * rect_width
  total_width <- x_end - x_start
  
  # Build color bar
  heatmap_rects <- paste0(
    "<rect x='", round(x_coords, 4),
    "' y='", round(y_pos, 2),
    "' width='", round(rect_width, 4),
    "' height='", round(height, 2),
    "' fill='", mapped_colors,
    "' stroke='none' />"
  )
  
  # Add outer border
  heatmap_border <- paste0(
    "<rect x='", round(x_start, 4),
    "' y='", round(y_pos, 2),
    "' width='", round(total_width, 4),
    "' height='", round(height, 2),
    "' fill='none' stroke='", border_color,
    "' stroke-width='", border_width, "' />"
  )
  
  heatmap_svg <- c("<g id='heatmap-bar'>", heatmap_rects, heatmap_border, "</g>")
  
  # Insert heatmap bar near top so existing black lines stay visible above it
  insert_at <- grep("<svg[^>]*>", svg_lines)[1]
  new_svg <- append(svg_lines, heatmap_svg, after = insert_at)
  
  # Create output folder
  if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
  
  # Write SVG
  file_base <- sub("\\.svg$", "", basename(svg_path))
  output_file <- file.path(output_folder, paste0(file_base, "_color.svg"))
  writeLines(new_svg, output_file)
  
  message("Heatmap bar added and saved to: ", output_file)
  invisible(output_file)
}