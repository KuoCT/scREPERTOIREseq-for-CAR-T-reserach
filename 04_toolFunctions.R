# Save images ----
saveSVG <- function(obj, filename, skip = FALSE, ...) {
  if (skip & file.exists(filename)) return ("Skip: file already exists.")
  svglite::svglite(file = filename, ...)
  
  on.exit(grDevices::dev.off(), add = TRUE)
  
  if (inherits(obj, "ggplot")) {
    print(obj)
  } else if (inherits(obj, c("Heatmap", "HeatmapList"))) {
    ComplexHeatmap::draw(obj)
  } else {
    stop(
      "Unsupported object type: ",
      paste(class(obj), collapse = ", "),
      call. = FALSE
    )
  }
  
  invisible(filename)
}


savePNG <- function(obj, filename, width = 6, height = 4, dpi = 300, skip = FALSE, ...) {
  if (skip & file.exists(filename)) return ("Skip: file already exists.")
  grDevices::png(filename = filename, width = width, height = height, units = "in", res = dpi, ...)
  
  on.exit(grDevices::dev.off(), add = TRUE)
  
  if (inherits(obj, "ggplot")) {
    print(obj)
  } else if (inherits(obj, c("Heatmap", "HeatmapList"))) {
    ComplexHeatmap::draw(obj)
  } else {
    stop(
      "Unsupported object type: ",
      paste(class(obj), collapse = ", "),
      call. = FALSE
    )
  }
  
  invisible(filename)
}


# Preset for ggplot ----
presetAdd <- function(...) {
  unlist(list(...), recursive = FALSE)
}


ggPreset <- function(style = "default") {
  default <- list(
    ggprism::theme_prism(base_size = 10),
    theme(
      plot.title = element_text(hjust = 0),
      plot.subtitle = element_text(hjust = 0),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.ticks = element_line(colour = "black"),
      legend.frame = element_rect(colour = "black")
    )
  )
  
  DimPlot <- list(
    guides(
      color = guide_legend(
        nrow = 2,
        byrow = TRUE,
        override.aes = list(size = 2.8)
      )
    )
  )
  
  DimPlot2 <- list(
    guides(
      color = guide_legend(
        nrow = 5,
        byrow = TRUE,
        override.aes = list(size = 2.8)
      )
    )
  )
  
  NoLegend <- list(
    theme(
      legend.position = "none"
    )
  )
  
  presets <- list(
    default = default,
    DimPlot = presetAdd(default, DimPlot),
    DimPlot2 = presetAdd(default, DimPlot2),
    NoLegend = presetAdd(default, NoLegend)
  )
  
  if (!style %in% names(presets)) {
    stop(
      "Unknown style: ", style,
      ". Available styles: ",
      paste(names(presets), collapse = ", "),
      call. = FALSE
    )
  }
  
  presets[[style]]
}


# PlotSet manager for grouping plots ----
# Internal store for the active plot set
.PLOTSET <- new.env(parent = emptyenv())
.PLOTSET$ps <- NULL
# Create or reuse the active PlotSet
newPlotSet <- function(overwrite_without_ask = FALSE) {
  has_existing <- !is.null(.PLOTSET$ps)
  
  if (!has_existing) {
    ps <- new.env(parent = emptyenv())
    ps$plots <- list()
    .PLOTSET$ps <- ps
    return(invisible(ps))
  }
  
  if (isTRUE(overwrite_without_ask)) {
    ps <- new.env(parent = emptyenv())
    ps$plots <- list()
    .PLOTSET$ps <- ps
    return(invisible(ps))
  }
  
  cat(
    "A PlotSet already exists.\n",
    "[1] overwrite\n",
    "[2] append\n",
    sep = ""
  )
  ans <- readline(prompt = "Select: ")
  ans <- trimws(ans)
  
  if (ans == "1") {
    ps <- new.env(parent = emptyenv())
    ps$plots <- list()
    .PLOTSET$ps <- ps
  } else if (ans == "2") {
    # Keep current PlotSet
  } else {
    stop("Invalid input. Please enter 1 or 2.", call. = FALSE)
  }
  
  invisible(.PLOTSET$ps)
}


# Return active PlotSet
activePlotSet <- function() {
  if (is.null(.PLOTSET$ps)) {
    stop(
      "No active PlotSet found. Please run newPlotSet() first.",
      call. = FALSE
    )
  }
  
  .PLOTSET$ps
}


# Add plot to active PlotSet and return the plot itself
.nextPlotFallbackName <- function(existing_names) {
  p_names <- grep("^p[0-9]+$", existing_names, value = TRUE)
  
  if (length(p_names) == 0) {
    return("p1")
  }
  
  p_ids <- as.integer(sub("^p", "", p_names))
  paste0("p", max(p_ids, na.rm = TRUE) + 1)
}


.makeUniquePlotName <- function(name, existing_names) {
  if (!name %in% existing_names) {
    return(name)
  }
  
  i <- 2
  repeat {
    candidate <- paste0(name, "_", i)
    if (!candidate %in% existing_names) {
      return(candidate)
    }
    i <- i + 1
  }
}


addPlot <- function(plot, name = NULL) {
  ps <- activePlotSet()
  existing_names <- names(ps$plots)
  
  if (is.null(name)) {
    auto_name <- NULL
    
    if (inherits(plot, "ggplot")) {
      auto_name <- tryCatch(
        ggplot2::get_labs(plot)$title,
        error = function(e) NULL
      )
    }
    
    if (!is.null(auto_name)) {
      auto_name <- tryCatch(
        as.character(auto_name),
        error = function(e) NULL
      )
    }
    
    if (
      is.null(auto_name) ||
      length(auto_name) == 0 ||
      !nzchar(trimws(auto_name[1]))
    ) {
      name <- .nextPlotFallbackName(existing_names)
    } else {
      name <- auto_name[1]
    }
  }
  
  name <- .makeUniquePlotName(name, existing_names)
  
  ps$plots[[name]] <- plot
  plot
}


# List plot names
listPlotSet <- function() {
  ps <- activePlotSet()
  names(ps$plots)
}


# Return the whole PlotSet as a named list
PlotSet <- function() {
  ps <- activePlotSet()
  ps$plots
}


# Internal plot drawing helper
.drawPlot <- function(plot) {
  if (inherits(plot, c("Heatmap", "HeatmapList"))) {
    ComplexHeatmap::draw(plot)
  } else {
    print(plot)
  }
  
  invisible(plot)
}


# Get plot by name or index
# getPlot() or getPlot(-1) returns the last plot
getPlot <- function(index = NULL, show = TRUE) {
  ps <- activePlotSet()
  
  if (length(ps$plots) == 0) {
    stop("The active PlotSet is empty.", call. = FALSE)
  }
  
  plot_names <- names(ps$plots)
  
  if (is.null(index) || identical(index, -1)) {
    plot <- ps$plots[[length(ps$plots)]]
  } else if (is.character(index) && length(index) == 1) {
    if (!index %in% plot_names) {
      stop("Plot name not found: ", index, call. = FALSE)
    }
    plot <- ps$plots[[index]]
  } else if (is.numeric(index) && length(index) == 1) {
    if (index < 1 || index > length(ps$plots)) {
      stop("Plot index out of range.", call. = FALSE)
    }
    plot <- ps$plots[[index]]
  } else {
    stop(
      "`index` must be NULL, -1, a plot name, or a single numeric index.",
      call. = FALSE
    )
  }
  
  if (show) {
    .drawPlot(plot)
  }
  
  invisible(plot)
}


# Remove active PlotSet
closePlotSet <- function() {
  .PLOTSET$ps <- NULL
  invisible(NULL)
}


# Improve UMAP-like scatter plots ----
improveUMAP <- function(
    plot,
    x = "umap_1",
    y = "umap_2",
    point_size = 1,
    outline_expand = 0.8,
    outline_colour = "black",
    rasterize = FALSE,
    raster_dpi = 600,
    raster_dev = "cairo",
    transparent_background = TRUE
) {
  # Add outline layer beneath existing point layers
  outline_layer <- ggplot2::geom_point(
    data = plot$data,
    mapping = ggplot2::aes(x = .data[[x]], y = .data[[y]]),
    inherit.aes = FALSE,
    colour = outline_colour,
    size = point_size + outline_expand
  )
  
  # Insert outline as the first layer
  plot$layers <- append(plot$layers, list(outline_layer), after = 0)
  
  # Rasterize all point layers if requested
  if (isTRUE(rasterize)) {
    plot$layers <- lapply(
      plot$layers,
      function(layer) {
        if (inherits(layer$geom, "GeomPoint")) {
          ggrastr::rasterise(
            layer,
            dpi = raster_dpi,
            dev = raster_dev
          )
        } else {
          layer
        }
      }
    )
  }
  
  # Make backgrounds transparent if requested
  if (isTRUE(transparent_background)) {
    plot <- plot +
      ggplot2::theme(
        panel.background = ggplot2::element_rect(
          fill = "transparent",
          colour = NA
        ),
        plot.background = ggplot2::element_rect(
          fill = "transparent",
          colour = NA
        ),
        legend.background = ggplot2::element_rect(
          fill = "transparent",
          colour = NA
        ),
        legend.box.background = ggplot2::element_rect(
          fill = "transparent",
          colour = NA
        )
      )
  }
  
  plot
}


# Data.frame helper ----
# Convert data.frame rownames to a column
df_r2c <- function(df, rowname_col) {
  if (missing(rowname_col)) {
    stop("Please provide 'rowname_col', e.g. df_r2c(df, 'gene').")
  }
  
  if (!is.data.frame(df)) {
    stop("'df' must be a data.frame.")
  }
  
  if (rowname_col %in% colnames(df)) {
    stop(paste0("Column '", rowname_col, "' already exists in df."))
  }
  
  out_df <- cbind(
    setNames(
      data.frame(rownames(df), check.names = FALSE),
      rowname_col
    ),
    df
  )
  
  rownames(out_df) <- NULL
  
  return(out_df)
}


# Convert one data.frame column to rownames
df_c2r <- function(df, rowname_col) {
  if (missing(rowname_col)) {
    stop("Please provide 'rowname_col', e.g. df_c2r(df, 'gene').")
  }
  
  if (!is.data.frame(df)) {
    stop("'df' must be a data.frame.")
  }
  
  if (!rowname_col %in% colnames(df)) {
    stop(paste0("Column '", rowname_col, "' not found in df."))
  }
  
  row_names <- df[[rowname_col]]
  
  if (any(is.na(row_names))) {
    stop("The selected rowname column contains NA values.")
  }
  
  if (anyDuplicated(row_names)) {
    stop("The selected rowname column contains duplicated values.")
  }
  
  out_df <- df[, setdiff(colnames(df), rowname_col), drop = FALSE]
  rownames(out_df) <- row_names
  
  return(out_df)
}


# Convert matrix rownames to a data.frame column ----
mat2df <- function(mat, rowname_col) {
  if (missing(rowname_col)) {
    stop("Please provide 'rowname_col', e.g. mat2df(mat, 'gene').")
  }
  
  if (!is.matrix(mat)) {
    stop("'mat' must be a matrix.")
  }
  
  if (is.null(rownames(mat))) {
    stop("'mat' must have rownames.")
  }
  
  out_df <- data.frame(
    mat,
    check.names = FALSE
  )
  
  out_df <- df_r2c(out_df, rowname_col)
  
  return(out_df)
}


# Convert data.frame to matrix ----
# If rowname_col is provided, use that column as matrix rownames.
# If rowname_col is not provided, keep existing data.frame rownames.
df2mat <- function(df, rowname_col = NULL) {
  if (!is.data.frame(df)) {
    stop("'df' must be a data.frame.")
  }
  
  if (!is.null(rowname_col)) {
    df <- df_c2r(df, rowname_col)
  }
  
  mat <- as.matrix(df)
  
  suppressWarnings(mode(mat) <- "numeric")
  
  if (any(is.na(mat))) {
    warning(
      "NA values were introduced during numeric conversion. ",
      "Please check whether non-numeric columns exist."
    )
  }
  
  return(mat)
}


# Save RDS with overwrite check ----
saveRDS2 <- function(object, file, overwrite = FALSE, ...) {
  if (file.exists(file) && !isTRUE(overwrite)) {
    message("File exists. Skip saving: ", file)
    return(invisible(FALSE))
  }
  
  saveRDS(object, file = file, ...)
  message("Saved: ", file)
  
  invisible(TRUE)
}