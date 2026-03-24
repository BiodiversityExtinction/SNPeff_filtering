#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(scales)
})

usage <- function() {
  cat(
    "Usage:\n",
    "  Rscript plot_snpeff_load.R --load FILE --groups FILE [options]\n\n",
    "Required arguments:\n",
    "  --load FILE       Load summary TSV from summarise_snpeff_load.sh\n",
    "  --groups FILE     Two-column sample/population file\n\n",
    "Optional arguments:\n",
    "  --outdir DIR      Output directory (default: current directory)\n",
    "  --prefix NAME     Output prefix (default: derived from --load)\n",
    "  --width NUM       Plot width in inches (default: 14)\n",
    "  --height NUM      Plot height in inches (default: 7)\n",
    "  --dpi NUM         PNG dpi (default: 300)\n",
    "  --help            Show this help\n\n",
    "What it does:\n",
    "  1. Filters the load table to samples present in the group file\n",
    "  2. Makes one LoF plot and one nonsynonymous plot\n",
    "  3. Plots realised values on the left y-axis and masked values on the right y-axis\n",
    "  4. Adds individual confidence intervals if present in the load table\n",
    "  5. Calculates pairwise population Rxy for the four main metrics\n",
    "  6. If there are two populations, makes a pairwise-value Rxy boxplot\n",
    "  7. If there are more than two populations, makes a faceted Rxy heatmap\n\n",
    "Expected group file format:\n",
    "  sample<TAB>population\n\n",
    sep = ""
  )
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0 || any(args %in% c("-h", "--help"))) {
  usage()
  quit(status = 0)
}

parse_args <- function(args) {
  out <- list(
    load = NULL,
    groups = NULL,
    outdir = ".",
    prefix = NULL,
    width = 14,
    height = 7,
    dpi = 300
  )

  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (key == "--load") {
      out$load <- args[[i + 1]]
      i <- i + 2
    } else if (key == "--groups") {
      out$groups <- args[[i + 1]]
      i <- i + 2
    } else if (key == "--outdir") {
      out$outdir <- args[[i + 1]]
      i <- i + 2
    } else if (key == "--prefix") {
      out$prefix <- args[[i + 1]]
      i <- i + 2
    } else if (key == "--width") {
      out$width <- as.numeric(args[[i + 1]])
      i <- i + 2
    } else if (key == "--height") {
      out$height <- as.numeric(args[[i + 1]])
      i <- i + 2
    } else if (key == "--dpi") {
      out$dpi <- as.numeric(args[[i + 1]])
      i <- i + 2
    } else {
      stop("Unknown argument: ", key, call. = FALSE)
    }
  }
  out
}

opt <- parse_args(args)

if (is.null(opt$load) || is.null(opt$groups)) {
  usage()
  stop("Both --load and --groups are required.", call. = FALSE)
}

if (!file.exists(opt$load)) {
  stop("Load file not found: ", opt$load, call. = FALSE)
}
if (!file.exists(opt$groups)) {
  stop("Groups file not found: ", opt$groups, call. = FALSE)
}

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

if (is.null(opt$prefix)) {
  prefix <- basename(opt$load)
  prefix <- sub("\\.tsv$", "", prefix)
  prefix <- sub("\\.txt$", "", prefix)
  opt$prefix <- prefix
}

detect_column <- function(df, choices) {
  found <- choices[choices %in% names(df)]
  if (length(found) == 0) {
    stop("Could not find any of these columns: ", paste(choices, collapse = ", "), call. = FALSE)
  }
  found[[1]]
}

load_df <- readr::read_tsv(opt$load, show_col_types = FALSE, progress = FALSE)
groups_df <- readr::read_tsv(
  opt$groups,
  col_names = c("sample", "population"),
  comment = "#",
  show_col_types = FALSE,
  progress = FALSE
)

groups_df <- groups_df %>%
  filter(!is.na(sample), !is.na(population)) %>%
  filter(sample != "sample")

load_df <- load_df %>%
  inner_join(groups_df, by = "sample")

if (nrow(load_df) == 0) {
  stop("No overlapping samples were found between the load table and the group file.", call. = FALSE)
}

cols <- list(
  lof_realised = detect_column(load_df, c("lof_realised_over_synonymous_total", "lof_realized_over_synonymous_total")),
  lof_masked = "lof_masked_over_synonymous_total",
  nonsyn_realised = detect_column(load_df, c("nonsynonymous_realised_over_synonymous_total", "nonsynonymous_realized_over_synonymous_total")),
  nonsyn_masked = "nonsynonymous_masked_over_synonymous_total",
  lof_realised_low = detect_column(load_df, c("lof_realised_over_synonymous_total_jackknife_ci_lower", "lof_realized_over_synonymous_total_jackknife_ci_lower")),
  lof_realised_high = detect_column(load_df, c("lof_realised_over_synonymous_total_jackknife_ci_upper", "lof_realized_over_synonymous_total_jackknife_ci_upper")),
  lof_masked_low = "lof_masked_over_synonymous_total_jackknife_ci_lower",
  lof_masked_high = "lof_masked_over_synonymous_total_jackknife_ci_upper",
  nonsyn_realised_low = detect_column(load_df, c("nonsynonymous_realised_over_synonymous_total_jackknife_ci_lower", "nonsynonymous_realized_over_synonymous_total_jackknife_ci_lower")),
  nonsyn_realised_high = detect_column(load_df, c("nonsynonymous_realised_over_synonymous_total_jackknife_ci_upper", "nonsynonymous_realized_over_synonymous_total_jackknife_ci_upper")),
  nonsyn_masked_low = "nonsynonymous_masked_over_synonymous_total_jackknife_ci_lower",
  nonsyn_masked_high = "nonsynonymous_masked_over_synonymous_total_jackknife_ci_upper"
)

load_df <- load_df %>%
  mutate(population = factor(population, levels = unique(population))) %>%
  arrange(population, sample) %>%
  mutate(sample = factor(sample, levels = sample))

make_metric_df <- function(df, title, realised_col, masked_col, realised_low, realised_high, masked_low, masked_high) {
  bind_rows(
    df %>%
      transmute(
        sample,
        population,
        state = "Realised",
        value = .data[[realised_col]],
        ci_lower = .data[[realised_low]],
        ci_upper = .data[[realised_high]]
      ),
    df %>%
      transmute(
        sample,
        population,
        state = "Masked",
        value = .data[[masked_col]],
        ci_lower = .data[[masked_low]],
        ci_upper = .data[[masked_high]]
      )
  ) %>%
    mutate(metric = title)
}

plot_dual_axis <- function(metric_df, metric_name, outfile_base, width, height, dpi) {
  realised_df <- metric_df %>% filter(state == "Realised")
  masked_df <- metric_df %>% filter(state == "Masked")

  realised_max <- max(c(realised_df$value, realised_df$ci_upper), na.rm = TRUE)
  masked_max <- max(c(masked_df$value, masked_df$ci_upper), na.rm = TRUE)
  scale_factor <- if (is.finite(masked_max) && masked_max > 0) realised_max / masked_max else 1
  if (!is.finite(scale_factor) || scale_factor <= 0) {
    scale_factor <- 1
  }

  realised_plot <- realised_df %>%
    mutate(plot_value = value, plot_low = ci_lower, plot_high = ci_upper)
  masked_plot <- masked_df %>%
    mutate(plot_value = value * scale_factor, plot_low = ci_lower * scale_factor, plot_high = ci_upper * scale_factor)

  plot_df <- bind_rows(realised_plot, masked_plot)

  p <- ggplot(plot_df, aes(x = sample, y = plot_value, colour = population, shape = state)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
    geom_linerange(aes(ymin = plot_low, ymax = plot_high), alpha = 0.7, position = position_dodge(width = 0.5)) +
    geom_point(size = 2.2, position = position_dodge(width = 0.5)) +
    scale_shape_manual(values = c("Realised" = 16, "Masked" = 17)) +
    scale_y_continuous(
      name = paste0(metric_name, " realised / synonymous_total"),
      labels = label_number(accuracy = 0.001),
      sec.axis = sec_axis(~ . / scale_factor, name = paste0(metric_name, " masked / synonymous_total"), labels = label_number(accuracy = 0.001))
    ) +
    labs(
      x = "Sample",
      colour = "Population",
      shape = "State",
      title = paste0(metric_name, " Load Per Individual")
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      panel.grid.major.x = element_blank(),
      plot.title = element_text(face = "bold")
    )

  ggsave(paste0(outfile_base, ".pdf"), p, width = width, height = height, units = "in")
  ggsave(paste0(outfile_base, ".png"), p, width = width, height = height, units = "in", dpi = dpi)
}

lof_df <- make_metric_df(
  load_df,
  title = "LoF",
  realised_col = cols$lof_realised,
  masked_col = cols$lof_masked,
  realised_low = cols$lof_realised_low,
  realised_high = cols$lof_realised_high,
  masked_low = cols$lof_masked_low,
  masked_high = cols$lof_masked_high
)

nonsyn_df <- make_metric_df(
  load_df,
  title = "Nonsynonymous",
  realised_col = cols$nonsyn_realised,
  masked_col = cols$nonsyn_masked,
  realised_low = cols$nonsyn_realised_low,
  realised_high = cols$nonsyn_realised_high,
  masked_low = cols$nonsyn_masked_low,
  masked_high = cols$nonsyn_masked_high
)

plot_dual_axis(
  lof_df,
  metric_name = "LoF",
  outfile_base = file.path(opt$outdir, paste0(opt$prefix, ".lof_dual_axis")),
  width = opt$width,
  height = opt$height,
  dpi = opt$dpi
)

plot_dual_axis(
  nonsyn_df,
  metric_name = "Nonsynonymous",
  outfile_base = file.path(opt$outdir, paste0(opt$prefix, ".nonsynonymous_dual_axis")),
  width = opt$width,
  height = opt$height,
  dpi = opt$dpi
)

rxy_metrics <- c(
  cols$lof_realised,
  cols$lof_masked,
  cols$nonsyn_realised,
  cols$nonsyn_masked
)

population_means <- load_df %>%
  group_by(population) %>%
  summarise(across(all_of(rxy_metrics), ~ mean(.x, na.rm = TRUE)), n_samples = n(), .groups = "drop")

rxy_rows <- list()
for (metric in rxy_metrics) {
  for (pop_x in population_means$population) {
    for (pop_y in population_means$population) {
      if (pop_x == pop_y) {
        next
      }
      mean_x <- population_means %>% filter(population == pop_x) %>% pull(all_of(metric))
      mean_y <- population_means %>% filter(population == pop_y) %>% pull(all_of(metric))
      rxy <- ifelse(is.na(mean_x) || is.na(mean_y) || mean_y == 0, NA_real_, mean_x / mean_y)
      rxy_rows[[length(rxy_rows) + 1]] <- tibble(
        metric = metric,
        population_x = as.character(pop_x),
        population_y = as.character(pop_y),
        mean_x = mean_x,
        mean_y = mean_y,
        rxy = rxy
      )
    }
  }
}

rxy_df <- bind_rows(rxy_rows)

metric_labels <- c(
  "lof_realised_over_synonymous_total" = "LoF realised /\nsynonymous",
  "lof_realized_over_synonymous_total" = "LoF realised /\nsynonymous",
  "lof_masked_over_synonymous_total" = "LoF masked /\nsynonymous",
  "nonsynonymous_realised_over_synonymous_total" = "Nonsynonymous realised /\nsynonymous",
  "nonsynonymous_realized_over_synonymous_total" = "Nonsynonymous realised /\nsynonymous",
  "nonsynonymous_masked_over_synonymous_total" = "Nonsynonymous masked /\nsynonymous"
)

metric_order <- c(
  "LoF realised /\nsynonymous",
  "LoF masked /\nsynonymous",
  "Nonsynonymous realised /\nsynonymous",
  "Nonsynonymous masked /\nsynonymous"
)

if (nlevels(load_df$population) == 2) {
  pop_levels <- levels(load_df$population)
  pop_x <- pop_levels[[1]]
  pop_y <- pop_levels[[2]]

  rxy_pair_rows <- list()
  for (metric in rxy_metrics) {
    x_vals <- load_df %>% filter(population == pop_x) %>% pull(all_of(metric))
    y_vals <- load_df %>% filter(population == pop_y) %>% pull(all_of(metric))
    pair_grid <- expand.grid(x_value = x_vals, y_value = y_vals, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    pair_grid$rxy <- ifelse(pair_grid$y_value == 0, NA_real_, pair_grid$x_value / pair_grid$y_value)
    rxy_pair_rows[[length(rxy_pair_rows) + 1]] <- tibble(
      metric = metric,
      metric_label = factor(recode(metric, !!!metric_labels), levels = metric_order),
      comparison = paste0(pop_x, " / ", pop_y),
      x_value = pair_grid$x_value,
      y_value = pair_grid$y_value,
      rxy = pair_grid$rxy
    )
  }
  rxy_pair_df <- bind_rows(rxy_pair_rows)

  observed_rxy_df <- rxy_df %>%
    filter(population_x == pop_x, population_y == pop_y) %>%
    mutate(
      metric_label = factor(recode(metric, !!!metric_labels), levels = metric_order),
      comparison = paste0(pop_x, " / ", pop_y)
    )

  rxy_plot <- ggplot(rxy_pair_df, aes(x = metric_label, y = rxy, fill = comparison)) +
    geom_hline(yintercept = 1, linewidth = 0.4, linetype = "dashed", colour = "grey35") +
    geom_boxplot(width = 0.6, alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.12, alpha = 0.15, size = 0.8, colour = "grey35") +
    geom_point(
      data = observed_rxy_df,
      aes(x = metric_label, y = rxy, colour = comparison),
      inherit.aes = FALSE,
      size = 3,
      stroke = 0.7,
      show.legend = TRUE
    ) +
    labs(
      x = NULL,
      y = "Rxy",
      fill = "Comparison",
      colour = "Comparison",
      title = paste0("Pairwise Rxy: ", pop_x, " / ", pop_y),
      subtitle = "Boxplots show all between-population pairwise sample ratios"
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 25, hjust = 1),
      panel.grid.major.x = element_blank(),
      plot.title = element_text(face = "bold"),
      legend.position = "bottom"
    )

  ggsave(
    file.path(opt$outdir, paste0(opt$prefix, ".rxy_boxplot.pdf")),
    rxy_plot,
    width = 8,
    height = 6,
    units = "in"
  )
  ggsave(
    file.path(opt$outdir, paste0(opt$prefix, ".rxy_boxplot.png")),
    rxy_plot,
    width = 8,
    height = 6,
    units = "in",
    dpi = opt$dpi
  )

  readr::write_tsv(rxy_pair_df, file.path(opt$outdir, paste0(opt$prefix, ".rxy_pairwise.tsv")))
} else {
  rxy_plot_df <- rxy_df %>%
    mutate(
      metric_label = recode(metric, !!!metric_labels),
      population_x = factor(population_x, levels = levels(load_df$population)),
      population_y = factor(population_y, levels = rev(levels(load_df$population)))
    )

  rxy_plot <- ggplot(rxy_plot_df, aes(x = population_y, y = population_x, fill = rxy)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = ifelse(is.na(rxy), "NA", sprintf("%.3f", rxy))), size = 3) +
    facet_wrap(~ metric_label, ncol = 2) +
    scale_fill_gradient2(
      low = "#2166ac",
      mid = "white",
      high = "#b2182b",
      midpoint = 1,
      na.value = "grey85"
    ) +
    labs(
      x = "Population Y",
      y = "Population X",
      fill = "Rxy",
      title = "Pairwise Population Rxy"
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank(),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold")
    )

  ggsave(
    file.path(opt$outdir, paste0(opt$prefix, ".rxy_heatmap.pdf")),
    rxy_plot,
    width = max(10, opt$width),
    height = max(8, opt$height),
    units = "in"
  )
  ggsave(
    file.path(opt$outdir, paste0(opt$prefix, ".rxy_heatmap.png")),
    rxy_plot,
    width = max(10, opt$width),
    height = max(8, opt$height),
    units = "in",
    dpi = opt$dpi
  )
}

readr::write_tsv(load_df %>% select(sample, population), file.path(opt$outdir, paste0(opt$prefix, ".samples_used.tsv")))
readr::write_tsv(population_means, file.path(opt$outdir, paste0(opt$prefix, ".population_means.tsv")))
readr::write_tsv(rxy_df, file.path(opt$outdir, paste0(opt$prefix, ".rxy.tsv")))

message("Finished.")
message("Plots:")
message("  ", file.path(opt$outdir, paste0(opt$prefix, ".lof_dual_axis.pdf")))
message("  ", file.path(opt$outdir, paste0(opt$prefix, ".nonsynonymous_dual_axis.pdf")))
if (nlevels(load_df$population) == 2) {
  message("  ", file.path(opt$outdir, paste0(opt$prefix, ".rxy_boxplot.pdf")))
} else {
  message("  ", file.path(opt$outdir, paste0(opt$prefix, ".rxy_heatmap.pdf")))
}
message("Tables:")
message("  ", file.path(opt$outdir, paste0(opt$prefix, ".samples_used.tsv")))
message("  ", file.path(opt$outdir, paste0(opt$prefix, ".population_means.tsv")))
message("  ", file.path(opt$outdir, paste0(opt$prefix, ".rxy.tsv")))
if (nlevels(load_df$population) == 2) {
  message("  ", file.path(opt$outdir, paste0(opt$prefix, ".rxy_pairwise.tsv")))
}
