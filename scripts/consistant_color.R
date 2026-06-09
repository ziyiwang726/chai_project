library(ggrepel)
library(stringr)

#################### color helpers #######################
blend_color <- function(col, target = "white", amount = 0.4) {
  a <- grDevices::col2rgb(col) / 255
  b <- grDevices::col2rgb(target) / 255
  out <- (1 - amount) * a + amount * b
  grDevices::rgb(out[1,], out[2,], out[3,])
}

variant_shades <- function(base_col, K = 12) {
  grDevices::colorRampPalette(c(
    blend_color(base_col, "white", 0.55),
    base_col,
    blend_color(base_col, "black", 0.25)
  ))(K)
}

make_method_palette <- function(methods, K = 12) {
  fixed <- c(
    chai        = "red",
    adapt_glm   = "#FF7F0E",
    adapt_gmm_z = "#2CA02C",
    adapt_gmm_p = "#9467BD",
    FDRreg      = "#8C564B",
    BH          = "black"
  )
  
  family_base <- c(
    OrderShapeEM = "#1F77B4",  # blue family
    IHW          = "#17BECF"   # teal family
  )
  
  ord_shades <- variant_shades(family_base["OrderShapeEM"], K = K)
  ihw_shades <- variant_shades(family_base["IHW"],          K = K)
  
  methods <- unique(as.character(methods))
  pal <- setNames(rep("#7F7F7F", length(methods)), methods)  # fallback grey
  
  for (m in methods) {
    if (m %in% names(fixed)) {
      pal[m] <- fixed[m]
      next
    }
    
    base <- sub("_x\\d+$", "", m)
    idx  <- suppressWarnings(as.integer(sub(".*_x", "", m)))
    
    if (base %in% names(family_base)) {
      if (is.na(idx)) {
        pal[m] <- family_base[base]  # e.g. "OrderShapeEM" with no suffix
      } else {
        idx <- max(1L, min(idx, K))
        pal[m] <- if (base == "OrderShapeEM") ord_shades[idx] else ihw_shades[idx]
      }
    }
  }
  
  pal
}

method_order <- function(methods) {
  methods <- unique(as.character(methods))
  core <- c("chai","adapt_glm","adapt_gmm_z","adapt_gmm_p","FDRreg")
  
  sort_family <- function(base) {
    fam <- methods[grepl(paste0("^", base, "(_x\\d+)?$"), methods)]
    base_present <- base %in% fam
    vars <- fam[grepl("_x\\d+$", fam)]
    if (length(vars)) {
      idx <- suppressWarnings(as.integer(sub(".*_x", "", vars)))
      vars <- vars[order(idx)]
    }
    c(if (base_present) base else character(0), vars)
  }
  
  out <- c(core, sort_family("OrderShapeEM"), sort_family("IHW"), "BH")
  out <- unique(out)
  out[out %in% methods]
}


################## plotting function (For simulations with informativeness parameter a) ##################
# For Simulation 1 and Simulation S1
plot_metric_aq <- function(long_table, metric_name, x_var = c("q","a"),
                           ylab, title,
                           K_variants = 12,
                           fdp_y = 0.2,
                           x_breaks = NULL,
                           q_level = 0.05,
                           method_font = 5.5) {
  
  x_var <- match.arg(x_var)
  
  dat <- long_table %>%
    dplyr::filter(Metric == metric_name) %>%
    mutate(Method = as.character(Method))
  
  ord <- method_order(dat$Method)
  dat <- dat %>% mutate(Method = factor(Method, levels = ord))
  pal <- make_method_palette(levels(dat$Method), K = K_variants)
  
  # endpoints for labels = last x (either max q or max a)
  endpoints <- dat %>%
    group_by(Method) %>%
    filter(.data[[x_var]] == max(.data[[x_var]], na.rm = TRUE)) %>%
    slice_tail(n = 1) %>%
    ungroup() %>%
    mutate(label = as.character(Method))
  
  y_max <- if (metric_name == "FDP") fdp_y else 1
  
  endpoints <- endpoints %>%
    mutate(
      Value_label = pmin(Value, y_max)   # clamp label anchor to the axis top
    )
  
  x_rng   <- range(dat[[x_var]], na.rm = TRUE)
  x_nudge <- diff(x_rng) * 0.03
  
  p <- ggplot(dat, aes(x = .data[[x_var]], y = Value, group = Method, color = Method)) +
    geom_line(linewidth = 0.9, na.rm = TRUE) +
    geom_point(size = 1.6, na.rm = TRUE) +
    ggrepel::geom_text_repel(
      data = endpoints,
      aes(x = .data[[x_var]], y = Value_label, label = label),
      size = method_font,
      nudge_x = x_nudge,
      hjust = 0,
      direction = "y",
      box.padding = 0.25,
      point.padding = 0.15,
      segment.size = 0.25,
      min.segment.length = 0,
      max.overlaps = Inf,
      seed = 1,
      show.legend = FALSE
    ) +
    scale_color_manual(values = pal, breaks = ord) +
    coord_cartesian(ylim = c(0, 1), clip = "off") +
    labs(x = x_var, y = ylab, title = title) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.margin = margin(5.5, 60, 5.5, 5.5),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      axis.text.x  = element_text(size = 12),
      axis.text.y  = element_text(size = 12)
    )
  
  if (!is.null(x_breaks)) {
    p <- p + scale_x_continuous(breaks = x_breaks,
                                expand = expansion(mult = c(0.02, 0.18)))
  } else {
    p <- p + scale_x_continuous(expand = expansion(mult = c(0.02, 0.18)))
  }
  
  if (metric_name == "FDP" && x_var == "q") {
    p <- p +
      annotate(
        "segment",
        x = 0, y = 0, xend = 0.1, yend = 0.1,
        linetype = "dashed",
        linewidth = 1.2,
        color = "black"
      ) +
      coord_cartesian(ylim = c(0, fdp_y))
  }
  
  if (metric_name == "FDP" && x_var == "a") {
    p <- p +
      geom_hline(
        yintercept = q_level,
        linetype = "dashed",
        linewidth = 1.2,
        color = "black"
      ) +
      coord_cartesian(ylim = c(0, fdp_y))
  }
  
  p
}

  
################## plotting function (For simulations without informativeness parameter a) ##################
# For Simulation 2 and Simulation S2
plot_metric <- function(long_table, metric_name, 
                        ylab, title, K_variants = 12, 
                        fdp_y = 0.2, method_font = 5.5){
  
  dat <- dplyr::filter(long_table, Metric == metric_name) %>%
    mutate(Method = as.character(Method))
  
  ord <- method_order(dat$Method)
  dat <- dat %>% mutate(Method = factor(Method, levels = ord))
  pal <- make_method_palette(levels(dat$Method), K = K_variants)
  
  # endpoints for labels
  endpoints <- dat %>%
    group_by(Method) %>%
    filter(q == max(q, na.rm = TRUE)) %>%
    slice_tail(n = 1) %>%
    ungroup() %>%
    mutate(
      # option A: wrap long names
      label = as.character(Method)
      # label = str_wrap(label, width = 12)
    )
  
  x_rng <- range(dat$q, na.rm = TRUE)
  x_nudge <- diff(x_rng) * 0.03  # more room to the right
  
  p <- ggplot(dat, aes(x = q, y = Value, group = Method, color = Method)) +
    geom_line(linewidth = 0.9, na.rm = TRUE) +
    geom_point(size = 1.6, na.rm = TRUE) +
    
    # repel labels to avoid overlap
    ggrepel::geom_text_repel(
      data = endpoints,
      aes(label = label),
      size = method_font,
      nudge_x = x_nudge,
      hjust = 0,
      direction = "y",
      box.padding = 0.25,
      point.padding = 0.15,
      segment.size = 0.25,
      min.segment.length = 0,
      max.overlaps = Inf,
      seed = 1,
      show.legend = FALSE
    ) +
    
    scale_color_manual(values = pal, breaks = ord) +
    scale_x_continuous(
      breaks = seq(0.01, 0.10, 0.01),
      expand = expansion(mult = c(0.02, 0.18))  # add right padding for labels
    ) +
    coord_cartesian(
      ylim = c(0, 1),
      clip = "off"   # <-- critical: do not clip labels
    ) +
    labs(x = "q level", y = ylab, title = title) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.margin = margin(5.5, 60, 5.5, 5.5),  # <-- more right margin
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      axis.text.x  = element_text(size = 12),
      axis.text.y  = element_text(size = 12)
    )
  
  if (metric_name == "FDP") {
    p <- p +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 1.2) +
      coord_cartesian(ylim = c(0, fdp_y))
  }
  
  p
}

################## plotting function (for real dataset) ###############
plot_rejections_vs_q <- function(long_table, # title = "Number of Discoveries vs q",
                                 K_variants = 12, x_breaks = NULL, 
                                 method_font = 5.5) {
  
  ord <- method_order(long_table$Method)
  long_table <- long_table %>% mutate(Method = factor(Method, levels = ord))
  pal <- make_method_palette(levels(long_table$Method), K = K_variants)
  
  endpoints <- long_table %>%
    group_by(Method) %>%
    filter(q == max(q, na.rm = TRUE)) %>%
    slice_tail(n = 1) %>%
    ungroup() %>%
    mutate(label = as.character(Method))
  
  x_rng   <- range(long_table$q, na.rm = TRUE)
  x_nudge <- diff(x_rng) * 0.03
  
  p <- ggplot(long_table, aes(x = q, y = Rejections, group = Method, color = Method)) +
    geom_line(linewidth = 0.9, na.rm = TRUE) +
    geom_point(size = 1.8, na.rm = TRUE) +
    ggrepel::geom_text_repel(
      data = endpoints,
      aes(label = label),
      size = method_font,
      nudge_x = x_nudge,
      hjust = 0,
      direction = "y",
      segment.size = 0.2,
      box.padding = 0.15,
      point.padding = 0.10,
      min.segment.length = 0,
      max.overlaps = Inf,
      seed = 1,
      show.legend = FALSE
    ) +
    scale_color_manual(values = pal, breaks = ord) +
    labs(# title = title, 
      x = "q level", y = "Number of Rejections") +
    coord_cartesian(clip = "off") +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.margin = margin(5.5, 60, 5.5, 5.5),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      axis.text.x  = element_text(size = 12),
      axis.text.y  = element_text(size = 12)
    )
  
  if (!is.null(x_breaks)) {
    p <- p + scale_x_continuous(breaks = x_breaks, expand = expansion(mult = c(0.02, 0.18)))
  } else {
    p <- p + scale_x_continuous(expand = expansion(mult = c(0.02, 0.18)))
  }
  
  p
}

