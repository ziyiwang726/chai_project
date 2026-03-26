# ---------- color helpers ----------
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
    DESeq2      = "#E377C2",
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

  out <- c(core, sort_family("OrderShapeEM"), sort_family("IHW"), "DESeq2", "BH")
  out <- unique(out)
  out[out %in% methods]
}


plot_rejections_vs_q <- function(long_table, title = "Number of Discoveries vs q",
                                 K_variants = 12, x_breaks = NULL) {

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
    labs(title = title, x = "q level", y = "Number of Rejections") +
    coord_cartesian(clip = "off") +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.margin = margin(5.5, 60, 5.5, 5.5)
    )

  if (!is.null(x_breaks)) {
    p <- p + scale_x_continuous(breaks = x_breaks, expand = expansion(mult = c(0.02, 0.18)))
  } else {
    p <- p + scale_x_continuous(expand = expansion(mult = c(0.02, 0.18)))
  }

  p
}
