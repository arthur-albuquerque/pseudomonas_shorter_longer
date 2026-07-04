# =========================================================================
# 01_main.R
# Two-stage IPD meta-analysis: 7- vs 14-day antibiotic course and all-cause
# mortality, with Pseudomonas aeruginosa as the effect modifier.
#
# Covers:
#   - Stage 1 (DATA): per-trial effect sizes (metafor::escalc) and
#     within-trial interactions. Figure 1: data-only forest plot.
#   - Stage 2 (FITS): Bayesian fixed-effects meta-analyses (brms).
#   - Figure 2: three pooled posteriors (single column).
#   - Pooled summaries (RR scale).
#
# See 02_diagnostics.R for posterior predictive checks.
# See 03_sensitivity.R for prior sensitivity analysis.
# =========================================================================
#
# The Stage-2 Bayesian fixed-effects meta-analysis uses the exact conjugate
# normal-normal closed form (no MCMC).  See:
#   Smith TC, Spiegelhalter DJ, Thomas A. Bayesian approaches to random-
#   effects meta-analysis: a comparative study. Statistics in Medicine,
#   1995, 14(24): 2685-2699. DOI: 10.1002/sim.4780142408

pacman::p_load(rio, dplyr, purrr, metafor, distributional,
               ggplot2, ggdist, patchwork, tibble, tidyr, MetBrewer, grid,
               ragg, ggtext, systemfonts)

dir.create("fits", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# ---- Publication typography ----
# One sans-serif family for both Figure 1 (base R) and Figure 2 / Concept A
# (ggplot2). Helvetica Neue is the de-facto standard for clinical-journal
# figures; we deliberately use sans (not the Tufte serif default) to match
# medical-journal convention. Rendered through the ragg device so font metrics
# are identical across base-R and ggplot output.
font_main <- "Helvetica Neue"

# Point sizes at the as-drawn figure dimensions (figures are intentionally
# large; these stay legible when a journal scales them to ~7 in column width).
fs_tag      <- 18    # Concept A corner tags (A / B)
fs_colhdr   <- 15    # Concept A column headers
fs_callout  <- 15    # bold P(%) headline number
fs_title    <- 12    # panel titles
fs_distlab  <- 11    # "Posterior" / "Prior"
fs_axis_ttl <- 10.5  # axis titles
fs_annot    <- 9.5   # "Primary Analysis"
fs_axis_txt <- 9     # numeric tick labels

# ---- Step 1: Load data, add SHORTEN2, normalise study labels ----
# The trial matrix -- 5 RCTs, two strata (Pseudomonas / Other Gram-negative):
#   Yahav 2019, BALANCE        -> BOTH strata (the only within-trial interaction info)
#   von Dach 2020, Molina 2022 -> Other Gram-negative stratum only
#   SHORTEN2                   -> Pseudomonas stratum only
# pseudomonas_data.xlsx has one row per study-subgroup: study, subgroup
# ("No Pseudomonas" / "Only Pseudomonas"), and 7-/14-day event & total counts.
dat <- import("pseudomonas_data.xlsx")

# SHORTEN2 contributes a Pseudomonas-only stratum (not in the xlsx).
shorten2 <- tibble(
  study = "SHORTEN2", subgroup = "Only Pseudomonas",
  seven_events = 6L, seven_total = 149L,
  fourteen_events = 7L, fourteen_total = 154L
)

dat <- bind_rows(dat, shorten2) |>
  mutate(
    study = case_when(
      grepl("Yahav",   study) ~ "Yahav 2019",
      grepl("BALANCE", study) ~ "BALANCE",
      TRUE ~ study
    ),
    # sublab: display label for Figure 1 only; the analysis keys off `subgroup`.
    sublab = ifelse(subgroup == "Only Pseudomonas", "Pseudomonas", "Other Gram-negative")
  )

# ---- Step 2: Stage-1 effect sizes (escalc) + within-trial interactions ----
# treat arm = 7-day (seven_*); control arm = 14-day (fourteen_*).
# mk_es(): Stage-1 effect sizes. metafor::escalc builds, for every study-subgroup
# row, the log effect size yi and sampling variance vi from the 2x2 table.
mk_es <- function(measure) {
  escalc(measure = measure,
         ai = seven_events, bi = seven_total - seven_events,
         ci = fourteen_events, di = fourteen_total - fourteen_events,
         data = dat)
}

# build_sets(): assemble the three Stage-1 inputs for one scale.
#   other = log effect sizes for the "No Pseudomonas" (Other Gram-negative) rows
#   pseud = log effect sizes for the "Only Pseudomonas" rows
#   inter = the WITHIN-TRIAL interaction (Yahav 2019, BALANCE only):
#               yi = (log effect | Pseudomonas) - (log effect | Other GN),
#               vi = vi_pseud + vi_other (independent strata, variances add).
build_sets <- function(measure) {
  es <- mk_es(measure)
  es$sei <- sqrt(es$vi)
  other <- es |> filter(subgroup == "No Pseudomonas")   |> select(study, yi, vi, sei)
  pseud <- es |> filter(subgroup == "Only Pseudomonas") |> select(study, yi, vi, sei)
  both  <- intersect(es$study[es$subgroup == "No Pseudomonas"],
                     es$study[es$subgroup == "Only Pseudomonas"])
  inter <- lapply(both, function(s) {
    p <- es[es$study == s & es$subgroup == "Only Pseudomonas", ]
    o <- es[es$study == s & es$subgroup == "No Pseudomonas", ]
    data.frame(study = s, yi = p$yi - o$yi, vi = p$vi + o$vi)
  }) |> bind_rows()
  inter$sei <- sqrt(inter$vi)
  list(other = as.data.frame(other),
       pseud = as.data.frame(pseud),
       inter = as.data.frame(inter))
}

# ---- Step 3: Figure 1 — two-panel DATA forest (metafor::forest) ----
# Left: per-trial, per-subgroup effect sizes grouped by trial, with 7d/14d
# event counts. Right: within-trial interaction per trial (circles, Fisher et
# al. style), "Insufficient data" where both strata are absent. No pooling.
draw_fig1 <- function() {
  es <- escalc(measure = "RR",
               ai = seven_events, bi = seven_total - seven_events,
               ci = fourteen_events, di = fourteen_total - fourteen_events,
               data = dat)
  es$sublab <- ifelse(es$subgroup == "Only Pseudomonas", "Pseudomonas", "Other Gram-negative")

  trial_order <- c("Yahav 2019", "von Dach 2020", "Molina 2022", "BALANCE", "SHORTEN2")
  sub_order   <- c("Pseudomonas", "Other Gram-negative")

  data_rows <- list(); header_y <- c(); trial_block_y <- c(); ycur <- 100
  for (tr in trial_order) {
    sub_here <- es[es$study == tr, ]
    sub_here <- sub_here[order(match(sub_here$sublab, sub_order)), ]
    header_y[tr] <- ycur; ycur <- ycur - 1; rys <- c()
    for (i in seq_len(nrow(sub_here))) {
      data_rows[[length(data_rows) + 1]] <- cbind(sub_here[i, ], y = ycur)
      rys <- c(rys, ycur); ycur <- ycur - 1
    }
    trial_block_y[tr] <- mean(rys); ycur <- ycur - 1
  }
  D <- do.call(rbind, data_rows)
  shift <- 2 - min(D$y)
  D$y <- D$y + shift; header_y <- header_y + shift; trial_block_y <- trial_block_y + shift
  ylim <- c(0, max(header_y) + 2.5); rows_vec <- D$y
  D$lab7  <- sprintf("%d/%d", D$seven_events, D$seven_total)
  D$lab14 <- sprintf("%d/%d", D$fourteen_events, D$fourteen_total)
  ci <- 1.96 * sqrt(D$vi)
  D$eff <- sprintf("%.2f (%.2f, %.2f)", exp(D$yi), exp(D$yi - ci), exp(D$yi + ci))

  both <- intersect(es$study[es$subgroup == "No Pseudomonas"],
                    es$study[es$subgroup == "Only Pseudomonas"])
  int_df <- data.frame()
  for (tr in trial_order) {
    if (tr %in% both) {
      p <- es[es$study == tr & es$subgroup == "Only Pseudomonas", ]
      o <- es[es$study == tr & es$subgroup == "No Pseudomonas", ]
      yi <- p$yi - o$yi; vi <- p$vi + o$vi; ic <- 1.96 * sqrt(vi)
      int_df <- rbind(int_df, data.frame(study = tr, yi = yi, vi = vi,
        y = trial_block_y[tr], est = TRUE,
        lab = sprintf("%.2f (%.2f, %.2f)", exp(yi), exp(yi - ic), exp(yi + ic))))
    } else {
      int_df <- rbind(int_df, data.frame(study = tr, yi = NA, vi = NA,
        y = trial_block_y[tr], est = FALSE, lab = NA))
    }
  }

  efflab <- "Risk Ratio"
  rorlab <- "Ratio of Risk Ratios"

  par(mar = c(5, 4, 2, 2), family = font_main)
  alim <- log(c(0.2, 5)); at_main <- log(c(0.2, 0.5, 1, 2, 5))
  xlim <- c(log(0.00002), log(110000))
  ilab.x  <- c(log(0.0025), log(0.02))
  annot.x <- log(15)
  rx0 <- log(80); rx1 <- log(1600)
  rannot.x <- log(13000)

  fres <- forest(D$yi, D$vi, slab = NA, rows = rows_vec, ylim = ylim,
         xlim = xlim, alim = alim, at = at_main, atransf = exp, refline = NA,
         annotate = FALSE, pch = 15, psize = 1.1, header = FALSE,
         xlab = paste0(efflab, " (7d vs 14d)"))
  segments(log(1), ylim[1], log(1), max(header_y) + 0.3, lty = 2, col = "gray50")

  text(xlim[1] + 1.0, rows_vec, D$sublab, pos = 4, cex = 0.78,
       font = ifelse(D$sublab == "Pseudomonas", 3, 1))
  text(xlim[1], header_y, trial_order, pos = 4, font = 2, cex = 0.9)
  htxt <- max(rows_vec) + 2.2
  text(xlim[1], htxt, "Study / Subgroup", pos = 4, font = 2, cex = 0.72)
  text(ilab.x[1], rows_vec, D$lab7,  cex = 0.72)
  text(ilab.x[2], rows_vec, D$lab14, cex = 0.72)
  text(ilab.x, htxt, c("7-day\nevents/total", "14-day\nevents/total"), cex = 0.72, font = 2)
  text(annot.x, htxt, paste0(efflab, "\n(95% CI)"), cex = 0.72, font = 2)
  text(annot.x, rows_vec, D$eff, cex = 0.72)

  ror_lo <- log(0.2); ror_hi <- log(5)
  mapx <- function(v) rx0 + (v - ror_lo) / (ror_hi - ror_lo) * (rx1 - rx0)
  ref_int <- mapx(log(1))
  segments(ref_int, ylim[1], ref_int, max(header_y) + 0.3, lty = 2, col = "gray50")
  ats <- c(0.2, 0.5, 1, 2, 5)
  axis(side = 1, at = mapx(log(ats)), labels = ats, cex.axis = fres$cex.axis, gap.axis = -1, line = 0)
  text(mean(c(rx0, rx1)), htxt, "Subgroup Difference\n(Interaction)", font = 2, cex = 0.78)
  mtext("Ratio of Risk Ratios", side = 1, at = mean(c(rx0, rx1)), line = 2.4, cex = fres$cex.lab)
  text(rannot.x, htxt, paste0(rorlab, "\n(95% CI)"), cex = 0.72, font = 2)
  for (i in seq_len(nrow(int_df))) {
    yy <- int_df$y[i]
    if (int_df$est[i]) {
      icw <- 1.96 * sqrt(int_df$vi[i])
      raw_lo <- mapx(int_df$yi[i] - icw); raw_hi <- mapx(int_df$yi[i] + icw)
      lo <- max(raw_lo, rx0); hi <- min(raw_hi, rx1)
      segments(lo, yy, hi, yy, col = "black", lwd = 1.4)
      aw <- (rx1 - rx0) * 0.03
      if (raw_hi > rx1) arrows(rx1 - aw, yy, rx1, yy, length = 0.05, angle = 30, code = 2, lwd = 1.4)
      if (raw_lo < rx0) arrows(rx0 + aw, yy, rx0, yy, length = 0.05, angle = 30, code = 2, lwd = 1.4)
      points(mapx(int_df$yi[i]), yy, pch = 21, bg = "white", col = "black", cex = 1.4)
      text(rannot.x, yy, int_df$lab[i], cex = 0.72)
    } else {
      text(mean(c(rx0, rx1)), yy, "Insufficient data", cex = 0.72, font = 3, col = "gray40")
    }
  }
}

ragg::agg_png("figures/figure1_forest_rr.png", width = 10.25, height = 6.5,
    units = "in", res = 300, background = "white"); draw_fig1(); dev.off()
ragg::agg_tiff("figures/figure1_forest_rr.tiff", width = 10.25, height = 6.5,
     units = "in", res = 300, background = "white"); draw_fig1(); dev.off()
message("Saved: figures/figure1_forest_rr.{png,tiff}")

# ---- Step 4: Stage-2 Bayesian fixed-effects meta-analyses (closed form) ----
# Model: yi | se(sei) ~ 1 with Intercept (pooled effect mu) ~ N(0, tau).
# Because se(sei) fixes the residual SD at the known Stage-1 SE, the only
# parameter is mu and its posterior is EXACTLY Gaussian (conjugate
# normal-normal). No MCMC is needed; densities, ribbons, and probability strips
# in Figure 2 all flow from this one model.
#
# Prior widths (tau) on the pooled effect:
#   Subgroup models (Other GN, Pseudomonas): tau = 0.82;
#     exp(+/-1.96*0.82) ~ (0.2, 5) covers the plausible RR range.
#   Interaction model: tau = 0.35, tighter; exp(+/-1.96*0.35) ~ (0.5, 2.0).
#     Effect-modification (RoRR) effects are expected smaller than main
#     effects, and only 2 trials inform the interaction, so it is regularised
#     more strongly toward the null.
tau_fixed <- 0.82   # subgroup models (Other GN, Pseudomonas)
tau_inter <- 0.35   # interaction model

# conj(): exact conjugate normal-normal posterior for the pooled effect mu.
# Given study log-effects yi with SEs sei and prior mu ~ N(0, tau):
#   A = sum(1/sei^2); B = sum(yi/sei^2)
#   v = 1/(1/tau^2 + A)  (posterior variance);  m = v*B  (posterior mean)
#   mu | data ~ Normal(m, sqrt(v))
conj <- function(d, tau) {
  w <- 1 / d$sei^2
  A <- sum(w); B <- sum(d$yi * w)
  v <- 1 / (1 / tau^2 + A)
  list(m = v * B, sd = sqrt(v))
}

S_rr <- build_sets("RR")

post_other <- conj(S_rr$other, tau_fixed)
post_pseud <- conj(S_rr$pseud, tau_fixed)
post_inter <- conj(S_rr$inter, tau_inter)

# ---- Step 5: Primary-analysis directional posterior probabilities ----
# Exact directional tail probability (%) of the Gaussian posterior, UNROUNDED
# (the strips reuse this for a smooth curve; the % callouts round it). Same model
# as everything else in Figure 2.
prob_dir <- function(post, direction) {
  z <- post$m / post$sd
  (if (direction == "above") pnorm(z) else pnorm(-z)) * 100
}
pp_other <- round(prob_dir(post_other, "below"))   # P(RR < 1)
pp_pseud <- round(prob_dir(post_pseud, "above"))   # P(RR > 1)
pp_inter <- round(prob_dir(post_inter, "above"))   # P(RoRR > 1)

# ---- Step 6: Shared theme + make_panel helper ----
fig_theme <- theme_bw(base_size = fs_axis_txt, base_family = font_main) +
  theme(
    plot.title       = element_markdown(face = "bold", size = fs_title, hjust = 0.5),
    axis.text.x      = element_text(size = fs_axis_txt, color = "gray20"),
    axis.ticks.y     = element_blank(),
    axis.title.x     = element_text(size = fs_axis_ttl),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(color = "gray85", linewidth = 0.3),
    legend.position  = "none",
    plot.margin      = margin(5, 3, 5, 3)
  )

rr_breaks   <- c(0.5, 0.7, 1, 1.5, 2)
rorr_breaks <- c(0.5, 0.7, 1, 1.5, 2)

colors          <- met.brewer(name = "Isfahan1", n = 8, type = "discrete")
col_no_pseudo   <- colors[2]
col_pseudo      <- colors[8]
col_interaction <- colors[5]

make_panel <- function(post, fill_col, threshold, direction,
                        prob_pct, x_breaks, x_labels, x_title, xlims,
                        panel_title) {
  post_dist <- dist_normal(post$m, post$sd)

  if (direction == "below") {
    p <- ggplot() +
      stat_slab(
        aes(xdist = post_dist, y = 0,
            fill = after_stat(x < threshold)),
        color = NA, slab_alpha = 0.8,
        normalize = "panels"
      ) +
      scale_fill_manual(values = c("TRUE" = fill_col, "FALSE" = "gray80"))
  } else {
    p <- ggplot() +
      stat_slab(
        aes(xdist = post_dist, y = 0,
            fill = after_stat(x > threshold)),
        color = NA, slab_alpha = 0.8,
        normalize = "panels"
      ) +
      scale_fill_manual(values = c("TRUE" = fill_col, "FALSE" = "gray80"))
  }

  p <- p +
    stat_pointinterval(
      aes(xdist = post_dist, y = 0),
      .width = 0.95,
      point_size = 2, interval_size_range = c(0.8, 1.2),
      color = "gray20"
    ) +
    geom_vline(xintercept = threshold, linetype = "dashed", color = "gray40") +
    annotate(
      "text",
      x = ifelse(direction == "below",
                   qnorm(0.025, post$m, post$sd),
                   qnorm(0.9,   post$m, post$sd)),
      y = 0.85,
      label = paste0(prob_pct, "%"),
      size  = fs_callout / .pt,
      family = font_main,
      fontface = "bold",
      hjust = ifelse(direction == "below", 1, 0),
      color = fill_col
    ) +
    scale_x_continuous(
      breaks = x_breaks,
      labels = x_labels,
      name   = x_title
    ) +
    scale_y_continuous(labels = NULL, name = NULL, expand = expansion(mult = c(0, 0.05))) +
    coord_cartesian(xlim = xlims, ylim = c(-0.05, NA)) +
    labs(title = panel_title) +
    fig_theme

  return(p)
}

# ---- Step 7: shared x-limits for the posterior panels (log RR scale) ----
# threshold = 0 is the log-scale null (ratio = 1); `direction` sets which tail
# is shaded: Pseudomonas/Interaction -> "above", Other GN -> "below".
# The single-column posteriors are not saved on their own; they form column A
# of the combined Figure 2 assembled in Step 9.
xlim2 <- log(c(0.4, 2.5))

# ---- Step 8: Print pooled summaries ----
summ <- function(label, post, prob, dir_lt) {
  z <- qnorm(0.975)
  cat(sprintf("\n%s: %.2f (95%% CrI %.2f - %.2f), P(%s)=%d%%",
      label, exp(post$m), exp(post$m - z * post$sd),
      exp(post$m + z * post$sd), ifelse(dir_lt, "<1", ">1"), prob))
}
cat("\n=== RR scale ===")
summ("Other GN RR (7 vs 14d)",    post_other, pp_other, TRUE)
summ("Pseudomonas RR (7 vs 14d)", post_pseud, pp_pseud, FALSE)
summ("Interaction RoRR",          post_inter, pp_inter, FALSE)
cat("\n")

# =========================================================================
# Step 9: Figure 2 — combined primary-analysis + prior-sensitivity composite.
# Two labelled columns: A = primary-analysis posteriors (the Step 6 panels,
# priors not drawn so the x-axis can zoom onto the posteriors), B = a single
# posterior 95%-CrI ribbon per model traced across the intercept-prior SD
# grid. The prior funnel is omitted because it is identical across models
# (N(0, sigma)) and added only visual redundancy; the primary-analysis
# sigma is marked with a dotted horizontal line and a directional-% label.
#
# Every element is the exact conjugate normal-normal posterior (Step 4 conj()),
# evaluated across the intercept-prior width grid below. No MCMC, no caching.
# =========================================================================

# ---- 9a: prior-width sensitivity grid (intercept prior SD) ----
# The same six prior widths sweep both Column-B elements (ribbons + strips).
# Vector order is irrelevant: every consumer orders points by sigma.
sens_sigmas <- seq(0.095, 1.18, 0.01)

# ---- 9b: rebuild the three LEFT posterior panels (same calls as Step 7) ----
p_pseud <- make_panel(
  post_pseud, fill_col = col_pseudo, threshold = 0,
  direction = "above", prob_pct = pp_pseud,
  x_breaks = log(rr_breaks), x_labels = rr_breaks,
  x_title = "Risk Ratio", xlims = xlim2, panel_title = "*Pseudomonas* Subgroup")
p_other <- make_panel(
  post_other, fill_col = col_no_pseudo, threshold = 0,
  direction = "below", prob_pct = pp_other,
  x_breaks = log(rr_breaks), x_labels = rr_breaks,
  x_title = "Risk Ratio", xlims = xlim2, panel_title = "Other Gram-negative Subgroup")
p_int <- make_panel(
  post_inter, fill_col = col_interaction, threshold = 0,
  direction = "above", prob_pct = pp_inter,
  x_breaks = log(rorr_breaks), x_labels = rorr_breaks,
  x_title = "Ratio of Risk Ratios", xlims = xlim2,
  panel_title = "Subgroup Difference (Interaction)")

# =========================================================================
# Column B ribbons — prior-vs-posterior stability ribbon.
# x = effect (shared with left column); y = prior width sigma, axis left
# UNLABELED — the widening light-gray prior funnel itself encodes prior width.
# Posterior 95% HDI ribbon (colored) overlaid on the prior 95% HDI ribbon
# (light gray), each with a median line; the primary analysis is highlighted with
# a ggdist point-interval on top of the ribbon (matching the left-panel interval)
# + a direct "Primary Analysis" label. Straight segments between fitted sigmas.
# Uses fig_theme directly (y-axis already blanked there).
# =========================================================================

# Per-sigma 95% interval + median, in closed form. `which` picks the posterior
# (conjugate N(m, sd) at each prior width) or the prior (N(0, sigma)). For a
# Gaussian the central 95% interval equals the 95% HDI.
sens_bands <- function(d, sigmas, which = "post") {
  z <- qnorm(0.975)
  bind_rows(lapply(sigmas, function(s) {
    if (which == "post") {
      p <- conj(d, s)
      tibble(sigma = s, lo95 = p$m - z * p$sd, m = p$m, hi95 = p$m + z * p$sd)
    } else {
      tibble(sigma = s, lo95 = -z * s, m = 0, hi95 = z * s)
    }
  }))
}
# Order a (sigma -> value) series by sigma (straight-line path; no spline).
sm <- function(x, y) { o <- order(x); tibble(sigma = x[o], val = y[o]) }
# Closed polygon for a horizontal interval band (straight edges).
ribbon_poly2 <- function(band, lo, hi) {
  a <- sm(band$sigma, band[[lo]]); b <- sm(band$sigma, band[[hi]])
  tibble(x = c(a$val, rev(b$val)), y = c(a$sigma, rev(b$sigma)))
}

make_sens_ribbon <- function(d, fill_col, primary_sigma, direction, prim_prob,
                             x_title, panel_title) {
  post_band <- sens_bands(d, sens_sigmas, "post")
  post_med  <- sm(post_band$sigma, post_band$m)
  prim      <- post_band[which.min(abs(post_band$sigma - primary_sigma)), ]
  prim_post <- conj(d, primary_sigma)

  clampx    <- function(x) pmin(pmax(x, log(0.26)), log(7.8))
  lab_off   <- 0.10
  post_side <- if (prim$m >= 0) 1 else -1
  top_band  <- post_band[which.max(post_band$sigma), ]
  x_postlab <- if (post_side > 0) clampx(top_band$hi95 + lab_off) else clampx(top_band$lo95 - lab_off)
  hj_post   <- if (post_side > 0) 0 else 1
  x_primlab <- clampx(prim$hi95 + lab_off)
  x_problab <- clampx(prim$lo95 - lab_off)

  ggplot() +
    # posterior 95% CrI ribbon across prior SDs + posterior median
    geom_polygon(data = ribbon_poly2(post_band, "lo95", "hi95"),
                 aes(x, y), fill = fill_col, alpha = 0.40) +
    geom_line(data = post_med, aes(val, sigma), orientation = "y",
              color = fill_col, linewidth = 0.9) +
    # no-effect reference
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.4) +
    # primary-analysis sigma reference
    geom_hline(yintercept = primary_sigma, linetype = "dotted",
               color = "gray50", linewidth = 0.4) +
    # primary analysis highlighted on top of the ribbon
    stat_pointinterval(
      aes(xdist = dist_normal(prim_post$m, prim_post$sd), y = primary_sigma),
      .width = 0.95,
      point_size = 2, interval_size_range = c(0.8, 1.2),
      color = "gray20"
    ) +
    annotate("text", x = x_primlab, y = primary_sigma, label = "Primary Analysis",
             hjust = 0, vjust = 0.5, size = fs_annot / .pt, lineheight = 0.9,
             family = font_main, fontface = "bold", color = "black") +
    # directional posterior probability % at the primary marker (ties B to A)
    annotate("text", x = x_problab, y = primary_sigma,
             label = paste0(round(prim_prob), "%"),
             hjust = 1, vjust = 0.5, size = fs_callout / .pt,
             family = font_main, fontface = "bold", color = fill_col) +
    scale_x_continuous(breaks = log(c(0.25, 0.5, 1, 2, 5)),
                      labels = c(0.25, 0.5, 1, 2, 5), name = x_title) +
    scale_y_continuous(
      breaks = c(0.1, 0.35, 0.82, 1.18),
      labels = c(0.1, 0.35, 0.82, 1.18),
      name   = "Prior Standard Deviation"
    ) +
    coord_cartesian(xlim = log(c(0.25, 6)),
                    ylim = c(min(sens_sigmas), max(sens_sigmas))) +
    labs(title = panel_title) +
    fig_theme +
    theme(axis.text.y  = element_text(size = fs_axis_txt, color = "gray20"),
          axis.title.y = element_text(size = fs_axis_ttl))
}

sA_pseud <- make_sens_ribbon(S_rr$pseud, col_pseudo,      tau_fixed, "above", prob_dir(post_pseud, "above"),
                             "Risk Ratio",           "*Pseudomonas* Subgroup")
sA_other <- make_sens_ribbon(S_rr$other, col_no_pseudo,   tau_fixed, "below", prob_dir(post_other, "below"),
                             "Risk Ratio",           "Other Gram-negative Subgroup")
sA_inter <- make_sens_ribbon(S_rr$inter, col_interaction, tau_inter, "above", prob_dir(post_inter, "above"),
                             "Ratio of Risk Ratios", "Subgroup Difference (Interaction)")

# ---- 9c: two labelled columns (A = Primary Analysis, B = Prior Sensitivity) ----
# Each column gets a large centred header; wrap_elements collapses each column
# into a single taggable unit so plot_annotation(tag_levels = "A") prints just
# "A" (left) and "B" (right) in the corners.
hdr_left  <- wrap_elements(grid::textGrob(
  "Primary Analysis",
  gp = grid::gpar(fontsize = fs_colhdr, fontface = "bold", fontfamily = font_main)))
hdr_right <- wrap_elements(grid::textGrob(
  "Prior Sensitivity Analyses",
  gp = grid::gpar(fontsize = fs_colhdr, fontface = "bold", fontfamily = font_main)))

col_left  <- hdr_left  / wrap_elements(p_pseud / p_other / p_int) +
  plot_layout(heights = c(0.045, 1))
col_right <- hdr_right / wrap_elements(
  wrap_plots(sA_pseud, sA_other, sA_inter, ncol = 1)
) +
  plot_layout(heights = c(0.045, 1))

figA <- (wrap_elements(col_left) | wrap_elements(col_right)) +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(tag_levels = "A",
                  theme = theme(plot.tag = element_text(size = fs_tag, face = "bold",
                                                        family = font_main),
                                plot.margin = margin(2, 1, 1, 1)))

ragg::agg_png("figures/figure2.png", width = 12.5, height = 8,
    units = "in", res = 300, background = "white"); print(figA); dev.off()
ragg::agg_tiff("figures/figure2.tiff", width = 12.5, height = 8,
     units = "in", res = 300, background = "white"); print(figA); dev.off()
message("Saved: figures/figure2.{png,tiff}")
