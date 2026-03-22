pacman::p_load(rio, dplyr, brms, metafor, tidybayes, ggplot2, ggdist,
               patchwork, tibble, tidyr)

# ---- Step 1: Informative priors from external studies ----

dat <- import("pseudomonas_data.xlsx")

external_no_pseudo <- dat |>
  filter(study %in% c("von Dach 2020", "Molina 2022", "Yahav 2019 no Pseudomonas"))

es_no_pseudo <- escalc(
  measure = "OR",
  ai = seven_events,
  bi = seven_total - seven_events,
  ci = fourteen_events,
  di = fourteen_total - fourteen_events,
  data = external_no_pseudo
)

re_no_pseudo <- rma(yi, vi, data = es_no_pseudo, test = "knha")
treat_prior_mean <- as.numeric(re_no_pseudo$beta)
treat_prior_sd   <- as.numeric(re_no_pseudo$se)

yahav <- dat |> filter(grepl("Yahav", study))

es_yahav <- escalc(
  measure = "OR",
  ai = seven_events,
  bi = seven_total - seven_events,
  ci = fourteen_events,
  di = fourteen_total - fourteen_events,
  data = yahav
)

yahav_no_pseudo  <- es_yahav |> filter(subgroup == "No Pseudomonas")
yahav_pseudo     <- es_yahav |> filter(subgroup == "Only Pseudomonas")
interaction_mean <- as.numeric(yahav_pseudo$yi - yahav_no_pseudo$yi)
interaction_sd   <- sqrt(as.numeric(yahav_pseudo$vi + yahav_no_pseudo$vi))

# ---- Step 2: Expand BALANCE data to Bernoulli (individual-level) format ----

balance <- dat |> filter(grepl("BALANCE", study))

expand_arm <- function(events, total, treat, pseudo) {
  data.frame(
    y     = c(rep(1L, events), rep(0L, total - events)),
    treat = treat,
    pseudo = pseudo
  )
}

balance_long <- bind_rows(
  expand_arm(
    balance$fourteen_events[balance$subgroup == "No Pseudomonas"],
    balance$fourteen_total[balance$subgroup == "No Pseudomonas"],
    treat = 0L, pseudo = 0L
  ),
  expand_arm(
    balance$seven_events[balance$subgroup == "No Pseudomonas"],
    balance$seven_total[balance$subgroup == "No Pseudomonas"],
    treat = 1L, pseudo = 0L
  ),
  expand_arm(
    balance$fourteen_events[balance$subgroup == "Only Pseudomonas"],
    balance$fourteen_total[balance$subgroup == "Only Pseudomonas"],
    treat = 0L, pseudo = 1L
  ),
  expand_arm(
    balance$seven_events[balance$subgroup == "Only Pseudomonas"],
    balance$seven_total[balance$subgroup == "Only Pseudomonas"],
    treat = 1L, pseudo = 1L
  )
)

# ---- Step 3: Fit INFORMATIVE prior model ----

treat_prior_str       <- sprintf("normal(%0.4f, %0.4f)", treat_prior_mean, treat_prior_sd)
interaction_prior_str <- sprintf("normal(%0.4f, %0.4f)", interaction_mean, interaction_sd)

info_priors <- c(
  prior(normal(0, 1.5), class = "Intercept"),
  prior_string(treat_prior_str,       class = "b", coef = "treat"),
  prior(normal(0, 1.5),               class = "b", coef = "pseudo"),
  prior_string(interaction_prior_str, class = "b", coef = "treat:pseudo")
)

info_fit <- brm(
  y ~ treat * pseudo,
  data    = balance_long,
  family  = bernoulli(link = "logit"),
  prior   = info_priors,
  sample_prior = "yes",
  chains  = 4,
  iter    = 4000,
  warmup  = 2000,
  cores   = 4,
  seed    = 123,
  backend = "cmdstanr",
  control = list(adapt_delta = 0.95),
  file = "fits/primary",
  file_refit = "on_change"
)

plot(info_fit)

# ---- Step 4: Fit WEAKLY INFORMATIVE prior model ----

weak_priors <- c(
  prior(normal(0, 1.5), class = "Intercept"),
  prior(normal(0, 0.82), class = "b", coef = "treat"),
  prior(normal(0, 1.5),  class = "b", coef = "pseudo"),
  prior(normal(0, 0.71), class = "b", coef = "treat:pseudo")
)

weak_fit <- brm(
  y ~ treat * pseudo,
  data    = balance_long,
  family  = bernoulli(link = "logit"),
  prior   = weak_priors,
  sample_prior = "yes",
  chains  = 4,
  iter    = 4000,
  warmup  = 2000,
  cores   = 4,
  seed    = 123,
  backend = "cmdstanr",
  control = list(adapt_delta = 0.95),
  file = "fits/sensitivity",
  file_refit = "on_change"
)

plot(weak_fit)

# ---- Step 5: Extract draws ----

info_post <- as_draws_df(info_fit)
weak_post <- as_draws_df(weak_fit)

# Informative posterior draws
info_no_pseudo   <- info_post$b_treat
info_pseudo      <- info_post$b_treat + info_post$`b_treat:pseudo`
info_interaction <- info_post$`b_treat:pseudo`

# Informative prior draws (sampled via sample_prior = "yes")
info_prior_no_pseudo   <- info_post$prior_b_treat
info_prior_pseudo      <- info_post$prior_b_treat + info_post$`prior_b_treat:pseudo`
info_prior_interaction <- info_post$`prior_b_treat:pseudo`

# Weakly informative posterior draws
weak_no_pseudo   <- weak_post$b_treat
weak_pseudo      <- weak_post$b_treat + weak_post$`b_treat:pseudo`
weak_interaction <- weak_post$`b_treat:pseudo`

# Weakly informative prior draws
weak_prior_no_pseudo   <- weak_post$prior_b_treat
weak_prior_pseudo      <- weak_post$prior_b_treat + weak_post$`prior_b_treat:pseudo`
weak_prior_interaction <- weak_post$`prior_b_treat:pseudo`

# ---- Step 6: Compute posterior probabilities ----

# Row 1 (No Pseudo): P(OR < 1) = P(log-OR < 0)
info_prob_no_pseudo <- round(mean(info_no_pseudo < 0) * 100, 0)
weak_prob_no_pseudo <- round(mean(weak_no_pseudo < 0) * 100, 0)

# Row 2 (Pseudo): P(OR > 1) = P(log-OR > 0)
info_prob_pseudo <- round(mean(info_pseudo > 0) * 100, 0)
weak_prob_pseudo <- round(mean(weak_pseudo > 0) * 100, 0)

# Row 3 (Interaction): P(ROR > 1) = P(interaction > 0)
info_prob_interaction <- round(mean(info_interaction > 0) * 100, 0)
weak_prob_interaction <- round(mean(weak_interaction > 0) * 100, 0)

# ---- Step 7: Shared theme and scales ----

fig_theme <- theme_bw(base_size = 10) +
  theme(
    plot.title       = element_text(face = "bold", size = 10, hjust = 0.5),
    axis.text.y      = element_blank(),
    axis.ticks.y     = element_blank(),
    axis.title.x     = element_text(size = 9),
    axis.title.y     = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(color = "gray85", linewidth = 0.3),
    legend.position  = "none",
    plot.margin      = margin(4, 8, 4, 8)
  )

or_breaks  <- c(0.1, 0.2, 0.5, 1, 2, 5, 10)
ror_breaks <- c(0.1, 0.2, 0.5, 1, 2, 5, 10)

# Colors for each row
col_no_pseudo   <- "#3182bd"  # blue
col_pseudo      <- "#e6550d"  # orange
col_interaction <- "#31a354"  # green

# ---- Helper function for panels ----

make_panel <- function(post_draws, prior_draws, fill_col, threshold, direction,
                       prob_pct, x_breaks, x_labels, x_title, xlims,
                       panel_title) {

  df_post  <- tibble(x = post_draws, dist = "Posterior")
  df_prior <- tibble(x = prior_draws, dist = "Prior")

  # Posterior slab with conditional fill
  if (direction == "below") {
    p <- ggplot() +
      # Prior slab: gray outline, no fill
      stat_slab(
        data = df_prior, aes(x = x, y = 0),
        fill = NA, color = "gray20", linewidth = 0.5, alpha = 0.3,
        normalize = "panels"
      ) +
      # Posterior slab: colored fill by threshold
      stat_slab(
        data = df_post, aes(x = x, y = 0,
                            fill = after_stat(x < threshold)),
        color = NA, slab_alpha = 0.8,
        normalize = "panels"
      ) +
      scale_fill_manual(values = c("TRUE" = fill_col, "FALSE" = "gray80"))
  } else {
    p <- ggplot() +
      stat_slab(
        data = df_prior, aes(x = x, y = 0),
        fill = NA, color = "gray20", linewidth = 0.5, alpha = 0.3,
        normalize = "panels"
      ) +
      stat_slab(
        data = df_post, aes(x = x, y = 0,
                            fill = after_stat(x > threshold)),
        color = NA, slab_alpha = 0.8,
        normalize = "panels"
      ) +
      scale_fill_manual(values = c("TRUE" = fill_col, "FALSE" = "gray80"))
  }

  # Add point + interval for posterior on top
  p <- p +
    stat_pointinterval(
      data = df_post, aes(x = x, y = 0),
      point_interval = median_hdi, .width = 0.95,
      point_size = 2, interval_size_range = c(0.8, 1.2),
      color = "gray20"
    ) +
    geom_vline(xintercept = threshold, linetype = "dashed", color = "gray40") +
    annotate(
      "text",
      x = ifelse(direction == "below",
                 quantile(post_draws, 0.015),
                 quantile(post_draws, 0.985)),
      y = 0.85,
      label = paste0(prob_pct, "%"),
      size  = 3.5,
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

# ---- Step 8: Build 6 panels ----

or_xlim  <- log(c(0.1, 10))
ror_xlim <- log(c(0.1, 10))

# Row 1: No Pseudomonas — P(OR < 1), fill below 0, blue
p_info_no_pseudo <- make_panel(
  info_no_pseudo, info_prior_no_pseudo,
  fill_col = col_no_pseudo, threshold = 0, direction = "below",
  prob_pct = info_prob_no_pseudo,
  x_breaks = log(or_breaks), x_labels = or_breaks,
  x_title = "Odds Ratio (log scale)", xlims = or_xlim,
  panel_title = "No Pseudomonas"
)

p_weak_no_pseudo <- make_panel(
  weak_no_pseudo, weak_prior_no_pseudo,
  fill_col = col_no_pseudo, threshold = 0, direction = "below",
  prob_pct = weak_prob_no_pseudo,
  x_breaks = log(or_breaks), x_labels = or_breaks,
  x_title = "Odds Ratio (log scale)", xlims = or_xlim,
  panel_title = "No Pseudomonas"
)

# Row 2: Pseudomonas — P(OR > 1), fill above 0, orange
p_info_pseudo <- make_panel(
  info_pseudo, info_prior_pseudo,
  fill_col = col_pseudo, threshold = 0, direction = "above",
  prob_pct = info_prob_pseudo,
  x_breaks = log(or_breaks), x_labels = or_breaks,
  x_title = "Odds Ratio (log scale)", xlims = or_xlim,
  panel_title = "Pseudomonas"
)

p_weak_pseudo <- make_panel(
  weak_pseudo, weak_prior_pseudo,
  fill_col = col_pseudo, threshold = 0, direction = "above",
  prob_pct = weak_prob_pseudo,
  x_breaks = log(or_breaks), x_labels = or_breaks,
  x_title = "Odds Ratio (log scale)", xlims = or_xlim,
  panel_title = "Pseudomonas"
)

# Row 3: Interaction — P(ROR > 1), fill above 0, green
p_info_interaction <- make_panel(
  info_interaction, info_prior_interaction,
  fill_col = col_interaction, threshold = 0, direction = "above",
  prob_pct = info_prob_interaction,
  x_breaks = log(ror_breaks), x_labels = ror_breaks,
  x_title = "Ratio of Odds Ratios (log scale)", xlims = ror_xlim,
  panel_title = "Interaction"
)

p_weak_interaction <- make_panel(
  weak_interaction, weak_prior_interaction,
  fill_col = col_interaction, threshold = 0, direction = "above",
  prob_pct = weak_prob_interaction,
  x_breaks = log(ror_breaks), x_labels = ror_breaks,
  x_title = "Ratio of Odds Ratios (log scale)", xlims = ror_xlim,
  panel_title = "Interaction"
)

# ---- Step 9: Assemble with patchwork ----

# Column header labels
library(grid)

header_primary <- wrap_elements(
  textGrob("Primary Analysis\n(Informative Priors)",
           gp = gpar(fontsize = 12, fontface = "bold"))
)
header_sensitivity <- wrap_elements(
  textGrob("Sensitivity Analysis\n(Weakly Informative Priors)",
           gp = gpar(fontsize = 12, fontface = "bold"))
)

# Layout: header row + 3 data rows, 2 columns
figure_letter <- (header_primary | header_sensitivity) /
  (p_info_no_pseudo | p_weak_no_pseudo) /
  (p_info_pseudo    | p_weak_pseudo) /
  (p_info_interaction | p_weak_interaction) +
  plot_layout(heights = c(0.8, 3, 3, 3))

ggsave(
  "figure_letter.png",
  plot   = figure_letter,
  width  = 8,
  height = 10,
  dpi    = 300,
  bg     = "white"
)

message("Saved: figure_letter.png")

# Print probabilities for letter text
cat("\n=== Values for letter text ===\n")
cat("No Pseudo OR:", round(exp(median(info_no_pseudo)), 2), "\n")
cat("No Pseudo 95% CrI:", round(exp(quantile(info_no_pseudo, 0.025)), 2), "-",
    round(exp(quantile(info_no_pseudo, 0.975)), 2), "\n")
cat("No Pseudo P(OR<1):", info_prob_no_pseudo, "%\n")
cat("Pseudo OR:", round(exp(median(info_pseudo)), 2), "\n")
cat("Pseudo 95% CrI:", round(exp(quantile(info_pseudo, 0.025)), 2), "-",
    round(exp(quantile(info_pseudo, 0.975)), 2), "\n")
cat("Pseudo P(OR>1):", info_prob_pseudo, "%\n")
cat("Interaction P(ROR>1):", info_prob_interaction, "%\n")
