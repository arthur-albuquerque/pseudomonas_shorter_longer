pacman::p_load(rio, dplyr, brms, metafor, tidybayes, ggplot2, ggdist,
               patchwork, tibble, tidyr, MetBrewer, gt, marginaleffects)


# Data

dat <- import("pseudomonas_data.xlsx")

study_table <- dat |>
  select(study, subgroup, seven_events, seven_total, fourteen_events, fourteen_total) |>
  mutate(
    subgroup = ifelse(subgroup == "No Pseudomonas", "Other Gram-negative bacteria", "Pseudomonas aeruginosa"),
    study = case_when(
      study == "Yahav 2019 no Pseudomonas" ~ "Yahav 2019",
      study == "Yahav 2019 Pseudomonas Only" ~ "Yahav 2019",
      study == "BALANCE no Pseudomonas" ~ "BALANCE",
      study == "BALANCE Pseudomonas only" ~ "BALANCE",
      TRUE ~ study
    ),
    seven_risk = seven_events / seven_total,
    fourteen_risk = fourteen_events / fourteen_total
  ) |>
  select(study, subgroup, seven_events, seven_total, seven_risk, fourteen_events, fourteen_total, fourteen_risk) |>
  arrange(desc(subgroup), study) |>
  gt(groupname_col = "subgroup") |>
  cols_label(
    study = "Study",
    seven_events = "Events",
    seven_total = "Total",
    seven_risk = "Event Risk",
    fourteen_events = "Events",
    fourteen_total = "Total",
    fourteen_risk = "Event Risk"
  ) |>
  fmt_percent(columns = c(seven_risk, fourteen_risk), decimals = 1) |>
  tab_spanner(label = "7 days", columns = c(seven_events, seven_total, seven_risk)) |>
  tab_spanner(label = "14 days", columns = c(fourteen_events, fourteen_total, fourteen_risk)) |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_row_groups()
  )

if (!dir.exists("figures")) dir.create("figures")
gtsave(study_table, filename = "figures/table_studies.html")

# Total events and patients per arm, aggregated by subgroup
subgroup_arm_totals <- dat |>
  group_by(subgroup) |>
  summarise(
    `7-day_events`     = sum(seven_events),
    `7-day_patients`   = sum(seven_total),
    `7-day_event_rate` = `7-day_events` / `7-day_patients` * 100,
    `14-day_events`     = sum(fourteen_events),
    `14-day_patients`   = sum(fourteen_total),
    `14-day_event_rate` = `14-day_events` / `14-day_patients` * 100,
    .groups = "drop"
  )

subgroup_arm_totals

# ---- Step 1: Historical priors from external studies ----

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

# Sanity check: verify expansion matches original counts
check <- balance_long |>
  group_by(pseudo) |>
  summarise(n = n(), events = sum(y), .groups = "drop")

check

orig <- balance |>
  group_by(subgroup) |>
  summarise(
    n      = sum(seven_total + fourteen_total),
    events = sum(seven_events + fourteen_events),
    .groups = "drop"
  ) |>
  mutate(pseudo = ifelse(subgroup == "Only Pseudomonas", 1L, 0L))

orig

stopifnot(
  all.equal(check$n[check$pseudo == 0],      orig$n[orig$pseudo == 0]),
  all.equal(check$n[check$pseudo == 1],      orig$n[orig$pseudo == 1]),
  all.equal(check$events[check$pseudo == 0], orig$events[orig$pseudo == 0]),
  all.equal(check$events[check$pseudo == 1], orig$events[orig$pseudo == 1])
)

# ---- Step 3: Fit HISTORICAL prior model ----

treat_prior_str       <- sprintf("normal(%0.4f, %0.4f)", treat_prior_mean, treat_prior_sd)
interaction_prior_str <- sprintf("normal(%0.4f, %0.4f)", interaction_mean, interaction_sd)

info_priors <- c(
  prior(normal(0, 1.5),               class = "b", coef = "Intercept"),
  prior_string(treat_prior_str,       class = "b", coef = "treat"),
  prior(normal(0, 1.5),               class = "b", coef = "pseudo"),
  prior_string(interaction_prior_str, class = "b", coef = "treat:pseudo")
)

info_fit <- brm(
  y ~ 0 + Intercept + treat * pseudo,
  data    = balance_long,
  family  = bernoulli(link = "logit"),
  prior   = info_priors,
  sample_prior = "yes",
  chains  = 4,
  iter    = 10000,
  warmup  = 5000,
  cores   = 4,
  seed    = 123,
  backend = "cmdstanr",
  control = list(adapt_delta = 0.95),
  file = "fits/primary",
  file_refit = "on_change"
)


# ---- Step 4: Fit WEAKLY INFORMATIVE prior model ----

weak_priors <- c(
  prior(normal(0, 1.5),  class = "b", coef = "Intercept"),
  prior(normal(0, 0.82), class = "b", coef = "treat"),
  prior(normal(0, 1.5),  class = "b", coef = "pseudo"),
  prior(normal(0, 0.71), class = "b", coef = "treat:pseudo")
)

weak_fit <- brm(
  y ~ 0 + Intercept + treat * pseudo,
  data    = balance_long,
  family  = bernoulli(link = "logit"),
  prior   = weak_priors,
  sample_prior = "yes",
  chains  = 4,
  iter    = 10000,
  warmup  = 5000,
  cores   = 4,
  seed    = 123,
  backend = "cmdstanr",
  control = list(adapt_delta = 0.95),
  file = "fits/sensitivity",
  file_refit = "on_change"
)


# ---- Step 4b: Fit SKEPTICAL prior model ----

skeptic_priors <- c(
  prior(normal(0, 1.5),               class = "b", coef = "Intercept"),
  prior(normal(0, 0.71),              class = "b", coef = "treat"),
  prior(normal(0, 1.5),               class = "b", coef = "pseudo"),
  prior(normal(0, 0.36),              class = "b", coef = "treat:pseudo")
)

skeptic_fit <- brm(
  y ~ 0 + Intercept + treat * pseudo,
  data    = balance_long,
  family  = bernoulli(link = "logit"),
  prior   = skeptic_priors,
  sample_prior = "yes",
  chains  = 4,
  iter    = 10000,
  warmup  = 5000,
  cores   = 4,
  seed    = 123,
  backend = "cmdstanr",
  control = list(adapt_delta = 0.95),
  file = "fits/skeptic",
  file_refit = "on_change"
)


# ---- Step 5: Extract draws ----

info_post    <- as_draws_df(info_fit)
weak_post    <- as_draws_df(weak_fit)
skeptic_post <- as_draws_df(skeptic_fit)

# HISTORICAL posterior draws
info_no_pseudo   <- info_post$b_treat
info_pseudo      <- info_post$b_treat + info_post$`b_treat:pseudo`
info_interaction <- info_post$`b_treat:pseudo`

# HISTORICAL prior draws (sampled via sample_prior = "yes")
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

# Skeptical posterior draws
skeptic_no_pseudo   <- skeptic_post$b_treat
skeptic_pseudo      <- skeptic_post$b_treat + skeptic_post$`b_treat:pseudo`
skeptic_interaction <- skeptic_post$`b_treat:pseudo`

# Skeptical prior draws
skeptic_prior_no_pseudo   <- skeptic_post$prior_b_treat
skeptic_prior_pseudo      <- skeptic_post$prior_b_treat + skeptic_post$`prior_b_treat:pseudo`
skeptic_prior_interaction <- skeptic_post$`prior_b_treat:pseudo`

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

# Skeptical probabilities
skeptic_prob_no_pseudo   <- round(mean(skeptic_no_pseudo < 0) * 100, 0)
skeptic_prob_pseudo      <- round(mean(skeptic_pseudo > 0) * 100, 0)
skeptic_prob_interaction <- round(mean(skeptic_interaction > 0) * 100, 0)

# ---- Step 7: Shared theme and scales ----

fig_theme <- theme_bw(base_size = 10) +
  theme(
    plot.title       = element_text(face = "bold", size = 13, hjust = 0.5),
    axis.text.y      = element_blank(),
    axis.ticks.y     = element_blank(),
    axis.title.x     = element_text(size = 12),
    axis.title.y     = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(color = "gray85", linewidth = 0.3),
    legend.position  = "none",
    plot.margin      = margin(4, 8, 4, 8)
  )

or_breaks  <- c(0.1, 0.2, 0.3, 0.5, 0.7, 1, 1.5, 2, 3, 5, 10)
ror_breaks  <- c(0.1, 0.2, 0.3, 0.5, 0.7, 1, 1.5, 2, 3, 5, 10)

colors = met.brewer(name="Isfahan1", n=8, type="discrete")

col_no_pseudo   <- colors[2]
col_pseudo      <- colors[8]
col_interaction <- colors[5]


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
                 quantile(post_draws, 0.025),
                 quantile(post_draws, 0.9)),
      y = 0.85,
      label = paste0(prob_pct, "%"),
      size  = 8,
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

or_xlim  <- log(c(0.2, 5))
ror_xlim <- log(c(0.2, 5))

# Row 1: No Pseudomonas — P(OR < 1), fill below 0, blue
p_info_no_pseudo <- make_panel(
  info_no_pseudo, info_prior_no_pseudo,
  fill_col = col_no_pseudo, threshold = 0, direction = "below",
  prob_pct = info_prob_no_pseudo,
  x_breaks = log(or_breaks), x_labels = or_breaks,
  x_title = "Odds Ratio (log scale)", xlims = or_xlim,
  panel_title = "Other Gram-negative Subgroup"
)

p_weak_no_pseudo <- make_panel(
  weak_no_pseudo, weak_prior_no_pseudo,
  fill_col = col_no_pseudo, threshold = 0, direction = "below",
  prob_pct = weak_prob_no_pseudo,
  x_breaks = log(or_breaks), x_labels = or_breaks,
  x_title = "Odds Ratio (log scale)", xlims = or_xlim,
  panel_title = "Other Gram-negative Subgroup"
)

# Row 2: Pseudomonas — P(OR > 1), fill above 0, orange
p_info_pseudo <- make_panel(
  info_pseudo, info_prior_pseudo,
  fill_col = col_pseudo, threshold = 0, direction = "above",
  prob_pct = info_prob_pseudo,
  x_breaks = log(or_breaks), x_labels = or_breaks,
  x_title = "Odds Ratio (log scale)", xlims = or_xlim,
  panel_title = "Pseudomonas Subgroup"
)

p_weak_pseudo <- make_panel(
  weak_pseudo, weak_prior_pseudo,
  fill_col = col_pseudo, threshold = 0, direction = "above",
  prob_pct = weak_prob_pseudo,
  x_breaks = log(or_breaks), x_labels = or_breaks,
  x_title = "Odds Ratio (log scale)", xlims = or_xlim,
  panel_title = "Pseudomonas Subgroup"
)

# Row 3: Interaction — P(ROR > 1), fill above 0, green
p_info_interaction <- make_panel(
  info_interaction, info_prior_interaction,
  fill_col = col_interaction, threshold = 0, direction = "above",
  prob_pct = info_prob_interaction,
  x_breaks = log(ror_breaks), x_labels = ror_breaks,
  x_title = "Ratio of Odds Ratios (log scale)", xlims = ror_xlim,
  panel_title = "Subgroup Difference (Interaction)"
)

p_weak_interaction <- make_panel(
  weak_interaction, weak_prior_interaction,
  fill_col = col_interaction, threshold = 0, direction = "above",
  prob_pct = weak_prob_interaction,
  x_breaks = log(ror_breaks), x_labels = ror_breaks,
  x_title = "Ratio of Odds Ratios (log scale)", xlims = ror_xlim,
  panel_title = "Subgroup Difference (Interaction)"
)

# Column 3: Skeptical prior panels
p_skeptic_no_pseudo <- make_panel(
  skeptic_no_pseudo, skeptic_prior_no_pseudo,
  fill_col = col_no_pseudo, threshold = 0, direction = "below",
  prob_pct = skeptic_prob_no_pseudo,
  x_breaks = log(or_breaks), x_labels = or_breaks,
  x_title = "Odds Ratio (log scale)", xlims = or_xlim,
  panel_title = "Other Gram-negative Subgroup"
)

p_skeptic_pseudo <- make_panel(
  skeptic_pseudo, skeptic_prior_pseudo,
  fill_col = col_pseudo, threshold = 0, direction = "above",
  prob_pct = skeptic_prob_pseudo,
  x_breaks = log(or_breaks), x_labels = or_breaks,
  x_title = "Odds Ratio (log scale)", xlims = or_xlim,
  panel_title = "Pseudomonas Subgroup"
)

p_skeptic_interaction <- make_panel(
  skeptic_interaction, skeptic_prior_interaction,
  fill_col = col_interaction, threshold = 0, direction = "above",
  prob_pct = skeptic_prob_interaction,
  x_breaks = log(ror_breaks), x_labels = ror_breaks,
  x_title = "Ratio of Odds Ratios (log scale)", xlims = ror_xlim,
  panel_title = "Subgroup Difference (Interaction)"
)

# ---- Step 9: Assemble with patchwork ----

# Column header labels
library(grid)

header_primary <- wrap_elements(
  textGrob("Primary Analysis\n(Historical Prior)",
           gp = gpar(fontsize = 12, fontface = "bold"))
)
header_sensitivity <- wrap_elements(
  textGrob("Sensitivity Analysis\n(Weakly Informative Prior)",
           gp = gpar(fontsize = 12, fontface = "bold"))
)
header_skeptic <- wrap_elements(
  textGrob("Sensitivity Analysis\n(Skeptical Prior)",
           gp = gpar(fontsize = 12, fontface = "bold"))
)

# Layout: header row + 3 data rows, 3 columns
figure_letter <- (header_primary | header_skeptic | header_sensitivity) /
  (p_info_pseudo    | p_skeptic_pseudo    | p_weak_pseudo) /
  (p_info_no_pseudo | p_skeptic_no_pseudo | p_weak_no_pseudo) /
  (p_info_interaction | p_skeptic_interaction | p_weak_interaction) +
  plot_layout(heights = c(0.8, 3, 3, 3))

tiff("figures/figure_letter.tiff", width = 12, height = 10, units = "in", res = 300, bg = "white")
print(figure_letter)
grid.lines(x = unit(c(1/3, 1/3), "npc"), y = unit(c(0, 1), "npc"),
           gp = gpar(col = "black", lwd = 0.5))
grid.lines(x = unit(c(2/3, 2/3), "npc"), y = unit(c(0, 1), "npc"),
           gp = gpar(col = "black", lwd = 0.5))
dev.off()

png("figures/figure_letter.png", width = 12, height = 10, units = "in", res = 300, bg = "white")
print(figure_letter)
grid.lines(x = unit(c(1/3, 1/3), "npc"), y = unit(c(0, 1), "npc"),
           gp = gpar(col = "black", lwd = 0.5))
grid.lines(x = unit(c(2/3, 2/3), "npc"), y = unit(c(0, 1), "npc"),
           gp = gpar(col = "black", lwd = 0.5))
dev.off()

message("Saved: figures/figure_letter.tiff and figures/figure_letter.png")

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

median_hdi(info_interaction) |> mutate(across(y:ymax, ~exp(.)))

# =========================================================================
# Step 10 — RR-scale pipeline
#
# Uses the same three fits (info_fit, weak_fit, skeptic_fit). Posterior log-RR
# draws come from marginaleffects::avg_comparisons() with
# comparison = "lnratioavg" and by = "pseudo": this returns
# log( mean(p | treat = 1) / mean(p | treat = 0) ) within each pseudo level.
# The log RoRR is formed as the difference of the two by-group log-RR draws.
# Prior log-RR draws come from joint prior coefficient draws.
# =========================================================================

# ---- Step 10.1: Helpers ----

get_post_logrr_draws <- function(fit) {
  cmp <- avg_comparisons(
    fit,
    variables  = "treat",
    by         = "pseudo",
    comparison = "lnratioavg"
  )
  d <- posterior_draws(cmp, shape = "long") |>
    select(drawid, pseudo, draw) |>
    pivot_wider(
      names_from  = pseudo,
      values_from = draw,
      names_glue  = "logrr_p{pseudo}"
    )
  log_rr_no_pseudo <- d$logrr_p0
  log_rr_pseudo    <- d$logrr_p1
  log_rorr         <- log_rr_pseudo - log_rr_no_pseudo
  list(log_rr_no_pseudo = log_rr_no_pseudo,
       log_rr_pseudo    = log_rr_pseudo,
       log_rorr         = log_rorr)
}

# Fixed-baseline prior on log-RR: convert the prior on
# log-OR (and the log-OR interaction) to the log-RR scale using the corresponding
# the BALANCE 14-day arm event rate in each subgroup. 
baseline_no_pseudo <- balance$fourteen_events[balance$subgroup == "No Pseudomonas"] /
                      balance$fourteen_total [balance$subgroup == "No Pseudomonas"]
baseline_pseudo    <- balance$fourteen_events[balance$subgroup == "Only Pseudomonas"] /
                      balance$fourteen_total [balance$subgroup == "Only Pseudomonas"]

get_prior_logrr_draws <- function(post_df,
                                  p0_no_pseudo = baseline_no_pseudo,
                                  p0_pseudo    = baseline_pseudo) {
  logit_p0_no_pseudo <- qlogis(p0_no_pseudo)
  logit_p0_pseudo    <- qlogis(p0_pseudo)
  p_t1_p0 <- plogis(logit_p0_no_pseudo + post_df$prior_b_treat)
  p_t1_p1 <- plogis(logit_p0_pseudo    + post_df$prior_b_treat +
                    post_df$`prior_b_treat:pseudo`)
  log_rr_no_pseudo <- log(p_t1_p0 / p0_no_pseudo)
  log_rr_pseudo    <- log(p_t1_p1 / p0_pseudo)
  log_rorr         <- log_rr_pseudo - log_rr_no_pseudo
  list(log_rr_no_pseudo = log_rr_no_pseudo,
       log_rr_pseudo    = log_rr_pseudo,
       log_rorr         = log_rorr)
}

# ---- Step 10.2: Posterior + prior log-RR draws for each fit ----

info_rr_post_draws    <- get_post_logrr_draws(info_fit)
weak_rr_post_draws    <- get_post_logrr_draws(weak_fit)
skeptic_rr_post_draws <- get_post_logrr_draws(skeptic_fit)

info_rr_prior_draws    <- get_prior_logrr_draws(info_post)
weak_rr_prior_draws    <- get_prior_logrr_draws(weak_post)
skeptic_rr_prior_draws <- get_prior_logrr_draws(skeptic_post)

info_rr_no_pseudo   <- info_rr_post_draws$log_rr_no_pseudo
info_rr_pseudo      <- info_rr_post_draws$log_rr_pseudo
info_rr_interaction <- info_rr_post_draws$log_rorr
info_prior_rr_no_pseudo   <- info_rr_prior_draws$log_rr_no_pseudo
info_prior_rr_pseudo      <- info_rr_prior_draws$log_rr_pseudo
info_prior_rr_interaction <- info_rr_prior_draws$log_rorr

weak_rr_no_pseudo   <- weak_rr_post_draws$log_rr_no_pseudo
weak_rr_pseudo      <- weak_rr_post_draws$log_rr_pseudo
weak_rr_interaction <- weak_rr_post_draws$log_rorr
weak_prior_rr_no_pseudo   <- weak_rr_prior_draws$log_rr_no_pseudo
weak_prior_rr_pseudo      <- weak_rr_prior_draws$log_rr_pseudo
weak_prior_rr_interaction <- weak_rr_prior_draws$log_rorr

skeptic_rr_no_pseudo   <- skeptic_rr_post_draws$log_rr_no_pseudo
skeptic_rr_pseudo      <- skeptic_rr_post_draws$log_rr_pseudo
skeptic_rr_interaction <- skeptic_rr_post_draws$log_rorr
skeptic_prior_rr_no_pseudo   <- skeptic_rr_prior_draws$log_rr_no_pseudo
skeptic_prior_rr_pseudo      <- skeptic_rr_prior_draws$log_rr_pseudo
skeptic_prior_rr_interaction <- skeptic_rr_prior_draws$log_rorr

# ---- Step 10.3: Posterior probabilities on RR scale ----

info_prob_rr_no_pseudo    <- round(mean(info_rr_no_pseudo    < 0) * 100, 0)
info_prob_rr_pseudo       <- round(mean(info_rr_pseudo       > 0) * 100, 0)
info_prob_rr_interaction  <- round(mean(info_rr_interaction  > 0) * 100, 0)
weak_prob_rr_no_pseudo    <- round(mean(weak_rr_no_pseudo    < 0) * 100, 0)
weak_prob_rr_pseudo       <- round(mean(weak_rr_pseudo       > 0) * 100, 0)
weak_prob_rr_interaction  <- round(mean(weak_rr_interaction  > 0) * 100, 0)
skeptic_prob_rr_no_pseudo   <- round(mean(skeptic_rr_no_pseudo   < 0) * 100, 0)
skeptic_prob_rr_pseudo      <- round(mean(skeptic_rr_pseudo      > 0) * 100, 0)
skeptic_prob_rr_interaction <- round(mean(skeptic_rr_interaction > 0) * 100, 0)

# ---- Step 10.4: Build nine RR panels with the existing make_panel() ----

rr_breaks   <- c(0.1, 0.2, 0.3, 0.5, 0.7, 1, 1.5, 2, 3, 5, 10)
rorr_breaks <- c(0.1, 0.2, 0.3, 0.5, 0.7, 1, 1.5, 2, 3, 5, 10)
rr_xlim     <- log(c(0.2, 5))
rorr_xlim   <- log(c(0.2, 5))

p_info_rr_no_pseudo <- make_panel(
  info_rr_no_pseudo, info_prior_rr_no_pseudo,
  fill_col = col_no_pseudo, threshold = 0, direction = "below",
  prob_pct = info_prob_rr_no_pseudo,
  x_breaks = log(rr_breaks), x_labels = rr_breaks,
  x_title = "Risk Ratio (log scale)", xlims = rr_xlim,
  panel_title = "Other Gram-negative Subgroup"
)
p_weak_rr_no_pseudo <- make_panel(
  weak_rr_no_pseudo, weak_prior_rr_no_pseudo,
  fill_col = col_no_pseudo, threshold = 0, direction = "below",
  prob_pct = weak_prob_rr_no_pseudo,
  x_breaks = log(rr_breaks), x_labels = rr_breaks,
  x_title = "Risk Ratio (log scale)", xlims = rr_xlim,
  panel_title = "Other Gram-negative Subgroup"
)
p_skeptic_rr_no_pseudo <- make_panel(
  skeptic_rr_no_pseudo, skeptic_prior_rr_no_pseudo,
  fill_col = col_no_pseudo, threshold = 0, direction = "below",
  prob_pct = skeptic_prob_rr_no_pseudo,
  x_breaks = log(rr_breaks), x_labels = rr_breaks,
  x_title = "Risk Ratio (log scale)", xlims = rr_xlim,
  panel_title = "Other Gram-negative Subgroup"
)

p_info_rr_pseudo <- make_panel(
  info_rr_pseudo, info_prior_rr_pseudo,
  fill_col = col_pseudo, threshold = 0, direction = "above",
  prob_pct = info_prob_rr_pseudo,
  x_breaks = log(rr_breaks), x_labels = rr_breaks,
  x_title = "Risk Ratio (log scale)", xlims = rr_xlim,
  panel_title = "Pseudomonas Subgroup"
)
p_weak_rr_pseudo <- make_panel(
  weak_rr_pseudo, weak_prior_rr_pseudo,
  fill_col = col_pseudo, threshold = 0, direction = "above",
  prob_pct = weak_prob_rr_pseudo,
  x_breaks = log(rr_breaks), x_labels = rr_breaks,
  x_title = "Risk Ratio (log scale)", xlims = rr_xlim,
  panel_title = "Pseudomonas Subgroup"
)
p_skeptic_rr_pseudo <- make_panel(
  skeptic_rr_pseudo, skeptic_prior_rr_pseudo,
  fill_col = col_pseudo, threshold = 0, direction = "above",
  prob_pct = skeptic_prob_rr_pseudo,
  x_breaks = log(rr_breaks), x_labels = rr_breaks,
  x_title = "Risk Ratio (log scale)", xlims = rr_xlim,
  panel_title = "Pseudomonas Subgroup"
)

p_info_rr_interaction <- make_panel(
  info_rr_interaction, info_prior_rr_interaction,
  fill_col = col_interaction, threshold = 0, direction = "above",
  prob_pct = info_prob_rr_interaction,
  x_breaks = log(rorr_breaks), x_labels = rorr_breaks,
  x_title = "Ratio of Risk Ratios (log scale)", xlims = rorr_xlim,
  panel_title = "Subgroup Difference (Interaction)"
)
p_weak_rr_interaction <- make_panel(
  weak_rr_interaction, weak_prior_rr_interaction,
  fill_col = col_interaction, threshold = 0, direction = "above",
  prob_pct = weak_prob_rr_interaction,
  x_breaks = log(rorr_breaks), x_labels = rorr_breaks,
  x_title = "Ratio of Risk Ratios (log scale)", xlims = rorr_xlim,
  panel_title = "Subgroup Difference (Interaction)"
)
p_skeptic_rr_interaction <- make_panel(
  skeptic_rr_interaction, skeptic_prior_rr_interaction,
  fill_col = col_interaction, threshold = 0, direction = "above",
  prob_pct = skeptic_prob_rr_interaction,
  x_breaks = log(rorr_breaks), x_labels = rorr_breaks,
  x_title = "Ratio of Risk Ratios (log scale)", xlims = rorr_xlim,
  panel_title = "Subgroup Difference (Interaction)"
)

# ---- Step 10.5: Assemble figure_letter_rr.tiff ----

figure_letter_rr <- (header_primary | header_skeptic | header_sensitivity) /
  (p_info_rr_pseudo      | p_skeptic_rr_pseudo      | p_weak_rr_pseudo) /
  (p_info_rr_no_pseudo   | p_skeptic_rr_no_pseudo   | p_weak_rr_no_pseudo) /
  (p_info_rr_interaction | p_skeptic_rr_interaction | p_weak_rr_interaction) +
  plot_layout(heights = c(0.8, 3, 3, 3))

tiff("figures/figure_letter_rr.tiff", width = 12, height = 10, units = "in", res = 300, bg = "white")
print(figure_letter_rr)
grid.lines(x = unit(c(1/3, 1/3), "npc"), y = unit(c(0, 1), "npc"),
           gp = gpar(col = "black", lwd = 0.5))
grid.lines(x = unit(c(2/3, 2/3), "npc"), y = unit(c(0, 1), "npc"),
           gp = gpar(col = "black", lwd = 0.5))
dev.off()

png("figures/figure_letter_rr.png", width = 12, height = 10, units = "in", res = 300, bg = "white")
print(figure_letter_rr)
grid.lines(x = unit(c(1/3, 1/3), "npc"), y = unit(c(0, 1), "npc"),
           gp = gpar(col = "black", lwd = 0.5))
grid.lines(x = unit(c(2/3, 2/3), "npc"), y = unit(c(0, 1), "npc"),
           gp = gpar(col = "black", lwd = 0.5))
dev.off()

message("Saved: figures/figure_letter_rr.tiff and figures/figure_letter_rr.png")

# ---- Step 10.6: Print RR-scale values for all three models ----

print_rr_block <- function(label, no_pseudo, pseudo, interaction,
                           prob_no_pseudo, prob_pseudo, prob_interaction) {
  cat("\n=== RR-scale values: ", label, " ===\n", sep = "")
  cat("No Pseudo RR:", round(exp(median(no_pseudo)), 2), "\n")
  cat("No Pseudo 95% CrI:", round(exp(quantile(no_pseudo, 0.025)), 2), "-",
      round(exp(quantile(no_pseudo, 0.975)), 2), "\n")
  cat("No Pseudo P(RR<1):", prob_no_pseudo, "%\n")
  cat("Pseudo RR:", round(exp(median(pseudo)), 2), "\n")
  cat("Pseudo 95% CrI:", round(exp(quantile(pseudo, 0.025)), 2), "-",
      round(exp(quantile(pseudo, 0.975)), 2), "\n")
  cat("Pseudo P(RR>1):", prob_pseudo, "%\n")
  cat("RoRR median:", round(exp(median(interaction)), 2), "\n")
  cat("RoRR 95% CrI:", round(exp(quantile(interaction, 0.025)), 2), "-",
      round(exp(quantile(interaction, 0.975)), 2), "\n")
  cat("Interaction P(RoRR>1):", prob_interaction, "%\n")
}

print_rr_block("Primary (historical prior)",
               info_rr_no_pseudo, info_rr_pseudo, info_rr_interaction,
               info_prob_rr_no_pseudo, info_prob_rr_pseudo, info_prob_rr_interaction)

print_rr_block("Sensitivity (weakly informative prior)",
               weak_rr_no_pseudo, weak_rr_pseudo, weak_rr_interaction,
               weak_prob_rr_no_pseudo, weak_prob_rr_pseudo, weak_prob_rr_interaction)

print_rr_block("Sensitivity (skeptical prior)",
               skeptic_rr_no_pseudo, skeptic_rr_pseudo, skeptic_rr_interaction,
               skeptic_prob_rr_no_pseudo, skeptic_prob_rr_pseudo, skeptic_prob_rr_interaction)

# =========================================================================
# Step 11 — Analytic prior summaries on the RR / RoRR scale (all 3 models)
#
# The model is fit on the log-OR scale; the paper reports on the risk-ratio
# scale. This block translates each of the three prior pairs (historical,
# skeptical, weakly informative) into median + 95% on the RR (treatment) and
# RoRR (interaction) scale.
#
# Conversion: Zhang & Yu (JAMA 1998; doi:10.1001/jama.280.19.1690) closed-form
#   RR = OR / (p0 * (OR - 1) + 1)
# advocated for fixed-baseline use in Bayesian / interval translation by
# Doi SA et al. (J Clin Epidemiol 2020; doi:10.1016/j.jclinepi.2020.08.019).
# Algebraically identical to RR = plogis(qlogis(p0) + log_OR) / p0 (used here
# and in get_prior_logrr_draws() above for the prior-draw pipeline).
#
# Baselines (BALANCE 14-day event rates):
#   baseline_no_pseudo for the treatment effect
#   baseline_pseudo    for the interaction (RoRR)
#
# For the interaction we plug in beta_treat = paired treat-prior mean
# (historical: treat_prior_mean; skeptical / weakly: 0), so the RoRR reflects
# the joint prior median actually assumed by each model.
# =========================================================================

# Helper: Normal(mu, sigma) on log-OR  ->  median + 95% on RR at fixed p0
or_prior_to_rr <- function(mu, sigma, p0) {
  q975 <- qnorm(0.975)
  list(
    med = plogis(qlogis(p0) + mu) / p0,
    lo  = plogis(qlogis(p0) + mu - q975 * sigma) / p0,
    hi  = plogis(qlogis(p0) + mu + q975 * sigma) / p0
  )
}

# Helper: Normal(mu_int, sigma_int) on log-ROR  ->  median + 95% on RoRR,
# evaluated with beta_treat held at mu_treat. With mu_treat = 0 this reduces
# to the simple form (RR_no_pseudo collapses to 1 and RoRR = RR_pseudo).
interaction_prior_to_rorr <- function(mu_int, sigma_int, mu_treat,
                                      p0_no_pseudo, p0_pseudo) {
  q975 <- qnorm(0.975)
  rr_no_pseudo <- plogis(qlogis(p0_no_pseudo) + mu_treat) / p0_no_pseudo
  rr_p <- function(b_int) plogis(qlogis(p0_pseudo) + mu_treat + b_int) / p0_pseudo
  list(
    med = rr_p(mu_int) / rr_no_pseudo,
    lo  = rr_p(mu_int - q975 * sigma_int) / rr_no_pseudo,
    hi  = rr_p(mu_int + q975 * sigma_int) / rr_no_pseudo
  )
}

# Historical: treat_prior_mean / sd from the meta-analysis (Step 1);
#             interaction_mean / sd from Yahav 2019 (Step 1).
hist_treat_rr        <- or_prior_to_rr(treat_prior_mean, treat_prior_sd, baseline_no_pseudo)
hist_interaction_rorr <- interaction_prior_to_rorr(
  interaction_mean, interaction_sd, treat_prior_mean,
  baseline_no_pseudo, baseline_pseudo
)

# Skeptical: Normal(0, 0.71) on log-OR for treat; Normal(0, 0.36) for interaction.
skeptic_treat_rr         <- or_prior_to_rr(0, 0.71, baseline_no_pseudo)
skeptic_interaction_rorr <- interaction_prior_to_rorr(
  0, 0.36, 0, baseline_no_pseudo, baseline_pseudo
)

# Weakly informative: Normal(0, 0.82) on log-OR for treat; Normal(0, 0.71) for interaction.
weak_treat_rr         <- or_prior_to_rr(0, 0.82, baseline_no_pseudo)
weak_interaction_rorr <- interaction_prior_to_rorr(
  0, 0.71, 0, baseline_no_pseudo, baseline_pseudo
)

# Wide-format summary table (3 rows x 3 columns) used by report.qmd.
fmt_med_ci <- function(x) sprintf("%.2f (%.2f - %.2f)", x$med, x$lo, x$hi)

prior_rr_summary <- tibble::tibble(
  Model = c("Historical (primary)", "Skeptical", "Weakly informative"),
  `Treatment effect, RR (median, 95% interval)` = c(
    fmt_med_ci(hist_treat_rr),
    fmt_med_ci(skeptic_treat_rr),
    fmt_med_ci(weak_treat_rr)
  ),
  `Interaction, RoRR (median, 95% interval)` = c(
    fmt_med_ci(hist_interaction_rorr),
    fmt_med_ci(skeptic_interaction_rorr),
    fmt_med_ci(weak_interaction_rorr)
  )
)

cat("\n=== Prior summaries on the RR / RoRR scale ===\n")
cat(sprintf("baseline_no_pseudo: %.4f\n", baseline_no_pseudo))
cat(sprintf("baseline_pseudo:    %.4f\n", baseline_pseudo))
print(prior_rr_summary)
