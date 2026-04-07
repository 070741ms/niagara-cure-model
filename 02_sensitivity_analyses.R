# =============================================================================
# 02_sensitivity_analyses.R
#
# Revision sensitivity analyses (Supplementary Figures S1-S5 and Tables S1-S5):
#   S1: Reconstruction uncertainty (perturbation)
#   S2: Follow-up truncation
#   S3: Empirical plateau assessment (Maller-Zhou approach)
#   S4: Akaike weights for cure regression models
#   S5: Conditional event-free survival
#
# Input:  data/efs_ripd_excel.csv
# Output: outputs/revision_outputs/
# =============================================================================

# =============================================================================
# NIAGARA Trial EFS — Revision Sensitivity Analyses  (v2)
# =============================================================================
#
# Addresses:
#   R1 Major #1  : Reconstruction uncertainty (random-noise perturbation)
#   R1 Major #2  : Follow-up truncation sensitivity + empirical plateau assessment
#   R1 Major #3  : Akaike weights for cure-regression models
#   R1 Suggested : Conditional survival analysis  CS(2 yr | s yr)
#
# Input:
#   efs_ripd_excel.csv
#
# Output (folders created automatically):
#   revision_outputs/
#     01_reconstruction_uncertainty/
#     02_followup_truncation/
#     03_maller_zhou/
#     04_akaike_weights/
#     05_conditional_survival/
#
# Dependencies: survival, ggplot2, gridExtra
#   (grid is used via grid:: without explicit library() call)
# =============================================================================


# =============================================================================
# 0.  Setup
# =============================================================================

required_packages <- c("survival", "ggplot2", "gridExtra")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  library(pkg, character.only = TRUE)
}

# ---- Output directories ------------------------------------------------------
dir.create("outputs", showWarnings = FALSE)
root_dir <- "outputs/revision_outputs"
subdirs  <- file.path(root_dir, sprintf("%02d_%s", 1:5,
                                        c("reconstruction_uncertainty", "followup_truncation",
                                          "maller_zhou", "akaike_weights", "conditional_survival")))
for (d in c(root_dir, subdirs)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# ---- Colour palette (matches original manuscript / NIAGARA paper) -----------
COL_CTL <- "#4A90D9"
COL_DUR <- "#F5A623"

# ---- Publication theme -------------------------------------------------------
theme_pub <- function(base_size = 13) {
  theme_classic(base_size = base_size) +
    theme(
      plot.background   = element_rect(fill = "white", colour = NA),
      panel.background  = element_rect(fill = "white", colour = NA),
      plot.title        = element_text(size = base_size + 1, face = "bold",   hjust = 0.5),
      plot.subtitle     = element_text(size = base_size - 1, colour = "gray35", hjust = 0.5),
      axis.title        = element_text(size = base_size,     face = "bold"),
      axis.text         = element_text(size = base_size - 1, colour = "black"),
      axis.line         = element_line(colour = "black", linewidth = 0.5),
      axis.ticks        = element_line(colour = "black", linewidth = 0.5),
      legend.position   = "top",
      legend.title      = element_blank(),
      legend.text       = element_text(size = base_size - 1),
      legend.background = element_rect(fill = "white", colour = NA),
      legend.key        = element_rect(fill = "white", colour = NA),
      plot.margin       = margin(15, 15, 15, 15)
    )
}

# Helper: signed percentage-point label (avoids "+-3 pp" for negative values)
pp_label <- function(x) ifelse(x >= 0, paste0("+", x, " pp"), paste0(x, " pp"))

cat("\n=================================================================\n")
cat("   NIAGARA Revision Sensitivity Analyses  (v2)\n")
cat("=================================================================\n\n")


# =============================================================================
# 1.  Data Loading & Validation
# =============================================================================

input_file <- "efs_ripd_excel.csv"
if (!file.exists(input_file)) {
  stop("'efs_ripd_excel.csv' not found.\n",
       "Ensure the working directory is the repository root.")
}

df        <- read.csv(input_file, stringsAsFactors = FALSE)
df$time   <- as.numeric(df$time)
df$status <- as.integer(df$status)

# ---- Input validation --------------------------------------------------------
if (!all(df$status %in% c(0L, 1L))) {
  stop("Column 'status' must contain only 0 (censored) and 1 (event).\n",
       "Found values: ", paste(sort(unique(df$status)), collapse = ", "))
}

valid_arms   <- c("Control", "Durvalumab")
unknown_arms <- setdiff(unique(df$arm), valid_arms)
if (length(unknown_arms) > 0L) {
  stop("Unknown arm values: ", paste(unknown_arms, collapse = ", "),
       ".\nExpected 'Control' and 'Durvalumab'.")
}

df$arm  <- factor(df$arm, levels = valid_arms)
df$time <- pmax(df$time, 0.001)   # log-normal requires t > 0

dfC <- df[df$arm == "Control",    ]
dfT <- df[df$arm == "Durvalumab", ]

cat(sprintf("Data loaded:  Total N = %d\n", nrow(df)))
cat(sprintf("  Control    : n = %d,  events = %d (%.1f%%)\n",
            nrow(dfC), sum(dfC$status), 100 * mean(dfC$status)))
cat(sprintf("  Durvalumab : n = %d,  events = %d (%.1f%%)\n\n",
            nrow(dfT), sum(dfT$status), 100 * mean(dfT$status)))


# =============================================================================
# 2.  Core Mixture Cure Model Functions  (log-normal)
# All sensitivity analyses use log-normal throughout.
# Rationale: log-normal was identified as best-fitting by AIC in 01_main_analysis.R
# (see outputs/tables/Table_ModelSelection.csv). Fixing the distribution isolates
# the sensitivity being tested from distributional uncertainty.
# =============================================================================
# All sensitivity analyses below use the log-normal mixture cure model.
# Rationale: log-normal was identified as the best-fitting distribution by AIC
# in the primary analysis (01_main_analysis.R, Table_ModelSelection.csv).
# Fixing the distribution across all sensitivity analyses avoids confounding
# the sensitivity results with distributional uncertainty, and is consistent
# with the approach described in the manuscript.
# =============================================================================

inv_logit <- function(z) 1 / (1 + exp(-z))

nll_lognormal_cure <- function(par, t, status) {
  pi    <- pmin(pmax(inv_logit(par[1]), 1e-9), 1 - 1e-9)
  mu    <- par[2]
  sigma <- exp(par[3])
  t     <- pmax(t, 1e-12)
  Su    <- 1 - plnorm(t, meanlog = mu, sdlog = sigma)
  fu    <- dlnorm(t,    meanlog = mu, sdlog = sigma)
  lik   <- ifelse(status == 1L, (1 - pi) * fu, pi + (1 - pi) * Su)
  -sum(log(pmax(lik, 1e-300)))
}

fit_cure <- function(t, status) {
  # Guard: need at least a few events for a meaningful cure model
  n_ev <- sum(status == 1L, na.rm = TRUE)
  if (n_ev == 0L) {
    warning("No events; cure model cannot be fitted.")
    return(NULL)
  }
  if (n_ev < 5L) {
    warning(sprintf("Very few events (%d); cure model estimates may be unreliable.", n_ev))
  }
  
  med_t <- max(median(t[status == 1L], na.rm = TRUE), 0.1)
  init  <- c(0, log(med_t), 0)
  
  opt <- tryCatch(
    optim(par = init, fn = nll_lognormal_cure,
          t = t, status = status,
          method = "L-BFGS-B",
          lower  = c(-10, -5, log(0.05)),
          upper  = c( 10, 10, log(5.00))),
    error = function(e) NULL
  )
  if (is.null(opt)) return(NULL)
  
  pi     <- inv_logit(opt$par[1])
  loglik <- -opt$value
  list(pi = pi, mu = opt$par[2], sigma = exp(opt$par[3]),
       a = opt$par[1], loglik = loglik,
       aic = -2 * loglik + 2 * 3)
}

# ---- Reference fit -----------------------------------------------------------
cat("Fitting reference cure model (log-normal, full data)...\n")
ref_C <- fit_cure(dfC$time, dfC$status)
ref_T <- fit_cure(dfT$time, dfT$status)
if (is.null(ref_C) || is.null(ref_T)) {
  stop("Reference cure model failed to converge. Check the input data.")
}
ref_delta <- ref_T$pi - ref_C$pi
cat(sprintf("  pi_Control = %.3f,  pi_Durvalumab = %.3f,  Delta-pi = +%.3f\n\n",
            ref_C$pi, ref_T$pi, ref_delta))


# =============================================================================
# 3.  ANALYSIS 1 — Reconstruction Uncertainty (Perturbation)
#     R1 Major Issue #1
# =============================================================================

cat("--- Analysis 1: Reconstruction Uncertainty ---\n")

perturb_analysis <- function(df, sd_noise, n_iter = 200, seed = 42) {
  set.seed(seed)
  out <- vector("list", n_iter)
  for (i in seq_len(n_iter)) {
    df_p      <- df
    df_p$time <- pmax(df$time + rnorm(nrow(df), 0, sd_noise), 0.001)
    dC_p <- df_p[df_p$arm == "Control",    ]
    dT_p <- df_p[df_p$arm == "Durvalumab", ]
    fC <- tryCatch(fit_cure(dC_p$time, dC_p$status), error = function(e) NULL)
    fT <- tryCatch(fit_cure(dT_p$time, dT_p$status), error = function(e) NULL)
    if (!is.null(fC) && !is.null(fT)) {
      out[[i]] <- data.frame(iter = i, sd = sd_noise,
                             piC = fC$pi, piT = fT$pi,
                             delta = fT$pi - fC$pi)
    }
  }
  do.call(rbind, Filter(Negate(is.null), out))
}

cat(sprintf("  Running 200 iterations at noise SD = 0.5 months...\n"))
pert05 <- perturb_analysis(df, sd_noise = 0.5, n_iter = 200, seed = 42)
cat(sprintf("  Running 200 iterations at noise SD = 1.0 months...\n"))
pert10 <- perturb_analysis(df, sd_noise = 1.0, n_iter = 200, seed = 43)

pert_all          <- rbind(pert05, pert10)
pert_all$sd_label <- factor(
  ifelse(pert_all$sd == 0.5, "Noise SD = 0.5 months", "Noise SD = 1.0 months"),
  levels = c("Noise SD = 0.5 months", "Noise SD = 1.0 months")
)

# Summary table
smry_pert <- function(x) {
  data.frame(Mean = round(mean(x$delta)*100, 2), SD = round(sd(x$delta)*100, 2),
             Min  = round(min(x$delta)*100, 2),  Max = round(max(x$delta)*100, 2),
             CI95_lo = round(quantile(x$delta, 0.025)*100, 2),
             CI95_hi = round(quantile(x$delta, 0.975)*100, 2),
             N_converged = nrow(x))
}
pert_table <- rbind(
  data.frame(Noise_SD_months = 0, Mean = round(ref_delta*100,2),
             SD=NA, Min=NA, Max=NA, CI95_lo=NA, CI95_hi=NA, N_converged=1),
  cbind(Noise_SD_months = 0.5, smry_pert(pert05)),
  cbind(Noise_SD_months = 1.0, smry_pert(pert10))
)
write.csv(pert_table,
          file.path(subdirs[1], "Table_perturbation_summary.csv"),
          row.names = FALSE)

# ---- Figure: violin + box ---------------------------------------------------
# Reference value is placed in the subtitle (not annotated inside the plot)
# because the dashed line passes through the widest part of both violins
ref_y <- ref_delta * 100

fig1 <- ggplot(pert_all, aes(x = sd_label, y = delta * 100, fill = sd_label)) +
  geom_violin(alpha = 0.50, colour = NA, trim = FALSE) +
  geom_boxplot(width = 0.13, colour = "gray20", fill = "white",
               outlier.size = 0.7, linewidth = 0.55) +
  geom_hline(yintercept = ref_y,
             linetype = "dashed", colour = "black", linewidth = 0.85) +
  scale_fill_manual(values = c("Noise SD = 0.5 months" = "#7FBBEA",
                               "Noise SD = 1.0 months" = COL_DUR),
                    guide = "none") +
  # Fix #5: signed label — no "+-" artefact
  scale_y_continuous(labels = function(x) pp_label(round(x, 0))) +
  labs(title    = "Reconstruction Uncertainty: Perturbation Sensitivity",
       subtitle  = paste0(
         "Distribution of \u0394\u03c0 across 200 datasets per noise level\n",
         sprintf("Dashed line = reference \u0394\u03c0 = +%.1f pp (unperturbed data)",
                 ref_y)),
       x = "Gaussian noise added to event / censoring times",
       y = "Cure fraction difference, \u0394\u03c0 (percentage points)") +
  theme_pub()

ggsave(file.path(subdirs[1], "Fig_perturbation_sensitivity.png"),
       fig1, width = 8, height = 6, dpi = 300, bg = "white")

cat(sprintf("  SD=0.5 mo: mean \u0394\u03c0 = +%.1f pp,  SD = %.1f pp,  95%% range [%s, %s]\n",
            mean(pert05$delta)*100, sd(pert05$delta)*100,
            pp_label(round(quantile(pert05$delta,.025)*100,1)),
            pp_label(round(quantile(pert05$delta,.975)*100,1))))
cat(sprintf("  SD=1.0 mo: mean \u0394\u03c0 = +%.1f pp,  SD = %.1f pp,  95%% range [%s, %s]\n",
            mean(pert10$delta)*100, sd(pert10$delta)*100,
            pp_label(round(quantile(pert10$delta,.025)*100,1)),
            pp_label(round(quantile(pert10$delta,.975)*100,1))))
cat("  >>> Saved: 01_reconstruction_uncertainty/\n\n")


# =============================================================================
# 4.  ANALYSIS 2 — Follow-up Truncation Sensitivity
#     R1 Major Issue #2 (part a)
# =============================================================================

cat("--- Analysis 2: Follow-up Truncation Sensitivity ---\n")

truncate_fit <- function(df, cutoff, B_boot = 300, seed = 100) {
  # Fix #12: explicit branch for full data (avoid relying on time > Inf == FALSE)
  if (is.finite(cutoff)) {
    df_t           <- df
    beyond         <- df_t$time > cutoff
    df_t$status[beyond] <- 0L
    df_t$time[beyond]   <- cutoff
    df_t$time           <- pmax(df_t$time, 0.001)
  } else {
    df_t <- df
  }
  
  dC <- df_t[df_t$arm == "Control",    ]
  dT <- df_t[df_t$arm == "Durvalumab", ]
  
  fC <- fit_cure(dC$time, dC$status)
  fT <- fit_cure(dT$time, dT$status)
  if (is.null(fC) || is.null(fT)) return(NULL)
  
  # Bootstrap — Fix #8: separate ok per arm
  set.seed(seed)
  piC_b <- piT_b <- dlt_b <- numeric(B_boot)
  for (b in seq_len(B_boot)) {
    idxC <- sample(nrow(dC), nrow(dC), replace = TRUE)
    idxT <- sample(nrow(dT), nrow(dT), replace = TRUE)
    fbC  <- tryCatch(fit_cure(dC$time[idxC], dC$status[idxC]), error = function(e) NULL)
    fbT  <- tryCatch(fit_cure(dT$time[idxT], dT$status[idxT]), error = function(e) NULL)
    piC_b[b] <- if (!is.null(fbC)) fbC$pi else NA_real_
    piT_b[b] <- if (!is.null(fbT)) fbT$pi else NA_real_
    dlt_b[b] <- if (!is.null(fbC) && !is.null(fbT)) fbT$pi - fbC$pi else NA_real_
  }
  # Each arm uses only its own converged replicates for its CI
  okC <- !is.na(piC_b)
  okT <- !is.na(piT_b)
  okD <- okC & okT          # delta requires both
  
  list(cutoff = cutoff,
       n_evC  = sum(dC$status), n_evT = sum(dT$status),
       piC    = fC$pi,  piT   = fT$pi,  delta  = fT$pi - fC$pi,
       piC_lo = quantile(piC_b[okC], 0.025),
       piC_hi = quantile(piC_b[okC], 0.975),
       piT_lo = quantile(piT_b[okT], 0.025),
       piT_hi = quantile(piT_b[okT], 0.975),
       dlt_lo = quantile(dlt_b[okD], 0.025),
       dlt_hi = quantile(dlt_b[okD], 0.975),
       n_valid = sum(okD))
}

cutoff_vals   <- c(30, 36, 42, 48, Inf)
cutoff_labels <- c("30 months", "36 months", "42 months", "48 months", "Full data")

trunc_res <- setNames(
  lapply(seq_along(cutoff_vals), function(i) {
    cat(sprintf("  Truncating at %s...\n", cutoff_labels[i]))
    truncate_fit(df, cutoff_vals[i], B_boot = 300, seed = 100 + i)
  }),
  cutoff_labels
)
trunc_res <- Filter(Negate(is.null), trunc_res)   # drop any failed fits

# Fix #2: derive all plot dimensions from the *actual* surviving entries
active_labels <- names(trunc_res)
label_map     <- c("30 months" = "30 mo", "36 months" = "36 mo",
                   "42 months" = "42 mo", "48 months" = "48 mo",
                   "Full data" = "Full")
x_labs_short  <- label_map[active_labels]
is_full_vec   <- grepl("Full", active_labels)

# Summary table
trunc_table <- data.frame(
  Truncation_point    = active_labels,
  Events_Control      = sapply(trunc_res, `[[`, "n_evC"),
  Events_Durvalumab   = sapply(trunc_res, `[[`, "n_evT"),
  pi_Control_pct      = round(sapply(trunc_res, `[[`, "piC")   * 100, 1),
  pi_Control_95CI     = sapply(trunc_res, function(x)
    sprintf("%.1f\u2013%.1f", x$piC_lo*100, x$piC_hi*100)),
  pi_Durvalumab_pct   = round(sapply(trunc_res, `[[`, "piT")   * 100, 1),
  pi_Durvalumab_95CI  = sapply(trunc_res, function(x)
    sprintf("%.1f\u2013%.1f", x$piT_lo*100, x$piT_hi*100)),
  Delta_pi_pp         = round(sapply(trunc_res, `[[`, "delta") * 100, 1),
  Delta_pi_95CI       = sapply(trunc_res, function(x)
    sprintf("%s to %s",
            pp_label(round(x$dlt_lo*100,1)),
            pp_label(round(x$dlt_hi*100,1))))
)
write.csv(trunc_table,
          file.path(subdirs[2], "Table_truncation_results.csv"),
          row.names = FALSE)

# ---- Plot data frames (derived from surviving entries only) -----------------
x_fac <- factor(x_labs_short, levels = x_labs_short)

trunc_df <- data.frame(
  xlb     = x_fac,
  delta   = sapply(trunc_res, `[[`, "delta") * 100,
  dlt_lo  = sapply(trunc_res, `[[`, "dlt_lo") * 100,
  dlt_hi  = sapply(trunc_res, `[[`, "dlt_hi") * 100,
  is_full = is_full_vec
)

trunc_arm_df <- rbind(
  data.frame(xlb   = x_fac,
             pi    = sapply(trunc_res, `[[`, "piC") * 100,
             pi_lo = sapply(trunc_res, `[[`, "piC_lo") * 100,
             pi_hi = sapply(trunc_res, `[[`, "piC_hi") * 100,
             arm   = "Control"),
  data.frame(xlb   = x_fac,
             pi    = sapply(trunc_res, `[[`, "piT") * 100,
             pi_lo = sapply(trunc_res, `[[`, "piT_lo") * 100,
             pi_hi = sapply(trunc_res, `[[`, "piT_hi") * 100,
             arm   = "Durvalumab")
)
trunc_arm_df$arm <- factor(trunc_arm_df$arm, levels = c("Control", "Durvalumab"))

fig2a <- ggplot(trunc_df, aes(x = xlb, y = delta, group = 1)) +
  geom_ribbon(aes(ymin = dlt_lo, ymax = dlt_hi), fill = "gray80", alpha = 0.55) +
  geom_line(linewidth = 1.1, colour = "black") +
  geom_point(aes(colour = is_full), size = 4.5) +
  scale_colour_manual(values = c("FALSE" = "gray30", "TRUE" = COL_DUR),
                      guide = "none") +
  # Fix #5: signed label
  scale_y_continuous(labels = function(x) pp_label(round(x, 0))) +
  labs(title    = "\u0394\u03c0 by Truncation Point",
       subtitle = "Shaded band = 95% bootstrap CI",
       x = "Data truncation point",
       y = "Cure fraction difference, \u0394\u03c0") +
  theme_pub()

fig2b <- ggplot(trunc_arm_df,
                aes(x = xlb, y = pi, colour = arm, fill = arm, group = arm)) +
  geom_ribbon(aes(ymin = pi_lo, ymax = pi_hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 3.5) +
  scale_colour_manual(values = c("Control" = COL_CTL, "Durvalumab" = COL_DUR),
                      labels = c("Control (GC)", "Durvalumab + GC")) +
  scale_fill_manual(values   = c("Control" = COL_CTL, "Durvalumab" = COL_DUR),
                    labels   = c("Control (GC)", "Durvalumab + GC")) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(title    = "Cure Fraction (\u03c0) by Truncation Point",
       subtitle = "Shaded band = 95% bootstrap CI",
       x = "Data truncation point",
       y = "Estimated cure fraction (%)") +
  theme_pub()

fig2 <- gridExtra::arrangeGrob(
  fig2a, fig2b, ncol = 2,
  top = grid::textGrob("Follow-up Truncation Sensitivity Analysis",
                       gp = grid::gpar(fontsize = 14, fontface = "bold"))
)
ggsave(file.path(subdirs[2], "Fig_truncation_sensitivity.png"),
       fig2, width = 12, height = 6, dpi = 300, bg = "white")
cat("  >>> Saved: 02_followup_truncation/\n\n")


# =============================================================================
# 5.  ANALYSIS 3 — Empirical Plateau Assessment (Maller-Zhou approach)
#     R1 Major Issue #2 (part b)
#
#     NOTE: This section implements an empirical approach inspired by
#     Maller & Zhou (Biometrics 1995), assessing whether the KM estimator
#     shows a stable non-zero plateau and whether censored observations
#     accumulate beyond the last event time.  It is NOT the exact parametric
#     test described in that paper; the formal statistic requires distributional
#     assumptions not applicable here.  Results are reported descriptively.
# =============================================================================

cat("--- Analysis 3: Empirical Plateau Assessment (Maller-Zhou approach) ---\n")

plateau_assessment <- function(t, status, arm_label, B = 1000, seed = 7) {
  t_last_ev    <- max(t[status == 1L])
  n_tail_cens  <- sum(t > t_last_ev & status == 0L)
  n_tail_total <- sum(t > t_last_ev)
  
  km0        <- survfit(Surv(t, status) ~ 1)
  S_at_ev    <- summary(km0, times = t_last_ev, extend = TRUE)$surv
  S_at_ev    <- S_at_ev[length(S_at_ev)]
  S_plateau  <- summary(km0, times = max(t),    extend = TRUE)$surv
  S_plateau  <- S_plateau[length(S_plateau)]
  
  # Bootstrap CI for KM plateau (value of KM at the end of follow-up)
  set.seed(seed)
  n      <- length(t)
  plat_b <- numeric(B)
  for (b in seq_len(B)) {
    idx_b  <- sample(n, n, replace = TRUE)
    km_b   <- tryCatch(survfit(Surv(t[idx_b], status[idx_b]) ~ 1),
                       error = function(e) NULL)
    if (is.null(km_b)) { plat_b[b] <- NA_real_; next }
    sm_b   <- summary(km_b, times = max(t[idx_b]), extend = TRUE)$surv
    plat_b[b] <- sm_b[length(sm_b)]
  }
  ok <- !is.na(plat_b)
  ci <- quantile(plat_b[ok], c(0.025, 0.975))
  
  cat(sprintf(
    "  [%s]  last event = %.1f mo | S(last event) = %.3f | ",
    arm_label, t_last_ev, S_at_ev))
  cat(sprintf(
    "censored beyond = %d | KM plateau = %.3f (95%% CI: %.3f\u2013%.3f)\n",
    n_tail_cens, S_plateau, ci[1], ci[2]))
  
  # Fix #7: report numbers; remove arbitrary threshold-based verdict
  data.frame(
    Arm                              = arm_label,
    Last_event_time_months           = round(t_last_ev, 1),
    KM_S_at_last_event               = round(S_at_ev,   4),
    N_censored_beyond_last_event     = n_tail_cens,
    N_total_obs_beyond_last_event    = n_tail_total,
    KM_plateau_end_of_followup       = round(S_plateau, 4),
    KM_plateau_95CI_lower            = round(ci[1],     4),
    KM_plateau_95CI_upper            = round(ci[2],     4),
    Prop_bootstrap_plateau_below_001 = round(mean(plat_b[ok] < 0.01), 4)
  )
}

mz_C <- plateau_assessment(dfC$time, dfC$status, "Control",    B = 1000, seed =  7)
mz_T <- plateau_assessment(dfT$time, dfT$status, "Durvalumab", B = 1000, seed = 13)
mz_table <- rbind(mz_C, mz_T)
write.csv(mz_table,
          file.path(subdirs[3], "Table_maller_zhou.csv"),
          row.names = FALSE)

# ---- Figure ------------------------------------------------------------------
km_full   <- survfit(Surv(time, status) ~ arm, data = df)
km_sum    <- summary(km_full)
km_df_all <- data.frame(
  time = km_sum$time,
  surv = km_sum$surv,
  arm  = sub("arm=", "", as.character(km_sum$strata))
)
km_df_all$arm <- factor(km_df_all$arm, levels = c("Control", "Durvalumab"))

last_ev_C <- mz_C$Last_event_time_months
last_ev_T <- mz_T$Last_event_time_months

# Position text labels so they do not overlap regardless of data
txt_y_C <- 8
txt_y_T <- max(txt_y_C + 10, mz_T$KM_plateau_end_of_followup * 100 * 0.2)

fig3 <- ggplot(km_df_all, aes(x = time, y = surv * 100, colour = arm)) +
  annotate("rect", xmin = last_ev_C, xmax = max(dfC$time),
           ymin = 0, ymax = 100, alpha = 0.07, fill = COL_CTL) +
  annotate("rect", xmin = last_ev_T, xmax = max(dfT$time),
           ymin = 0, ymax = 100, alpha = 0.07, fill = COL_DUR) +
  geom_step(linewidth = 1.0) +
  geom_vline(xintercept = last_ev_C, colour = COL_CTL,
             linetype = "dashed", linewidth = 0.75, alpha = 0.85) +
  geom_vline(xintercept = last_ev_T, colour = COL_DUR,
             linetype = "dashed", linewidth = 0.75, alpha = 0.85) +
  annotate("text", x = last_ev_C - 0.5, y = txt_y_C,
           label = sprintf("Last event\n%.1f mo", last_ev_C),
           colour = COL_CTL, size = 3.5, hjust = 1, fontface = "italic") +
  annotate("text", x = last_ev_T + 0.5, y = txt_y_T,
           label = sprintf("Last event\n%.1f mo", last_ev_T),
           colour = COL_DUR, size = 3.5, hjust = 0, fontface = "italic") +
  scale_colour_manual(
    values = c("Control" = COL_CTL, "Durvalumab" = COL_DUR),
    labels = c("Control (GC alone)", "Durvalumab + GC")) +
  scale_x_continuous(breaks = seq(0, 60, 12), limits = c(0, 66)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20),
                     labels = function(x) paste0(x, "%")) +
  labs(title    = "Empirical Plateau Assessment (Maller-Zhou approach)",
       subtitle  = paste0("Shaded area = region beyond last observed event per arm\n",
                          "Dashed line = last event time per arm"),
       x = "Time since randomisation (months)",
       y = "Event-free survival (%)") +
  theme_pub()

ggsave(file.path(subdirs[3], "Fig_maller_zhou_KM.png"),
       fig3, width = 9, height = 6, dpi = 300, bg = "white")
cat("  >>> Saved: 03_maller_zhou/\n\n")


# =============================================================================
# 6.  ANALYSIS 4 — Akaike Weights for Cure-Regression Models
#     R1 Major Issue #3
# =============================================================================

cat("--- Analysis 4: Akaike Weights ---\n")

# NLL — log-normal cure regression
# Parameter convention:
#   "cure_only"    → [beta0, beta1, alpha0,         log_sigma]  (4 par)
#   "latency_only" → [beta0,        alpha0,  alpha1, log_sigma]  (4 par)
#   "both"         → [beta0, beta1, alpha0,  alpha1, log_sigma]  (5 par)

nll_cure_reg <- function(par, t, status, x, model) {
  t <- pmax(as.numeric(t), 1e-12)
  if (model == "cure_only") {
    beta0 <- par[1]; beta1 <- par[2]
    alpha0 <- par[3]; alpha1 <- 0; log_sigma <- par[4]
  } else if (model == "latency_only") {
    beta0 <- par[1]; beta1 <- 0
    alpha0 <- par[2]; alpha1 <- par[3]; log_sigma <- par[4]
  } else {
    beta0 <- par[1]; beta1 <- par[2]
    alpha0 <- par[3]; alpha1 <- par[4]; log_sigma <- par[5]
  }
  pi    <- pmin(pmax(inv_logit(beta0 + beta1 * x), 1e-9), 1 - 1e-9)
  mu    <- alpha0 + alpha1 * x
  sigma <- exp(log_sigma)
  t_    <- pmax(t, 1e-12)
  Su    <- 1 - plnorm(t_, meanlog = mu, sdlog = sigma)
  fu    <- dlnorm(t_,    meanlog = mu, sdlog = sigma)
  lik   <- ifelse(status == 1L, (1 - pi) * fu, pi + (1 - pi) * Su)
  -sum(log(pmax(lik, 1e-300)))
}

fit_reg_model <- function(model, t, status, x, b0, b1, a0, a1, sg) {
  if (model == "cure_only") {
    init  <- c(b0, b1, a0, sg)
    lower <- c(-10, -10, -10, log(0.05)); upper <- c(10, 10, 10, log(5))
    n_par <- 4L
  } else if (model == "latency_only") {
    init  <- c(b0, a0, a1, sg)
    lower <- c(-10, -10, -10, log(0.05)); upper <- c(10, 10, 10, log(5))
    n_par <- 4L
  } else {
    init  <- c(b0, b1, a0, a1, sg)
    lower <- c(-10, -10, -10, -10, log(0.05)); upper <- c(10, 10, 10, 10, log(5))
    n_par <- 5L
  }
  opt <- tryCatch(
    optim(par = init, fn = nll_cure_reg,
          t = t, status = status, x = x, model = model,
          method = "L-BFGS-B", lower = lower, upper = upper),
    error = function(e) NULL
  )
  if (is.null(opt)) return(NULL)
  loglik <- -opt$value
  aic    <- -2 * loglik + 2 * n_par
  if (model == "latency_only") {
    piC <- inv_logit(opt$par[1]); piT <- piC
  } else {
    piC <- inv_logit(opt$par[1]); piT <- inv_logit(opt$par[1] + opt$par[2])
  }
  list(model = model, n_par = n_par, loglik = loglik, aic = aic,
       piC = piC, piT = piT, delta_pi = piT - piC)
}

t_all <- df$time; s_all <- df$status
x_all <- as.integer(df$arm == "Durvalumab")

b0_i <- ref_C$a;              b1_i <- ref_T$a - ref_C$a
a0_i <- ref_C$mu;             a1_i <- ref_T$mu - ref_C$mu
sg_i <- log(mean(c(ref_C$sigma, ref_T$sigma)))

cat("  Fitting cure-only model...\n")
r_cure <- fit_reg_model("cure_only",    t_all, s_all, x_all, b0_i,b1_i,a0_i,a1_i,sg_i)
cat("  Fitting latency-only model...\n")
r_lat  <- fit_reg_model("latency_only", t_all, s_all, x_all, b0_i,b1_i,a0_i,a1_i,sg_i)
cat("  Fitting cure+latency model...\n")
r_both <- fit_reg_model("both",         t_all, s_all, x_all, b0_i,b1_i,a0_i,a1_i,sg_i)

# Fix #1: build aic_table from whatever actually converged (no fixed 3-row assumption)
model_display <- c(cure_only    = "Cure effect only",
                   latency_only = "Latency effect only",
                   both         = "Cure + Latency (both)")

models_ok <- Filter(Negate(is.null), list(r_cure, r_lat, r_both))
if (length(models_ok) == 0L) stop("All cure-regression models failed to converge.")

aic_vals  <- sapply(models_ok, `[[`, "aic")
delta_aic <- aic_vals - min(aic_vals)
w_aic     <- exp(-0.5 * delta_aic) / sum(exp(-0.5 * delta_aic))

aic_table <- data.frame(
  Model             = model_display[sapply(models_ok, `[[`, "model")],
  N_parameters      = sapply(models_ok, `[[`, "n_par"),
  Log_likelihood    = round(sapply(models_ok, `[[`, "loglik"), 2),
  AIC               = round(aic_vals, 1),
  Delta_AIC         = round(delta_aic, 1),
  Akaike_weight     = round(w_aic, 4),
  Akaike_weight_pct = paste0(round(w_aic * 100, 1), "%"),
  pi_Control_pct    = round(sapply(models_ok, `[[`, "piC") * 100, 1),
  pi_Durvalumab_pct = round(sapply(models_ok, `[[`, "piT") * 100, 1),
  Delta_pi_pp       = round(sapply(models_ok, `[[`, "delta_pi") * 100, 1),
  row.names         = NULL
)
write.csv(aic_table,
          file.path(subdirs[4], "Table_akaike_weights.csv"),
          row.names = FALSE)

cat("\n  Akaike Weights:\n")
print(aic_table[, c("Model", "AIC", "Delta_AIC", "Akaike_weight_pct")])

# ---- Figure ------------------------------------------------------------------
aic_table$Model_f <- factor(aic_table$Model,
                            levels = rev(c("Cure effect only",
                                           "Cure + Latency (both)",
                                           "Latency effect only")))
aic_table$is_best <- aic_table$Delta_AIC == 0

fig4 <- ggplot(aic_table,
               aes(x = Model_f, y = Akaike_weight * 100, fill = is_best)) +
  geom_col(colour = "black", linewidth = 0.4, width = 0.55) +
  geom_text(
    aes(label = sprintf("%s\n(AIC = %.1f; \u0394AIC = +%.1f)",
                        Akaike_weight_pct, AIC, Delta_AIC)),
    vjust = -0.3, size = 4.2, fontface = "bold", lineheight = 0.95) +
  scale_fill_manual(values = c("TRUE" = COL_DUR, "FALSE" = "gray70"),
                    guide = "none") +
  scale_y_continuous(limits = c(0, 105),
                     labels = function(x) paste0(x, "%"),
                     expand = c(0, 0)) +
  labs(title    = "Akaike Weights: Cure-Regression Model Comparison",
       subtitle  = "Higher weight reflects greater relative empirical support",
       x = NULL, y = "Akaike weight (%)") +
  theme_pub() +
  theme(axis.text.x = element_text(size = 12))

ggsave(file.path(subdirs[4], "Fig_akaike_weights.png"),
       fig4, width = 8, height = 6, dpi = 300, bg = "white")
cat("  >>> Saved: 04_akaike_weights/\n\n")


# =============================================================================
# 7.  ANALYSIS 5 — Conditional Survival Analysis
#     R1 Suggested Additional Analysis
# =============================================================================

cat("--- Analysis 5: Conditional Survival Analysis ---\n")
cat("  CS(2 yr | s yr) = S(s + 24 mo) / S(s)  via Kaplan-Meier\n")

s_months  <- c(0, 12, 24, 36, 48)
t_add_mo  <- 24

km_C <- survfit(Surv(time, status) ~ 1, data = dfC)
km_T <- survfit(Surv(time, status) ~ 1, data = dfT)

km_surv_at <- function(km_fit, t_q) {
  sm <- summary(km_fit, times = t_q, extend = TRUE)$surv
  sm[length(sm)]
}

cs_point <- function(km_fit, s, t_add = 24) {
  S_s    <- km_surv_at(km_fit, s)
  S_plus <- km_surv_at(km_fit, s + t_add)
  if (is.na(S_s) || S_s <= 0) return(NA_real_)
  S_plus / S_s
}

cs_C_pt <- sapply(s_months, cs_point, km_fit = km_C, t_add = t_add_mo)
cs_T_pt <- sapply(s_months, cs_point, km_fit = km_T, t_add = t_add_mo)

# Bootstrap CI — Fix #8 style: separate per arm
cs_boot_arm <- function(dat, s_mo_vec, t_add = 24, B = 500, seed = 77) {
  set.seed(seed)
  n   <- nrow(dat)
  mat <- matrix(NA_real_, nrow = B, ncol = length(s_mo_vec))
  for (b in seq_len(B)) {
    idx  <- sample(n, n, replace = TRUE)
    km_b <- tryCatch(survfit(Surv(time, status) ~ 1, data = dat[idx, ]),
                     error = function(e) NULL)
    if (is.null(km_b)) next
    for (j in seq_along(s_mo_vec)) {
      mat[b, j] <- tryCatch(cs_point(km_b, s_mo_vec[j], t_add),
                            error = function(e) NA_real_)
    }
  }
  mat
}

cat("  Bootstrap Control    (B = 500)...\n")
boot_C_mat <- cs_boot_arm(dfC, s_months, t_add_mo, B = 500, seed = 77)
cat("  Bootstrap Durvalumab (B = 500)...\n")
boot_T_mat <- cs_boot_arm(dfT, s_months, t_add_mo, B = 500, seed = 78)

ci_boot <- function(mat, prob) apply(mat, 2, quantile, probs = prob, na.rm = TRUE)

cs_C_lo <- ci_boot(boot_C_mat, 0.025)
cs_C_hi <- ci_boot(boot_C_mat, 0.975)
cs_T_lo <- ci_boot(boot_T_mat, 0.025)
cs_T_hi <- ci_boot(boot_T_mat, 0.975)

# Fix #3: accurate description of KM extension (not parametric extrapolation)
max_obs <- min(max(dfC$time), max(dfT$time))   # conservative: shorter of the two
cs_note <- ifelse(
  (s_months + t_add_mo) > max_obs,
  "KM extended beyond last observed event (step function held constant at plateau)",
  "Within observed follow-up"
)

cs_table <- data.frame(
  Conditioning_s_years  = s_months / 12,
  CS_Control_pct        = round(cs_C_pt * 100, 1),
  CS_Control_95CI       = sprintf("%.1f\u2013%.1f", cs_C_lo*100, cs_C_hi*100),
  CS_Durvalumab_pct     = round(cs_T_pt * 100, 1),
  CS_Durvalumab_95CI    = sprintf("%.1f\u2013%.1f", cs_T_lo*100, cs_T_hi*100),
  Difference_pp         = round((cs_T_pt - cs_C_pt) * 100, 1),
  Note                  = cs_note
)
write.csv(cs_table,
          file.path(subdirs[5], "Table_conditional_survival.csv"),
          row.names = FALSE)

cat("\n  CS(2 yr | s yr):\n")
print(cs_table[, c("Conditioning_s_years", "CS_Control_pct",
                   "CS_Durvalumab_pct", "Difference_pp", "Note")])

# ---- Figure ------------------------------------------------------------------
cs_df <- rbind(
  data.frame(s_yr = s_months/12, cs = cs_C_pt*100,
             lo = cs_C_lo*100, hi = cs_C_hi*100, arm = "Control"),
  data.frame(s_yr = s_months/12, cs = cs_T_pt*100,
             lo = cs_T_lo*100, hi = cs_T_hi*100, arm = "Durvalumab")
)
cs_df$arm      <- factor(cs_df$arm, levels = c("Control", "Durvalumab"))
cs_df$km_ext   <- (cs_df$s_yr * 12 + t_add_mo) > max_obs   # KM extension flag

# Fix #4: dynamic y-axis limits
y_lo_raw <- floor(min(cs_df$lo, na.rm = TRUE) / 5) * 5
y_lo     <- max(0, y_lo_raw - 5)
y_hi     <- 100

fig5 <- ggplot(cs_df, aes(x = s_yr, y = cs, colour = arm, fill = arm)) +
  # CI ribbons — solid for within-FU, lighter for KM-extended region
  geom_ribbon(data = subset(cs_df, !km_ext),
              aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_ribbon(data = subset(cs_df, km_ext),
              aes(ymin = lo, ymax = hi), alpha = 0.07, colour = NA) +
  geom_line(data = subset(cs_df, !km_ext), linewidth = 1.3) +
  geom_line(data = subset(cs_df, km_ext),  linewidth = 1.3, linetype = "dashed") +
  geom_point(aes(shape = km_ext), size = 4.5) +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 1),
                     labels = c("Within observed follow-up",
                                "KM extended beyond last observed event"),
                     name = NULL) +
  scale_colour_manual(
    values = c("Control" = COL_CTL, "Durvalumab" = COL_DUR),
    labels = c("Control (GC alone)", "Durvalumab + GC")) +
  scale_fill_manual(
    values = c("Control" = COL_CTL, "Durvalumab" = COL_DUR),
    labels = c("Control (GC alone)", "Durvalumab + GC")) +
  scale_x_continuous(breaks = s_months / 12,
                     labels = paste0(s_months / 12, " yr")) +
  # Fix #4: dynamic limits
  scale_y_continuous(limits = c(y_lo, y_hi),
                     breaks = seq(y_lo, y_hi, 10),
                     labels = function(x) paste0(x, "%")) +
  labs(title    = "Conditional Event-free Survival",
       subtitle  = paste0(
         "Probability of remaining event-free for an additional 2 years\n",
         "given already event-free at s years | Shaded band = 95% bootstrap CI\n",
         "Open circle / dashed = KM step function extended beyond last observed event"),
       x = "Years already event-free (s)",
       y = "CS(2 yr | s) (%)") +
  theme_pub() +
  theme(legend.position = "bottom", legend.box = "vertical")

ggsave(file.path(subdirs[5], "Fig_conditional_survival.png"),
       fig5, width = 8.5, height = 6.5, dpi = 300, bg = "white")
cat("  >>> Saved: 05_conditional_survival/\n\n")


# =============================================================================
# 8.  Final Summary
# =============================================================================

cat("=================================================================\n")
cat("             REVISION ANALYSES COMPLETE  (v2)                    \n")
cat("=================================================================\n\n")
cat(sprintf("Output folder: %s/\n\n", root_dir))

cat("  01_reconstruction_uncertainty/\n")
cat(sprintf("     SD=0.5 mo: mean \u0394\u03c0 = %s,  SD = %.1f pp\n",
            pp_label(round(mean(pert05$delta)*100,1)), sd(pert05$delta)*100))
cat(sprintf("     SD=1.0 mo: mean \u0394\u03c0 = %s,  SD = %.1f pp\n\n",
            pp_label(round(mean(pert10$delta)*100,1)), sd(pert10$delta)*100))

cat("  02_followup_truncation/\n")
for (nm in names(trunc_res)) {
  r <- trunc_res[[nm]]
  cat(sprintf("     %-12s : \u0394\u03c0 = %s (95%% CI %s to %s)\n",
              nm, pp_label(round(r$delta*100,1)),
              pp_label(round(r$dlt_lo*100,1)), pp_label(round(r$dlt_hi*100,1))))
}

cat("\n  03_maller_zhou/\n")
cat(sprintf("     Control    : last event %.1f mo, n_censored_tail=%d, plateau=%.3f\n",
            mz_C$Last_event_time_months, mz_C$N_censored_beyond_last_event,
            mz_C$KM_plateau_end_of_followup))
cat(sprintf("     Durvalumab : last event %.1f mo, n_censored_tail=%d, plateau=%.3f\n",
            mz_T$Last_event_time_months, mz_T$N_censored_beyond_last_event,
            mz_T$KM_plateau_end_of_followup))

cat("\n  04_akaike_weights/\n")
for (i in seq_len(nrow(aic_table))) {
  cat(sprintf("     %-28s AIC = %.1f,  weight = %s\n",
              aic_table$Model[i], aic_table$AIC[i], aic_table$Akaike_weight_pct[i]))
}

cat("\n  05_conditional_survival/\n")
for (i in seq_along(s_months)) {
  ext <- if (cs_note[i] != "Within observed follow-up") "  [KM extended]" else ""
  cat(sprintf("     CS(2 yr | s=%d yr):  Control %.1f%%,  Durvalumab %.1f%%  (diff %s)%s\n",
              s_months[i]/12, cs_C_pt[i]*100, cs_T_pt[i]*100,
              pp_label(round((cs_T_pt[i]-cs_C_pt[i])*100,1)), ext))
}

cat("\n=================================================================\n")
cat("NOTE:  All cure-model analyses use the log-normal distribution\n")
cat("       (best fit by AIC in the original analysis).\n")
cat("=================================================================\n\n")