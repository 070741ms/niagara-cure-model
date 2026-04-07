# =============================================================================
# 01_main_analysis.R
#
# Main analysis: Kaplan-Meier, hazard ratio, mixture cure model,
# bootstrap confidence intervals, and cure regression.
#
# Input:  efs_ripd_excel.csv
# Output: outputs/tables/  outputs/figures/
# =============================================================================

required_pkgs <- c("survival", "ggplot2", "gridExtra", "scales")
for (p in required_pkgs) {
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cloud.r-project.org")
  library(p, character.only = TRUE)
}

dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/tables",  recursive = TRUE, showWarnings = FALSE)

COL_CTL  <- "#4A90D9"
COL_DURV <- "#F5A623"

theme_pub <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(plot.background  = element_rect(fill = "white", colour = NA),
          panel.background = element_rect(fill = "white", colour = NA),
          plot.title       = element_text(face = "bold", hjust = 0.5),
          axis.title       = element_text(face = "bold"),
          axis.text        = element_text(colour = "black"),
          legend.background = element_rect(fill = "white", colour = NA),
          legend.key       = element_rect(fill = "white", colour = NA),
          plot.margin      = margin(12, 12, 12, 12))
}


# =============================================================================
# 1. Data
# =============================================================================
df <- read.csv("efs_ripd_excel.csv", stringsAsFactors = FALSE)
df$time   <- as.numeric(df$time)
df$status <- as.integer(df$status)
df$arm    <- factor(df$arm, levels = c("Control", "Durvalumab"))

df$time <- pmax(df$time, 0.001)   # guard against t<=0 (consistent with reconstruction)
stopifnot(
  all(df$status %in% c(0L, 1L)),
  all(levels(df$arm) %in% c("Control", "Durvalumab")),
  all(df$time > 0)
)

dfC <- df[df$arm == "Control",    ]
dfT <- df[df$arm == "Durvalumab", ]

cat(sprintf("Data loaded: N=%d (Control=%d, Durvalumab=%d)\n\n",
            nrow(df), nrow(dfC), nrow(dfT)))


# =============================================================================
# 2. Kaplan-Meier and Cox regression
# =============================================================================
km  <- survfit(Surv(time, status) ~ arm, data = df)
cox <- coxph(Surv(time, status) ~ arm, data = df)
lrt <- survdiff(Surv(time, status) ~ arm, data = df)

hr    <- exp(coef(cox)["armDurvalumab"])
hr_ci <- exp(confint(cox)["armDurvalumab", ])
lrp   <- 1 - pchisq(lrt$chisq, df = 1L)

sm24  <- summary(km, times = 24)
efs24 <- setNames(sm24$surv, sub("arm=", "", as.character(sm24$strata)))
med   <- setNames(summary(km)$table[, "median"],
                  sub("arm=", "", rownames(summary(km)$table)))

cat(sprintf("HR: %.2f (95%% CI: %.2f-%.2f), log-rank p%s\n",
            hr, hr_ci[1], hr_ci[2],
            ifelse(lrp < 0.001, "<0.001", sprintf("=%.4f", lrp))))

# KM summary table
km_tbl <- data.frame(
  Arm       = c("Control", "Durvalumab"),
  N         = c(nrow(dfC), nrow(dfT)),
  Events    = c(sum(dfC$status), sum(dfT$status)),
  EFS_24mo  = round(efs24[c("Control", "Durvalumab")] * 100, 1),
  Median_EFS = c(round(med["Control"], 1),
                 ifelse(is.na(med["Durvalumab"]) | is.infinite(med["Durvalumab"]),
                        NA, round(med["Durvalumab"], 1)))
)
write.csv(km_tbl, "outputs/tables/Table_KM.csv", row.names = FALSE)

# KM plot (base Figure 1 -- cure model overlay added after model fitting)
km_sm  <- summary(km)
km_df  <- data.frame(
  time = km_sm$time,
  surv = km_sm$surv,
  arm  = sub("arm=", "", as.character(km_sm$strata))
)
km_df$arm <- factor(km_df$arm, levels = c("Control", "Durvalumab"))

# Prepend time=0 for each arm
km_df <- rbind(
  data.frame(time = 0, surv = 1, arm = "Control"),
  data.frame(time = 0, surv = 1, arm = "Durvalumab"),
  km_df
)
km_df$arm <- factor(km_df$arm, levels = c("Control", "Durvalumab"))


# =============================================================================
# 3. Mixture cure model
# =============================================================================

inv_logit <- function(z) 1 / (1 + exp(-z))

# Negative log-likelihood functions
nll_weibull <- function(par, t, ev) {
  pi    <- pmin(pmax(inv_logit(par[1]), 1e-9), 1 - 1e-9)
  shape <- exp(par[2]); scale <- exp(par[3]); t <- pmax(t, 1e-12)
  Su <- exp(-(t / scale)^shape)
  fu <- (shape / scale) * (t / scale)^(shape - 1) * Su
  lik <- ifelse(ev, (1 - pi) * fu, pi + (1 - pi) * Su)
  -sum(log(pmax(lik, 1e-300)))
}

nll_loglogistic <- function(par, t, ev) {
  pi    <- pmin(pmax(inv_logit(par[1]), 1e-9), 1 - 1e-9)
  shape <- exp(par[2]); scale <- exp(par[3]); t <- pmax(t, 1e-12)
  Su <- 1 / (1 + (t / scale)^shape)
  fu <- (shape / scale) * (t / scale)^(shape - 1) / (1 + (t / scale)^shape)^2
  lik <- ifelse(ev, (1 - pi) * fu, pi + (1 - pi) * Su)
  -sum(log(pmax(lik, 1e-300)))
}

nll_lognormal <- function(par, t, ev) {
  pi    <- pmin(pmax(inv_logit(par[1]), 1e-9), 1 - 1e-9)
  mu    <- par[2]; sigma <- exp(par[3]); t <- pmax(t, 1e-12)
  Su <- 1 - plnorm(t, mu, sigma)
  fu <- dlnorm(t, mu, sigma)
  lik <- ifelse(ev, (1 - pi) * fu, pi + (1 - pi) * Su)
  -sum(log(pmax(lik, 1e-300)))
}

fit_cure <- function(t, ev, dist) {
  med_t <- max(median(t[ev == 1L], na.rm = TRUE), 0.1)
  if (dist == "weibull") {
    fn  <- nll_weibull
    par <- c(0, log(1), log(med_t))
    lo  <- c(-10, log(0.1), log(0.1))
    hi  <- c( 10, log(10),  log(200))
  } else if (dist == "loglogistic") {
    fn  <- nll_loglogistic
    par <- c(0, log(2), log(med_t))
    lo  <- c(-10, log(0.1), log(0.1))
    hi  <- c( 10, log(10),  log(200))
  } else {
    fn  <- nll_lognormal
    par <- c(0, log(med_t), log(1))
    lo  <- c(-10, -5, log(0.05))
    hi  <- c( 10, 10, log(5))
  }
  opt <- tryCatch(
    optim(par, fn, t = t, ev = ev, method = "L-BFGS-B",
          lower = lo, upper = hi),
    error = function(e) NULL)
  if (is.null(opt)) return(NULL)
  pi  <- inv_logit(opt$par[1])
  ll  <- -opt$value
  list(dist = dist, pi = pi, par = opt$par,
       aic_arm = -2 * ll + 2 * 3)   # per-arm AIC (3 params)
}

# Survival prediction
S_pred <- function(t, fit) {
  pi <- fit$pi
  p  <- fit$par
  Su <- switch(fit$dist,
    weibull     = exp(-(t / exp(p[3]))^exp(p[2])),
    loglogistic = 1 / (1 + (t / exp(p[3]))^exp(p[2])),
    lognormal   = 1 - plnorm(pmax(t, 1e-12), p[2], exp(p[3]))
  )
  pi + (1 - pi) * Su
}

cat("\nFitting mixture cure models (3 distributions)...\n")
dists <- c("weibull", "loglogistic", "lognormal")
fits  <- list()

for (dist in dists) {
  fC <- fit_cure(dfC$time, dfC$status, dist)
  fT <- fit_cure(dfT$time, dfT$status, dist)
  if (!is.null(fC) && !is.null(fT)) {
    fits[[dist]] <- list(C = fC, T = fT,
                         delta = fT$pi - fC$pi,
                         aic   = fC$aic_arm + fT$aic_arm)
    cat(sprintf("  %-12s  pi_Control=%.3f  pi_Durv=%.3f  delta=+%.3f  AIC=%.1f\n",
                dist, fC$pi, fT$pi, fT$pi - fC$pi,
                fC$aic_arm + fT$aic_arm))
  }
}

best_dist <- names(which.min(sapply(fits, `[[`, "aic")))
best      <- fits[[best_dist]]
cat(sprintf("\nBest distribution (AIC): %s\n", best_dist))

# Model selection table
sel_tbl <- do.call(rbind, lapply(names(fits), function(d) {
  f <- fits[[d]]
  data.frame(Distribution   = d,
             pi_Control_pct = round(f$C$pi * 100, 1),
             pi_Durv_pct    = round(f$T$pi * 100, 1),
             Delta_pi_pp    = round(f$delta * 100, 1),
             AIC            = round(f$aic, 1),
             Best           = d == best_dist)
}))
write.csv(sel_tbl, "outputs/tables/Table_ModelSelection.csv", row.names = FALSE)


# =============================================================================
# 4. Bootstrap CI (B = 1000)
# =============================================================================
cat("\nBootstrap CI (B=1000)...\n")

set.seed(123)
B       <- 1000L
pi_C_b  <- numeric(B); pi_T_b <- numeric(B); dlt_b <- numeric(B)

for (b in seq_len(B)) {
  idxC <- sample(nrow(dfC), replace = TRUE)
  idxT <- sample(nrow(dfT), replace = TRUE)
  fC_b <- tryCatch(fit_cure(dfC$time[idxC], dfC$status[idxC], best_dist),
                   error = function(e) NULL)
  fT_b <- tryCatch(fit_cure(dfT$time[idxT], dfT$status[idxT], best_dist),
                   error = function(e) NULL)
  pi_C_b[b] <- if (!is.null(fC_b)) fC_b$pi else NA_real_
  pi_T_b[b] <- if (!is.null(fT_b)) fT_b$pi else NA_real_
  dlt_b[b]  <- if (!is.null(fC_b) && !is.null(fT_b)) fT_b$pi - fC_b$pi else NA_real_
}

ok <- !is.na(dlt_b)
cat(sprintf("  Converged: %d / %d\n", sum(ok), B))

ci <- function(x, valid) quantile(x[valid], c(0.025, 0.975))
ci_C   <- ci(pi_C_b, !is.na(pi_C_b))
ci_T   <- ci(pi_T_b, !is.na(pi_T_b))
ci_dlt <- ci(dlt_b, ok)

cat(sprintf("  pi_Control:    %.3f (95%% CI: %.3f-%.3f)\n",
            best$C$pi, ci_C[1], ci_C[2]))
cat(sprintf("  pi_Durvalumab: %.3f (95%% CI: %.3f-%.3f)\n",
            best$T$pi, ci_T[1], ci_T[2]))
cat(sprintf("  Delta-pi:      +%.3f (95%% CI: +%.3f to +%.3f)\n",
            best$delta, ci_dlt[1], ci_dlt[2]))

cure_tbl <- data.frame(
  Group            = c("Control", "Durvalumab", "Difference (Delta-pi)"),
  CureFraction_pct = round(c(best$C$pi, best$T$pi, best$delta) * 100, 1),
  CI95_lower_pct   = round(c(ci_C[1], ci_T[1], ci_dlt[1]) * 100, 1),
  CI95_upper_pct   = round(c(ci_C[2], ci_T[2], ci_dlt[2]) * 100, 1),
  Distribution     = best_dist
)
write.csv(cure_tbl, "outputs/tables/Table_CureFraction.csv", row.names = FALSE)


# =============================================================================
# 5. Cure regression
# =============================================================================
# Cure regression latency component: log-normal (plnorm/dlnorm) throughout.
# This matches the published analysis in which log-normal was identified as the
# best-fitting distribution above. If best_dist != "lognormal" in your data,
# creg_C/creg_T supply valid log-normal initial values (par[2]=mu, par[3]=log(sigma))
# so the optimiser starts in the correct parameter space.
if (best_dist == "lognormal") {
  creg_C <- best$C;  creg_T <- best$T
} else {
  creg_C <- fit_cure(dfC$time, dfC$status, "lognormal")
  creg_T <- fit_cure(dfT$time, dfT$status, "lognormal")
  if (is.null(creg_C) || is.null(creg_T))
    stop("Log-normal fit required for cure regression initialisation failed.")
  message(sprintf(
    "Note: cure regression uses log-normal latency; best_dist = '%s' in main analysis.",
    best_dist))
}
cat("\nCure regression (log-normal latency, 3 nested models)...\n")

nll_creg <- function(par, t, ev, x, model) {
  # Latency: log-normal (mu = a0 + a1*x, sigma = exp(ls)). Intentionally fixed.
  t <- pmax(as.numeric(t), 1e-12)
  if (model == "cure_only") {
    b0 <- par[1]; b1 <- par[2]; a0 <- par[3]; a1 <- 0;      ls <- par[4]; np <- 4L
  } else if (model == "latency_only") {
    b0 <- par[1]; b1 <- 0;      a0 <- par[2]; a1 <- par[3]; ls <- par[4]; np <- 4L
  } else {
    b0 <- par[1]; b1 <- par[2]; a0 <- par[3]; a1 <- par[4]; ls <- par[5]; np <- 5L
  }
  pi <- pmin(pmax(inv_logit(b0 + b1 * x), 1e-9), 1 - 1e-9)
  mu <- a0 + a1 * x; sigma <- exp(ls)
  Su <- 1 - plnorm(t, mu, sigma); fu <- dlnorm(t, mu, sigma)
  lik <- ifelse(ev, (1 - pi) * fu, pi + (1 - pi) * Su)
  (-sum(log(pmax(lik, 1e-300))))
}

fit_creg <- function(model) {
  x   <- as.integer(df$arm == "Durvalumab")
  b0i <- qlogis(creg_C$pi); b1i <- qlogis(creg_T$pi) - b0i
  a0i <- creg_C$par[2];     a1i <- creg_T$par[2] - a0i
  lsi <- log(mean(c(exp(creg_C$par[3]), exp(creg_T$par[3]))))

  if (model == "cure_only") {
    init <- c(b0i, b1i, a0i, lsi)
    lo   <- c(-10, -10, -10, log(0.05)); hi <- c(10, 10, 10, log(5))
  } else if (model == "latency_only") {
    init <- c(b0i, a0i, a1i, lsi)
    lo   <- c(-10, -10, -10, log(0.05)); hi <- c(10, 10, 10, log(5))
  } else {
    init <- c(b0i, b1i, a0i, a1i, lsi)
    lo   <- c(-10, -10, -10, -10, log(0.05)); hi <- c(10, 10, 10, 10, log(5))
  }

  obj <- function(par) nll_creg(par, df$time, df$status, x, model)
  opt <- tryCatch(
    optim(init, obj, method = "L-BFGS-B", lower = lo, upper = hi),
    error = function(e) NULL)
  if (is.null(opt)) return(NULL)

  np  <- if (model == "both") 5L else 4L
  ll  <- -opt$value
  aic <- -2 * ll + 2 * np

  par <- opt$par
  if (model == "cure_only") {
    piC <- inv_logit(par[1]); piT <- inv_logit(par[1] + par[2])
  } else if (model == "latency_only") {
    piC <- inv_logit(par[1]); piT <- piC
  } else {
    piC <- inv_logit(par[1]); piT <- inv_logit(par[1] + par[2])
  }
  list(model = model, piC = piC, piT = piT, delta = piT - piC,
       aic = aic, np = np, ll = ll)
}

r_co  <- fit_creg("cure_only")
r_lo  <- fit_creg("latency_only")
r_bo  <- fit_creg("both")

models_ok <- Filter(Negate(is.null), list(r_co, r_lo, r_bo))
if (length(models_ok) == 0L) stop("All cure regression models failed.")

aic_v     <- sapply(models_ok, `[[`, "aic")
daic      <- aic_v - min(aic_v)
w_aic     <- exp(-0.5 * daic) / sum(exp(-0.5 * daic))

creg_tbl <- data.frame(
  Model             = sapply(models_ok, `[[`, "model"),
  pi_Control_pct    = round(sapply(models_ok, `[[`, "piC") * 100, 1),
  pi_Durvalumab_pct = round(sapply(models_ok, `[[`, "piT") * 100, 1),
  Delta_pi_pp       = round(sapply(models_ok, `[[`, "delta") * 100, 1),
  AIC               = round(aic_v, 1),
  Delta_AIC         = round(daic, 1),
  Akaike_weight_pct = round(w_aic * 100, 1)
)
write.csv(creg_tbl, "outputs/tables/Table_CureRegression.csv", row.names = FALSE)

cat("\nCure regression results:\n")
print(creg_tbl[, c("Model", "AIC", "Delta_AIC", "Akaike_weight_pct")])
if (length(daic) >= 2L)
  cat(sprintf("\n  Note: \u0394AIC between best and second-best = %.1f\n",
              sort(daic)[2L]))


# =============================================================================
# 6. Figures
# =============================================================================

# -- Figure 1: KM with cure model overlay ------------------------------------
grid_t <- seq(0.01, 65, length.out = 500)
mod_df <- rbind(
  data.frame(time = grid_t,
             surv = S_pred(grid_t, best$C),
             arm  = "Control"),
  data.frame(time = grid_t,
             surv = S_pred(grid_t, best$T),
             arm  = "Durvalumab")
)
mod_df$arm <- factor(mod_df$arm, levels = c("Control", "Durvalumab"))

fig1 <- ggplot() +
  geom_step(data = km_df,
            aes(x = time, y = surv * 100, colour = arm),
            linewidth = 1.0) +
  geom_line(data = mod_df,
            aes(x = time, y = surv * 100, colour = arm),
            linetype = "dashed", linewidth = 0.9, alpha = 0.85) +
  geom_hline(yintercept = best$C$pi * 100,
             linetype = "dotted", colour = COL_CTL, linewidth = 0.7) +
  geom_hline(yintercept = best$T$pi * 100,
             linetype = "dotted", colour = COL_DURV, linewidth = 0.7) +
  annotate("text", x = 3, y = best$C$pi * 100 + 2.5, hjust = 0, size = 4,
           colour = COL_CTL,
           label = sprintf("\u03c0 = %.1f%%", best$C$pi * 100)) +
  annotate("text", x = 3, y = best$T$pi * 100 + 2.5, hjust = 0, size = 4,
           colour = COL_DURV,
           label = sprintf("\u03c0 = %.1f%%", best$T$pi * 100)) +
  annotate("label", x = 36, y = 12, hjust = 0.5, size = 3.8,
           fill = "white", label.size = 0.4,
           label = sprintf("\u0394\u03c0 = +%.1f pp (95%% CI: +%.1f to +%.1f)",
                           best$delta * 100,
                           ci_dlt[1] * 100, ci_dlt[2] * 100)) +
  scale_colour_manual(values = c(Control = COL_CTL, Durvalumab = COL_DURV),
                      labels = c("Control (GC alone)", "Durvalumab + GC")) +
  scale_x_continuous(breaks = seq(0, 60, 12), limits = c(0, 66),
                     expand = c(0.01, 0)) +
  scale_y_continuous(limits = c(0, 100),
                     labels = function(x) paste0(x, "%"),
                     expand = c(0.01, 0)) +
  labs(x = "Time (months)", y = "Event-Free Survival (%)", colour = NULL,
       title = "Event-Free Survival with Mixture Cure Model",
       subtitle = sprintf("NIAGARA rIPD — %s distribution", best_dist)) +
  theme_pub() +
  theme(legend.position = c(0.78, 0.88))

ggsave("outputs/figures/Fig1_KM_CureModel.png", fig1,
       width = 10, height = 7, dpi = 300, bg = "white")

# -- Figure 2: Cure regression model comparison ------------------------------
creg_df <- creg_tbl
creg_df$Model_f <- factor(
  creg_df$Model,
  levels = c("cure_only", "latency_only", "both"),
  labels = c("Cure effect\nonly", "Latency effect\nonly", "Both\neffects"))

fig2a <- ggplot(creg_df, aes(x = Model_f, y = Delta_AIC,
                              colour = Delta_AIC == 0)) +
  geom_segment(aes(xend = Model_f, y = 0, yend = Delta_AIC), linewidth = 1.5) +
  geom_point(size = 5) +
  geom_text(aes(label = sprintf("AIC=%.1f\n(\u0394+%.1f)", AIC, Delta_AIC)),
            vjust = -0.5, size = 3.5, colour = "black") +
  scale_colour_manual(values = c("TRUE" = COL_DURV, "FALSE" = "gray55"),
                      guide = "none") +
  scale_y_continuous(limits = c(-1, max(creg_df$Delta_AIC) + 4),
                     expand = c(0, 0)) +
  labs(x = NULL, y = "\u0394AIC from best model",
       title = "Cure regression model comparison") +
  theme_pub()

arm_long <- rbind(
  data.frame(Model_f = creg_df$Model_f,
             pi_pct  = creg_df$pi_Control_pct, Arm = "Control"),
  data.frame(Model_f = creg_df$Model_f,
             pi_pct  = creg_df$pi_Durvalumab_pct, Arm = "Durvalumab"))
arm_long$Arm <- factor(arm_long$Arm, levels = c("Control", "Durvalumab"))

fig2b <- ggplot(arm_long, aes(x = Model_f, y = pi_pct, fill = Arm)) +
  geom_col(position = position_dodge(0.7), width = 0.6,
           colour = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.1f%%", pi_pct)),
            position = position_dodge(0.7), vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = c(Control = COL_CTL, Durvalumab = COL_DURV)) +
  scale_y_continuous(limits = c(0, max(arm_long$pi_pct) + 12), expand = c(0, 0),
                     labels = function(x) paste0(x, "%")) +
  labs(x = NULL, y = "Cure fraction (%)", fill = NULL,
       title = "Estimated cure fractions by model") +
  theme_pub() + theme(legend.position = "top")

fig2 <- gridExtra::arrangeGrob(fig2a, fig2b, ncol = 2)
ggsave("outputs/figures/Fig2_CureRegression.png", fig2,
       width = 11, height = 6, dpi = 300, bg = "white")

cat("\nOutputs saved to outputs/tables/ and outputs/figures/\n")
