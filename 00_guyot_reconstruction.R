# =============================================================================
# 00_guyot_reconstruction.R
#
# Reconstructs individual patient data (rIPD) from digitized Kaplan-Meier
# curves using the Guyot algorithm.
#
# Reference:
#   Guyot P et al. BMC Med Res Methodol. 2012;12:9.
#
# Input  (place in same folder or set path below):
#   digitized_coordinates_control.csv    -- two columns, no header:
#   digitized_coordinates_durvalumab.csv    col1 = time (months)
#                                           col2 = survival (%)
#
# Output:
#   efs_ripd_excel.csv
#   columns: id, arm, time, status (1=event / 0=censored), endpoint
#
# Numbers at risk:
#   Powles T et al. NEJM 2024;391:1773-1786, Figure A.
#   3-month intervals (t = 0, 3, 6, ..., 57 months).
#
# Key algorithmic note:
#   In the final interval where n drops to 0, censored patients are extended
#   to max(t_km) rather than the interval endpoint, reflecting that these
#   patients were followed until the administrative data cutoff, not merely
#   until the last numbers-at-risk time point.
# =============================================================================

# ---------------------------------------------------------------------------
# Paths -- edit if files are in a different location
# ---------------------------------------------------------------------------
PATH_CTRL <- "digitized_coordinates_control.csv"
PATH_DURV <- "digitized_coordinates_durvalumab.csv"
PATH_OUT  <- "efs_ripd_excel.csv"

# ---------------------------------------------------------------------------
if (!requireNamespace("survival", quietly = TRUE))
  install.packages("survival", repos = "https://cloud.r-project.org")
library(survival)


# =============================================================================
# 1.  Numbers at risk
#     Source: Powles et al. NEJM 2024, Figure A (Event-Free Survival)
#     Times (months): 0  3  6  9 12 15 18 21 24 27 30 33 36 39 42 45 48 51 54 57
# =============================================================================
NAR_TIMES <- c( 0,  3,  6,  9, 12, 15, 18, 21, 24, 27,
               30, 33, 36, 39, 42, 45, 48, 51, 54, 57)

NAR_CTRL  <- c(530, 437, 381, 343, 313, 296, 281, 264, 228, 214,
               172, 132,  94,  69,  62,  24,  18,  16,   2,   0)

NAR_DURV  <- c(533, 475, 424, 386, 356, 344, 330, 315, 282, 255,
               202, 141, 115,  86,  81,  32,  20,  20,   1,   0)

stopifnot(
  length(NAR_TIMES) == 20L,
  length(NAR_CTRL)  == 20L,
  length(NAR_DURV)  == 20L,
  NAR_CTRL[1]       == 530L,
  NAR_DURV[1]       == 533L
)


# =============================================================================
# 2.  Load digitized KM coordinates
# =============================================================================
load_coords <- function(path) {
  if (!file.exists(path))
    stop("File not found: ", path,
         "\nSet PATH_CTRL / PATH_DURV at the top of this script.")

  d <- read.csv(path, header = FALSE,
                col.names  = c("t", "pct"),
                colClasses = c("numeric", "numeric"))

  # WebPlotDigitizer can return values slightly outside [0, 100] due to
  # pixel rounding; clip to valid range.
  d$pct <- pmin(pmax(d$pct, 0), 100)
  d$S   <- d$pct / 100

  d <- d[order(d$t), ]
  d$t <- pmax(d$t, 0)

  # Ensure the curve begins at S(0) = 1 exactly.
  if (d$t[1] > 0) {
    d <- rbind(data.frame(t = 0, pct = 100, S = 1), d)
  } else {
    d$S[1] <- 1
  }

  # Drop exact duplicate time points (keep first).
  d[!duplicated(round(d$t, 3)), ]
}

cat("Loading digitized coordinates...\n")
crd_C <- load_coords(PATH_CTRL)
crd_T <- load_coords(PATH_DURV)
cat(sprintf("  Control:    %d points,  t = %.2f to %.2f months\n",
            nrow(crd_C), min(crd_C$t), max(crd_C$t)))
cat(sprintf("  Durvalumab: %d points,  t = %.2f to %.2f months\n\n",
            nrow(crd_T), min(crd_T$t), max(crd_T$t)))


# =============================================================================
# 3.  Guyot reconstruction
# =============================================================================
reconstruct <- function(crd, nar_t, nar_n, arm_label, id_prefix) {

  t_km  <- crd$t
  S_km  <- crd$S
  t_max <- max(t_km)   # end of digitized curve = proxy for data cutoff

  # Interpolate KM survival at each numbers-at-risk time point.
  # Linear interpolation approximates the KM step function (Guyot et al. 2012);
  # this introduces minor error at event-dense intervals and should be noted as
  # a limitation when reporting reconstruction uncertainty.
  S_nar <- pmin(pmax(
    approx(t_km, S_km, xout = nar_t, method = "linear", rule = 2)$y,
    0), 1)

  t_out <- numeric(0)
  s_out <- integer(0)

  for (i in seq_len(length(nar_t) - 1L)) {
    t_lo <- nar_t[i];    t_hi <- nar_t[i + 1L]
    n_i  <- nar_n[i];    n_i1 <- nar_n[i + 1L]
    S_i  <- S_nar[i];    S_i1 <- S_nar[i + 1L]

    leaving <- n_i - n_i1
    if (leaving == 0L) next

    # Expected events in this interval (Guyot 2012, eq. 1).
    d_hat <- if (S_i > 1e-9) n_i * (1 - S_i1 / S_i) else 0
    d     <- as.integer(min(round(max(d_hat, 0)), leaving))
    k     <- as.integer(leaving - d)   # censored

    # Digitized KM step-down points strictly inside (t_lo, t_hi].
    in_int <- which(t_km > t_lo & t_km <= t_hi)

    # ---- Event times --------------------------------------------------------
    if (d > 0L) {
      if (length(in_int) > 0L) {
        # Reference survival: last digitized point at or before t_lo.
        prev_i <- max(which(t_km <= t_lo))
        drops  <- pmax(0, -diff(S_km[c(prev_i, in_int)]))
        if (sum(drops) < 1e-12) drops[] <- 1   # flat segment: uniform

        ev  <- as.integer(round(d * drops / sum(drops)))
        err <- d - sum(ev)
        if (err != 0L)
          ev[which.max(drops)] <- ev[which.max(drops)] + as.integer(err)
        ev <- pmax(0L, ev)

        for (j in seq_along(in_int)) {
          if (ev[j] > 0L) {
            t_out <- c(t_out, rep(t_km[in_int[j]], ev[j]))
            s_out <- c(s_out, rep(1L, ev[j]))
          }
        }
      } else {
        # No digitized steps in interval: place events uniformly.
        t_out <- c(t_out, seq(t_lo + 0.01, t_hi, length.out = d))
        s_out <- c(s_out, rep(1L, d))
      }
    }

    # ---- Censored times -----------------------------------------------------
    # KEY FIX: for the final interval (n_i1 == 0), patients who are censored
    # were followed until the data cutoff, not merely until t_hi.
    # Extend their placement to the end of the digitized curve (t_max) as a
    # proxy for the administrative censoring date.
    if (k > 0L) {
      t_cens_hi <- if (n_i1 == 0L) t_max else t_hi

      if (k == 1L) {
        # seq(from, to, length.out=1) returns `from`, not `to`.
        # Place the single patient at t_cens_hi explicitly.
        t_out <- c(t_out, t_cens_hi)
      } else {
        t_out <- c(t_out, seq(t_lo + 0.01, t_cens_hi, length.out = k))
      }
      s_out <- c(s_out, rep(0L, k))
    }
  }

  # ---- Patients beyond the last NAR time point ----------------------------
  # n_tail = 0 for both arms in NIAGARA (all accounted for by t=57),
  # so this block is a no-op here but is kept for generality.
  n_tail <- nar_n[length(nar_n)]
  if (n_tail > 0L) {
    t_tail <- seq(nar_t[length(nar_t)] + 0.01, t_max, length.out = n_tail)
    t_out  <- c(t_out, t_tail)
    s_out  <- c(s_out, rep(0L, n_tail))
  }

  # ---- Assemble output ----------------------------------------------------
  out <- data.frame(time = pmax(t_out, 0.001), status = s_out)
  out <- out[order(out$time), ]
  rownames(out) <- NULL

  out$id       <- sprintf("%s%04d", id_prefix, seq_len(nrow(out)))
  out$arm      <- arm_label
  out$endpoint <- "EFS"
  out[, c("id", "arm", "time", "status", "endpoint")]
}


# =============================================================================
# 4.  Run reconstruction
# =============================================================================
cat("Running Guyot reconstruction...\n")
ripd_C <- reconstruct(crd_C, NAR_TIMES, NAR_CTRL, "Control",    "C")
ripd_T <- reconstruct(crd_T, NAR_TIMES, NAR_DURV, "Durvalumab", "D")

cat(sprintf("  Control:    N = %d,  events = %d,  censored = %d,  max t = %.2f mo\n",
            nrow(ripd_C), sum(ripd_C$status), sum(ripd_C$status == 0L),
            max(ripd_C$time)))
cat(sprintf("  Durvalumab: N = %d,  events = %d,  censored = %d,  max t = %.2f mo\n",
            nrow(ripd_T), sum(ripd_T$status), sum(ripd_T$status == 0L),
            max(ripd_T$time)))
cat("  (Published: Control 246 events; Durvalumab 187 events)\n\n")

df_ripd <- rbind(ripd_C, ripd_T)


# =============================================================================
# 5.  Validation against published NIAGARA results
# =============================================================================
cat("Validation:\n")

km  <- survfit(Surv(time, status) ~ arm, data = df_ripd)
cox <- coxph(Surv(time, status) ~ arm, data = df_ripd)
lrt <- survdiff(Surv(time, status) ~ arm, data = df_ripd)

hr    <- exp(coef(cox)["armDurvalumab"])
ci    <- exp(confint(cox)["armDurvalumab", ])
lrp   <- 1 - pchisq(lrt$chisq, df = 1L)

sm24  <- summary(km, times = 24)
efs24 <- setNames(sm24$surv, sub("arm=", "", as.character(sm24$strata)))
med   <- setNames(summary(km)$table[, "median"],
                  sub("arm=", "", rownames(summary(km)$table)))

fmt_med <- function(x)
  ifelse(is.na(x) | is.infinite(x), "NR", sprintf("%.1f mo", x))

cat(sprintf("  HR:              %.2f (%.2f-%.2f)   [published: 0.68 (0.56-0.82)]\n",
            hr, ci[1], ci[2]))
cat(sprintf("  Log-rank p:      %s              [published: <0.001]\n",
            ifelse(lrp < 0.001, "<0.001", sprintf("%.4f", lrp))))
cat(sprintf("  24-mo EFS Ctrl:  %.1f%%             [published: 59.8%%]\n",
            efs24["Control"] * 100))
cat(sprintf("  24-mo EFS Durv:  %.1f%%             [published: 67.8%%]\n",
            efs24["Durvalumab"] * 100))
cat(sprintf("  Median EFS Ctrl: %s            [published: 46.1 mo]\n",
            fmt_med(med["Control"])))
cat(sprintf("  Median EFS Durv: %s               [published: NR]\n\n",
            fmt_med(med["Durvalumab"])))


# =============================================================================
# 6.  Save
# =============================================================================
out_dir <- dirname(PATH_OUT)
if (nchar(out_dir) > 0 && !dir.exists(out_dir))
  dir.create(out_dir, recursive = TRUE)

write.csv(df_ripd, PATH_OUT, row.names = FALSE)
cat(sprintf("Saved: %s  (N = %d,  max t = %.2f months)\n",
            PATH_OUT, nrow(df_ripd), max(df_ripd$time)))
