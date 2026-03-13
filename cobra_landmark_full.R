# ============================================================
# COBRA: Comprehensive Early Divergence & Landmark Analysis
# X ekseni = GÜN (0-90)
# ============================================================

library(dplyr)
library(survival)
library(survminer)
library(ggplot2)
library(ggpubr)

# IPD'yi oku
allIPD <- read.csv("COBRA/cobra_allIPD.csv", stringsAsFactors = FALSE)
allIPD$arm <- factor(allIPD$arm, levels = c("Rivaroxaban", "Apixaban"))

col_rivaro <- "#7B2D8E"
col_apix   <- "#2166AC"

cat("=== COBRA Early Divergence Analysis ===\n")
cat("Total N:", nrow(allIPD), "\n")
cat("Events - Rivaroxaban:", sum(allIPD$status[allIPD$arm=="Rivaroxaban"]==1),
    ", Apixaban:", sum(allIPD$status[allIPD$arm=="Apixaban"]==1), "\n\n")

# ============================================================
# A) ABSOLUTE RISK DIFFERENCE AT KEY TIMEPOINTS
# ============================================================
fit <- survfit(Surv(time, status) ~ arm, data = allIPD)

timepoints <- c(7, 14, 21, 30, 60, 90)
surv_at <- summary(fit, times = timepoints)

# Her arm için survival extract
arms <- levels(allIPD$arm)
rd_table <- data.frame(Day = timepoints)

for (a in arms) {
  idx <- surv_at$strata == paste0("arm=", a)
  rd_table[[paste0("CumInc_", a)]] <- round((1 - surv_at$surv[idx]) * 100, 2)
}

rd_table$RiskDiff_pct <- round(rd_table$CumInc_Rivaroxaban - rd_table$CumInc_Apixaban, 2)

cat("=== A) Absolute Risk Difference (%) ===\n")
print(rd_table)

# A3) Risk difference over time (fine grid)
tgrid <- seq(1, 90, by = 0.5)
surv_grid <- summary(fit, times = tgrid)

rivaro_ci <- numeric(length(tgrid))
apix_ci   <- numeric(length(tgrid))
for (i in seq_along(tgrid)) {
  idx_r <- surv_grid$strata == "arm=Rivaroxaban"
  idx_a <- surv_grid$strata == "arm=Apixaban"
  rivaro_ci[i] <- 1 - surv_grid$surv[idx_r][i]
  apix_ci[i]   <- 1 - surv_grid$surv[idx_a][i]
}

# Fix: use summary properly
surv_r <- summary(survfit(Surv(time, status) ~ 1, data = allIPD[allIPD$arm=="Rivaroxaban",]), times = tgrid)
surv_a <- summary(survfit(Surv(time, status) ~ 1, data = allIPD[allIPD$arm=="Apixaban",]), times = tgrid)

df_diff <- data.frame(
  day = surv_r$time,
  ci_rivaro = (1 - surv_r$surv) * 100,
  ci_apix   = (1 - surv_a$surv) * 100
)
df_diff$diff <- df_diff$ci_rivaro - df_diff$ci_apix

# Risk difference plot
p_diff <- ggplot(df_diff, aes(x = day, y = diff)) +
  geom_line(linewidth = 1, color = "darkred") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = c(7, 14, 21), linetype = "dotted", color = "grey50") +
  annotate("text", x = c(7, 14, 21), y = max(df_diff$diff) * 0.95,
           label = c("Day 7", "Day 14", "Day 21"), hjust = -0.1, size = 3.5) +
  scale_x_continuous(breaks = seq(0, 90, 10)) +
  labs(title = "Cumulative Risk Difference (Rivaroxaban - Apixaban)",
       x = "Days", y = "Risk Difference (%)") +
  theme_classic(base_size = 13)

# KM plot with vertical lines
fit_km <- survfit(Surv(time, status) ~ arm, data = allIPD)
p_km <- ggsurvplot(
  fit_km, data = allIPD, fun = "event",
  conf.int = FALSE, censor = FALSE,
  xlim = c(0, 90), ylim = c(0, 0.10),
  break.time.by = 10,
  palette = c(col_rivaro, col_apix),
  legend.title = "", legend.labs = c("Rivaroxaban", "Apixaban"),
  legend = c(0.25, 0.85),
  xlab = "Days", ylab = "Cumulative Incidence",
  risk.table = TRUE, risk.table.title = "No. at Risk",
  ggtheme = theme_classic(base_size = 13)
)
p_km$plot <- p_km$plot +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 0.10), breaks = seq(0, 0.10, 0.02)) +
  geom_vline(xintercept = c(7, 14, 21), linetype = "dotted", color = "grey40") +
  annotate("text", x = c(7, 14, 21), y = 0.095,
           label = c("D7", "D14", "D21"), hjust = -0.1, size = 3.5)

# ============================================================
# B) LANDMARK ANALYSES
# ============================================================
landmark_windows <- list(
  c(0, 7), c(0, 14), c(0, 21), c(21, 90), c(30, 90)
)

lm_results <- data.frame(
  Window = character(),
  N_Rivaro = integer(), N_Apix = integer(),
  Events_Rivaro = integer(), Events_Apix = integer(),
  CumRisk_Rivaro = numeric(), CumRisk_Apix = numeric(),
  HR = numeric(), HR_LCL = numeric(), HR_UCL = numeric(),
  P = numeric(),
  RiskDiff_pct = numeric(),
  stringsAsFactors = FALSE
)

for (w in landmark_windows) {
  t0 <- w[1]; t1 <- w[2]
  label <- paste0("Day ", t0, "-", t1)

  # Filter: landmark at t0, censor at t1
  d <- allIPD %>% filter(time > t0)
  d$time <- d$time - t0
  # Administrative censoring at t1-t0
  censor_at <- t1 - t0
  d$status[d$time > censor_at] <- 0
  d$time[d$time > censor_at] <- censor_at

  n_r <- sum(d$arm == "Rivaroxaban")
  n_a <- sum(d$arm == "Apixaban")
  ev_r <- sum(d$arm == "Rivaroxaban" & d$status == 1)
  ev_a <- sum(d$arm == "Apixaban" & d$status == 1)

  cr_r <- round(ev_r / n_r * 100, 2)
  cr_a <- round(ev_a / n_a * 100, 2)
  rd <- round(cr_r - cr_a, 2)

  # Cox HR
  if (ev_r + ev_a >= 2) {
    cx <- tryCatch({
      coxph(Surv(time, status) ~ arm, data = d)
    }, error = function(e) NULL)

    if (!is.null(cx)) {
      hr_val <- exp(coef(cx))[1]
      ci_val <- exp(confint(cx))[1, ]
      p_val  <- summary(cx)$coefficients[1, "Pr(>|z|)"]
    } else {
      hr_val <- NA; ci_val <- c(NA, NA); p_val <- NA
    }
  } else {
    hr_val <- NA; ci_val <- c(NA, NA); p_val <- NA
  }

  lm_results <- rbind(lm_results, data.frame(
    Window = label,
    N_Rivaro = n_r, N_Apix = n_a,
    Events_Rivaro = ev_r, Events_Apix = ev_a,
    CumRisk_Rivaro = cr_r, CumRisk_Apix = cr_a,
    HR = round(hr_val, 2), HR_LCL = round(ci_val[1], 2), HR_UCL = round(ci_val[2], 2),
    P = round(p_val, 4),
    RiskDiff_pct = rd,
    stringsAsFactors = FALSE
  ))
}

cat("\n=== B) Landmark Results ===\n")
print(lm_results)

# ============================================================
# C) PIECEWISE COX MODEL
# ============================================================
cat("\n=== C) Piecewise Cox Model ===\n")

# Create time-split dataset
cuts <- c(7, 14, 21)
allIPD$id <- seq_len(nrow(allIPD))
split_data <- survSplit(Surv(time, status) ~ ., data = allIPD, cut = cuts, episode = "period")
split_data$period <- factor(split_data$period, labels = c("0-7", "7-14", "14-21", "21-90"))

# Piecewise HR
pw_results <- data.frame(Period = character(), HR = numeric(),
                          LCL = numeric(), UCL = numeric(), P = numeric(),
                          stringsAsFactors = FALSE)

for (p_lev in levels(split_data$period)) {
  d_sub <- split_data[split_data$period == p_lev, ]
  ev_total <- sum(d_sub$status)
  if (ev_total >= 2) {
    cx <- tryCatch(coxph(Surv(tstart, time, status) ~ arm, data = d_sub), error = function(e) NULL)
    if (!is.null(cx) && !any(is.na(coef(cx)))) {
      hr_val <- exp(coef(cx))[1]
      ci_val <- exp(confint(cx))[1, ]
      p_val  <- summary(cx)$coefficients[1, "Pr(>|z|)"]
      pw_results <- rbind(pw_results, data.frame(
        Period = p_lev, HR = round(hr_val, 2),
        LCL = round(ci_val[1], 2), UCL = round(ci_val[2], 2),
        P = round(p_val, 4), stringsAsFactors = FALSE))
    } else {
      pw_results <- rbind(pw_results, data.frame(
        Period = p_lev, HR = NA, LCL = NA, UCL = NA, P = NA, stringsAsFactors = FALSE))
    }
  } else {
    pw_results <- rbind(pw_results, data.frame(
      Period = p_lev, HR = NA, LCL = NA, UCL = NA, P = NA, stringsAsFactors = FALSE))
  }
}

cat("Piecewise HRs (Apixaban vs Rivaroxaban):\n")
print(pw_results)

# Time-treatment interaction
cx_int <- coxph(Surv(tstart, time, status) ~ arm * period, data = split_data)
cat("\nTime-Treatment Interaction (ANOVA):\n")
print(anova(cx_int))

# PH test on overall model
cx_overall <- coxph(Surv(time, status) ~ arm, data = allIPD)
ph <- cox.zph(cx_overall)
cat("\nPH Assumption Test (overall):\n")
print(ph)

# ============================================================
# D) AREA UNDER SEPARATION
# ============================================================
cat("\n=== D) Area Under Separation ===\n")

# Cumulative area of risk difference
df_diff$cum_area <- cumsum(df_diff$diff * 0.5)  # 0.5 day increments
total_area <- tail(df_diff$cum_area, 1)

area_at <- sapply(c(7, 14, 21, 30), function(d) {
  idx <- which.min(abs(df_diff$day - d))
  df_diff$cum_area[idx]
})

area_pct <- round(area_at / total_area * 100, 1)
area_table <- data.frame(
  Day = c(7, 14, 21, 30),
  CumArea = round(area_at, 2),
  PctOfTotal = area_pct
)
cat("Total area (90 days):", round(total_area, 2), "\n")
print(area_table)

# ============================================================
# E) EARLY EVENT CONCENTRATION
# ============================================================
cat("\n=== E) Early Event Concentration ===\n")

event_conc <- function(arm_name) {
  ev <- allIPD %>% filter(arm == arm_name, status == 1)
  total_ev <- nrow(ev)
  sapply(c(7, 14, 21, 30), function(d) {
    n <- sum(ev$time <= d)
    c(n = n, pct = round(n / total_ev * 100, 1))
  })
}

ec_r <- event_conc("Rivaroxaban")
ec_a <- event_conc("Apixaban")

ec_table <- data.frame(
  Day = c(7, 14, 21, 30),
  Rivaro_n = ec_r["n", ], Rivaro_pct = ec_r["pct", ],
  Apix_n = ec_a["n", ], Apix_pct = ec_a["pct", ]
)
cat("Event concentration by arm:\n")
print(ec_table)

# Interval event rates
cat("\nInterval event rates:\n")
intervals <- list(c(0,7), c(7,14), c(14,21), c(21,30), c(30,60), c(60,90))
int_rates <- data.frame(Interval = character(),
                         Rivaro_events = integer(), Apix_events = integer(),
                         Rivaro_rate_per1000pd = numeric(), Apix_rate_per1000pd = numeric(),
                         stringsAsFactors = FALSE)

for (iv in intervals) {
  t0 <- iv[1]; t1 <- iv[2]; width <- t1 - t0
  ev_r <- sum(allIPD$arm == "Rivaroxaban" & allIPD$status == 1 & allIPD$time > t0 & allIPD$time <= t1)
  ev_a <- sum(allIPD$arm == "Apixaban" & allIPD$status == 1 & allIPD$time > t0 & allIPD$time <= t1)
  n_r <- sum(allIPD$arm == "Rivaroxaban" & allIPD$time > t0)
  n_a <- sum(allIPD$arm == "Apixaban" & allIPD$time > t0)
  rate_r <- round(ev_r / (n_r * width) * 1000, 2)
  rate_a <- round(ev_a / (n_a * width) * 1000, 2)
  int_rates <- rbind(int_rates, data.frame(
    Interval = paste0(t0, "-", t1),
    Rivaro_events = ev_r, Apix_events = ev_a,
    Rivaro_rate_per1000pd = rate_r, Apix_rate_per1000pd = rate_a,
    stringsAsFactors = FALSE))
}
print(int_rates)

# ============================================================
# F) SENSITIVITY — Perturbation
# ============================================================
cat("\n=== F) Sensitivity to Reconstruction Noise ===\n")

set.seed(42)
n_perturb <- 200
perturb_hrs <- numeric(n_perturb)

for (i in seq_len(n_perturb)) {
  d_pert <- allIPD
  # Add noise: ±0.5 day to event times
  noise <- runif(nrow(d_pert), -0.5, 0.5)
  d_pert$time <- pmax(d_pert$time + noise, 0.01)

  cx <- tryCatch(
    coxph(Surv(time, status) ~ arm, data = d_pert),
    error = function(e) NULL
  )
  if (!is.null(cx)) perturb_hrs[i] <- exp(coef(cx))[1]
  else perturb_hrs[i] <- NA
}

perturb_hrs <- na.omit(perturb_hrs)
cat("Perturbation HR (±0.5 day noise, 200 iterations):\n")
cat("  Median:", round(median(perturb_hrs), 3), "\n")
cat("  2.5th–97.5th:", round(quantile(perturb_hrs, 0.025), 3), "–",
    round(quantile(perturb_hrs, 0.975), 3), "\n")
cat("  All iterations HR < 1:", sum(perturb_hrs < 1) == length(perturb_hrs), "\n")

# ============================================================
# SAVE ALL PLOTS TO PDF
# ============================================================

# 1) KM with vertical lines
combined_km <- ggarrange(p_km$plot, p_km$table, ncol = 1, nrow = 2, heights = c(3, 1))
ggsave("COBRA/cobra_km_annotated.pdf", plot = combined_km, width = 9, height = 7)

# 2) Risk difference over time
ggsave("COBRA/cobra_risk_diff.pdf", plot = p_diff, width = 9, height = 5)

cat("\n=== PDF'ler kaydedildi ===\n")
cat("  COBRA/cobra_km_annotated.pdf\n")
cat("  COBRA/cobra_risk_diff.pdf\n")
