# ============================================================
# COBRA: Forensic Robustness / BEST-WORST Case Analysis
# ============================================================
library(dplyr)
library(survival)
library(survminer)
library(ggplot2)
library(ggpubr)

allIPD <- read.csv("COBRA/cobra_allIPD.csv", stringsAsFactors = FALSE)
allIPD$arm <- factor(allIPD$arm, levels = c("Rivaroxaban", "Apixaban"))

col_rivaro <- "#7B2D8E"
col_apix   <- "#2166AC"

# ============================================================
# LEVEL 1: REPORTED ANALYSIS REPLICATION
# ============================================================
cat("====== LEVEL 1: REPORTED ANALYSIS REPLICATION ======\n\n")

n_apix <- sum(allIPD$arm == "Apixaban")
n_rivaro <- sum(allIPD$arm == "Rivaroxaban")
ev_apix <- sum(allIPD$arm == "Apixaban" & allIPD$status == 1)
ev_rivaro <- sum(allIPD$arm == "Rivaroxaban" & allIPD$status == 1)

cat("Reconstructed: Apixaban", ev_apix, "/", n_apix,
    ", Rivaroxaban", ev_rivaro, "/", n_rivaro, "\n")
cat("Published:     Apixaban 44/1345, Rivaroxaban 96/1355\n\n")

# Cox HR
fit_cox <- coxph(Surv(time, status) ~ arm, data = allIPD)
hr1 <- exp(coef(fit_cox))[1]
ci1 <- exp(confint(fit_cox))[1, ]
p1 <- summary(fit_cox)$coefficients[1, "Pr(>|z|)"]

cat(sprintf("Reconstructed HR (Apix vs Rivaro): %.3f (95%% CI %.3f-%.3f), P = %.4f\n", hr1, ci1[1], ci1[2], p1))

# Risk at day 90
fit_km <- survfit(Surv(time, status) ~ arm, data = allIPD)
s90 <- summary(fit_km, times = 90)
risk_rivaro <- round((1 - s90$surv[s90$strata == "arm=Rivaroxaban"]) * 100, 2)
risk_apix <- round((1 - s90$surv[s90$strata == "arm=Apixaban"]) * 100, 2)
rd_90 <- round(risk_rivaro - risk_apix, 2)

cat(sprintf("90-day cumulative incidence: Rivaro %.2f%%, Apix %.2f%%, RD = %.2f%%\n", risk_rivaro, risk_apix, rd_90))

# RMST
rmst_diff <- function(ipd, tau = 90) {
  fit <- survfit(Surv(time, status) ~ arm, data = ipd)
  s <- summary(fit, rmean = tau)
  rmst_vals <- s$table[, "rmean"]
  # Rivaroxaban first, Apixaban second
  diff <- rmst_vals[2] - rmst_vals[1]  # Apix - Rivaro
  diff
}
rmst_d <- rmst_diff(allIPD)
cat(sprintf("RMST difference (Apix - Rivaro, tau=90): %.2f days\n\n", rmst_d))

# ============================================================
# LEVEL 2: BINARY ENDPOINT EXTREME SCENARIOS
# ============================================================
cat("====== LEVEL 2: BINARY ENDPOINT SCENARIOS ======\n\n")

# Helper: RR, RD, OR, Fisher P from 2x2
binary_analysis <- function(ev_a, n_a, ev_r, n_r, label) {
  risk_a <- ev_a / n_a
  risk_r <- ev_r / n_r
  rr <- risk_a / risk_r
  rd <- (risk_a - risk_r) * 100
  or_val <- (ev_a / (n_a - ev_a)) / (ev_r / (n_r - ev_r))

  # Fisher exact
  mat <- matrix(c(ev_a, n_a - ev_a, ev_r, n_r - ev_r), nrow = 2)
  ft <- fisher.test(mat)

  data.frame(
    Scenario = label,
    Apix_events = ev_a, Apix_N = n_a, Apix_risk = round(risk_a * 100, 2),
    Rivaro_events = ev_r, Rivaro_N = n_r, Rivaro_risk = round(risk_r * 100, 2),
    RR = round(rr, 3), OR = round(or_val, 3),
    RD_pct = round(rd, 2),
    Fisher_P = signif(ft$p.value, 4),
    stringsAsFactors = FALSE
  )
}

# Published
res_pub <- binary_analysis(44, 1345, 96, 1355, "Published (mITT)")

# Scenario A: Best-case for Apixaban
res_A <- binary_analysis(44, 1370, 131, 1390, "A: Best-case Apix")

# Scenario B: Worst-case for Apixaban
res_B <- binary_analysis(69, 1370, 96, 1390, "B: Worst-case Apix")

# Scenario C: Equal-risk missingness
imp_apix_C <- round(44 / 1345 * 25)   # ~1
imp_rivaro_C <- round(96 / 1355 * 35)  # ~2
res_C <- binary_analysis(44 + imp_apix_C, 1370, 96 + imp_rivaro_C, 1390,
                          "C: Equal-risk missing")

# Scenario D: Informative missingness (1.25x, 1.5x, 2.0x)
multipliers <- c(1.25, 1.5, 2.0)
res_D_list <- list()
for (m_a in multipliers) {
  for (m_r in multipliers) {
    obs_risk_a <- 44 / 1345
    obs_risk_r <- 96 / 1355
    imp_a <- round(obs_risk_a * m_a * 25)
    imp_r <- round(obs_risk_r * m_r * 35)
    label <- sprintf("D: Apix %.2fx, Rivaro %.2fx", m_a, m_r)
    res_D_list[[length(res_D_list) + 1]] <- binary_analysis(44 + imp_a, 1370, 96 + imp_r, 1390, label)
  }
}
res_D <- do.call(rbind, res_D_list)

# Scenario E: Tipping point
cat("=== Scenario E: Tipping Point ===\n")
# Start from published: Apix 44/1345+25=1370, Rivaro 96/1355+35=1390
# Add k events to Apix missing (0 to Rivaro missing)
tipping <- data.frame(
  Added_Apix = integer(), Added_Rivaro = integer(),
  RR = numeric(), RD = numeric(), P = numeric(),
  stringsAsFactors = FALSE
)

for (k_a in 0:25) {
  for (k_r in 0:35) {
    ea <- 44 + k_a
    er <- 96 + k_r
    rr <- (ea / 1370) / (er / 1390)
    rd <- (ea / 1370 - er / 1390) * 100
    mat <- matrix(c(ea, 1370 - ea, er, 1390 - er), nrow = 2)
    pv <- fisher.test(mat)$p.value
    tipping <- rbind(tipping, data.frame(
      Added_Apix = k_a, Added_Rivaro = k_r,
      RR = round(rr, 3), RD = round(rd, 2), P = round(pv, 4)))
  }
}

# Find tipping: RR crosses 1 with k_r = 0
tip_null <- tipping %>% filter(Added_Rivaro == 0, RR >= 1) %>% slice_min(Added_Apix)
cat("Tipping point (RR>=1, no Rivaro events added):\n")
cat("  Need", tip_null$Added_Apix[1], "extra Apix events (of 25 missing)\n")

# Find tipping: P > 0.05 with k_r = 0
tip_p05 <- tipping %>% filter(Added_Rivaro == 0, P > 0.05) %>% slice_min(Added_Apix)
cat("Tipping point (P>0.05, no Rivaro events added):\n")
cat("  Need", tip_p05$Added_Apix[1], "extra Apix events\n\n")

# Combine all Level 2
all_scenarios <- rbind(res_pub, res_A, res_B, res_C)
cat("=== Level 2 Summary ===\n")
print(all_scenarios)
cat("\nScenario D (informative missingness):\n")
print(res_D)

# ============================================================
# LEVEL 3: CLINICALLY STRUCTURED SENSITIVITY
# ============================================================
cat("\n====== LEVEL 3: CLINICALLY STRUCTURED SENSITIVITY ======\n\n")

# Group definitions
groups <- list(
  "All excluded (pooled)" = list(apix = 25, rivaro = 35),
  "Received drug, not analyzed" = list(apix = 18, rivaro = 27),
  "LTFU + Death" = list(apix = 14, rivaro = 17),
  "Excl post-rand ineligible" = list(apix = 20, rivaro = 33),  # 25-5=20, 35-2=33
  "Incl post-rand ineligible (strict ITT)" = list(apix = 25, rivaro = 35)
)

level3_results <- list()
for (gname in names(groups)) {
  g <- groups[[gname]]
  # Best-case Apix: 0 extra Apix events, all missing Rivaro have events
  best <- binary_analysis(44, 1345 + g$apix, 96 + g$rivaro, 1355 + g$rivaro,
                           paste0(gname, " — Best Apix"))
  # Worst-case Apix: all missing Apix have events, 0 extra Rivaro
  worst <- binary_analysis(44 + g$apix, 1345 + g$apix, 96, 1355 + g$rivaro,
                            paste0(gname, " — Worst Apix"))
  # Equal risk
  imp_a <- round(44 / 1345 * g$apix)
  imp_r <- round(96 / 1355 * g$rivaro)
  equal <- binary_analysis(44 + imp_a, 1345 + g$apix, 96 + imp_r, 1355 + g$rivaro,
                            paste0(gname, " — Equal risk"))
  level3_results[[length(level3_results) + 1]] <- best
  level3_results[[length(level3_results) + 1]] <- worst
  level3_results[[length(level3_results) + 1]] <- equal
}
level3_df <- do.call(rbind, level3_results)
print(level3_df)

# ============================================================
# TIME-TO-EVENT SENSITIVITY
# ============================================================
cat("\n====== TIME-TO-EVENT SENSITIVITY ======\n\n")

# Add missing patients with different event time assumptions
tte_scenarios <- list(
  "Early event (day 5)" = 5,
  "Mid event (day 30)" = 30,
  "Late event (day 75)" = 75,
  "Censored at day 45" = NA  # censored
)

tte_results <- data.frame(
  Scenario = character(), HR = numeric(), LCL = numeric(), UCL = numeric(), P = numeric(),
  stringsAsFactors = FALSE
)

# Worst-case TTE: all 25 Apix missing have events, 0 Rivaro missing events
for (sname in names(tte_scenarios)) {
  t_val <- tte_scenarios[[sname]]
  d_new <- allIPD

  if (!is.na(t_val)) {
    # Add 25 Apix events at t_val
    apix_add <- data.frame(arm = "Apixaban", time = t_val, status = 1)
    # Add 35 Rivaro censored at t_val
    rivaro_add <- data.frame(arm = "Rivaroxaban", time = t_val, status = 0)
  } else {
    # Censored
    apix_add <- data.frame(arm = "Apixaban", time = 45, status = 0)
    rivaro_add <- data.frame(arm = "Rivaroxaban", time = 45, status = 0)
  }

  d_new <- rbind(d_new,
                 apix_add[rep(1, 25), ],
                 rivaro_add[rep(1, 35), ])
  d_new$arm <- factor(d_new$arm, levels = c("Rivaroxaban", "Apixaban"))

  cx <- coxph(Surv(time, status) ~ arm, data = d_new)
  hr_v <- exp(coef(cx))[1]
  ci_v <- exp(confint(cx))[1, ]
  p_v <- summary(cx)$coefficients[1, "Pr(>|z|)"]

  tte_results <- rbind(tte_results, data.frame(
    Scenario = paste0("Worst-case Apix: ", sname),
    HR = round(hr_v, 3), LCL = round(ci_v[1], 3), UCL = round(ci_v[2], 3),
    P = signif(p_v, 4), stringsAsFactors = FALSE
  ))
}

# Best-case TTE: all 35 Rivaro missing have events, 0 Apix missing events
for (sname in names(tte_scenarios)) {
  t_val <- tte_scenarios[[sname]]
  d_new <- allIPD

  if (!is.na(t_val)) {
    rivaro_add <- data.frame(arm = "Rivaroxaban", time = t_val, status = 1)
    apix_add <- data.frame(arm = "Apixaban", time = t_val, status = 0)
  } else {
    rivaro_add <- data.frame(arm = "Rivaroxaban", time = 45, status = 0)
    apix_add <- data.frame(arm = "Apixaban", time = 45, status = 0)
  }

  d_new <- rbind(d_new,
                 apix_add[rep(1, 25), ],
                 rivaro_add[rep(1, 35), ])
  d_new$arm <- factor(d_new$arm, levels = c("Rivaroxaban", "Apixaban"))

  cx <- coxph(Surv(time, status) ~ arm, data = d_new)
  hr_v <- exp(coef(cx))[1]
  ci_v <- exp(confint(cx))[1, ]
  p_v <- summary(cx)$coefficients[1, "Pr(>|z|)"]

  tte_results <- rbind(tte_results, data.frame(
    Scenario = paste0("Best-case Apix: ", sname),
    HR = round(hr_v, 3), LCL = round(ci_v[1], 3), UCL = round(ci_v[2], 3),
    P = signif(p_v, 4), stringsAsFactors = FALSE
  ))
}

cat("Time-to-event sensitivity (Cox HR):\n")
print(tte_results)

# ============================================================
# DEATH HANDLING
# ============================================================
cat("\n====== DEATH HANDLING ======\n\n")

# Deaths: Apix 1, Rivaro 2
# Scenario 1: death as censoring (default — already in IPD)
cat("1) Death as censoring: same as published analysis (HR =", round(hr1, 3), ")\n")

# Scenario 2: death as competing risk — composite
d_death <- allIPD
# Add deaths as events
death_apix <- data.frame(arm = "Apixaban", time = 45, status = 1)  # mid-FU
death_rivaro <- data.frame(arm = "Rivaroxaban", time = c(30, 60), status = c(1, 1))
d_death <- rbind(d_death, death_apix, death_rivaro)
d_death$arm <- factor(d_death$arm, levels = c("Rivaroxaban", "Apixaban"))
cx_death <- coxph(Surv(time, status) ~ arm, data = d_death)
hr_death <- exp(coef(cx_death))[1]
ci_death <- exp(confint(cx_death))[1, ]
p_death <- summary(cx_death)$coefficients[1, "Pr(>|z|)"]
cat(sprintf("2) Death as composite event: HR = %.3f (%.3f-%.3f), P = %.4f\n",
            hr_death, ci_death[1], ci_death[2], p_death))

# Scenario 3: death as unfavorable — assign to worst outcome
# Apix death → event at last time, Rivaro deaths → event at last time
d_unfav <- allIPD
unfav_apix <- data.frame(arm = "Apixaban", time = 85, status = 1)
unfav_rivaro <- data.frame(arm = "Rivaroxaban", time = c(85, 85), status = c(1, 1))
d_unfav <- rbind(d_unfav, unfav_apix, unfav_rivaro)
d_unfav$arm <- factor(d_unfav$arm, levels = c("Rivaroxaban", "Apixaban"))
cx_unfav <- coxph(Surv(time, status) ~ arm, data = d_unfav)
hr_unfav <- exp(coef(cx_unfav))[1]
ci_unfav <- exp(confint(cx_unfav))[1, ]
p_unfav <- summary(cx_unfav)$coefficients[1, "Pr(>|z|)"]
cat(sprintf("3) Death as unfavorable (late event): HR = %.3f (%.3f-%.3f), P = %.4f\n",
            hr_unfav, ci_unfav[1], ci_unfav[2], p_unfav))

# ============================================================
# RECONSTRUCTION UNCERTAINTY (Perturbation)
# ============================================================
cat("\n====== RECONSTRUCTION UNCERTAINTY ======\n\n")

set.seed(123)
n_boot <- 500
boot_hrs <- numeric(n_boot)

for (i in seq_len(n_boot)) {
  d_pert <- allIPD
  noise <- runif(nrow(d_pert), -0.5, 0.5)
  d_pert$time <- pmax(d_pert$time + noise, 0.01)
  cx <- tryCatch(coxph(Surv(time, status) ~ arm, data = d_pert), error = function(e) NULL)
  if (!is.null(cx)) boot_hrs[i] <- exp(coef(cx))[1] else boot_hrs[i] <- NA
}

boot_hrs <- na.omit(boot_hrs)
cat("Perturbation (+-0.5 day, 500 iter):\n")
cat("  Median HR:", round(median(boot_hrs), 3), "\n")
cat("  2.5-97.5%:", round(quantile(boot_hrs, 0.025), 3), "-",
    round(quantile(boot_hrs, 0.975), 3), "\n")
cat("  All HR<1:", all(boot_hrs < 1), "\n")

# Cross-check: worst-case + perturbation
worst_hrs <- numeric(200)
for (i in seq_len(200)) {
  d_w <- allIPD
  noise <- runif(nrow(d_w), -0.5, 0.5)
  d_w$time <- pmax(d_w$time + noise, 0.01)
  # Add worst-case
  apix_w <- data.frame(arm = "Apixaban", time = runif(25, 1, 80), status = 1)
  rivaro_w <- data.frame(arm = "Rivaroxaban", time = runif(35, 1, 80), status = 0)
  d_w <- rbind(d_w, apix_w, rivaro_w)
  d_w$arm <- factor(d_w$arm, levels = c("Rivaroxaban", "Apixaban"))
  cx <- tryCatch(coxph(Surv(time, status) ~ arm, data = d_w), error = function(e) NULL)
  if (!is.null(cx)) worst_hrs[i] <- exp(coef(cx))[1] else worst_hrs[i] <- NA
}
worst_hrs <- na.omit(worst_hrs)
cat("\nWorst-case + perturbation (200 iter):\n")
cat("  Median HR:", round(median(worst_hrs), 3), "\n")
cat("  2.5-97.5%:", round(quantile(worst_hrs, 0.025), 3), "-",
    round(quantile(worst_hrs, 0.975), 3), "\n")
cat("  % iterations HR<1:", round(sum(worst_hrs < 1) / length(worst_hrs) * 100, 1), "%\n")

# ============================================================
# OVERLAY KM PLOT: Published + Worst + Best
# ============================================================

# Worst-case IPD
d_worst <- allIPD
d_worst <- rbind(d_worst,
  data.frame(arm = "Apixaban", time = runif(25, 5, 60), status = 1),
  data.frame(arm = "Rivaroxaban", time = runif(35, 5, 60), status = 0))
d_worst$arm <- factor(d_worst$arm, levels = c("Rivaroxaban", "Apixaban"))
d_worst$scenario <- "Worst-case Apix"

# Best-case IPD
d_best <- allIPD
d_best <- rbind(d_best,
  data.frame(arm = "Apixaban", time = runif(25, 5, 60), status = 0),
  data.frame(arm = "Rivaroxaban", time = runif(35, 5, 60), status = 1))
d_best$arm <- factor(d_best$arm, levels = c("Rivaroxaban", "Apixaban"))
d_best$scenario <- "Best-case Apix"

allIPD$scenario <- "Published"

# Compute cumulative incidence for each scenario
compute_ci <- function(d, scen) {
  fit <- survfit(Surv(time, status) ~ arm, data = d)
  tgrid <- seq(0.5, 90, 0.5)
  out <- data.frame()
  for (a in c("Rivaroxaban", "Apixaban")) {
    s <- summary(survfit(Surv(time, status) ~ 1, data = d[d$arm == a, ]), times = tgrid)
    out <- rbind(out, data.frame(day = s$time, ci = (1 - s$surv) * 100, arm = a, scenario = scen))
  }
  out
}

ci_pub <- compute_ci(allIPD, "Published")
ci_worst <- compute_ci(d_worst, "Worst-case Apix")
ci_best <- compute_ci(d_best, "Best-case Apix")
ci_all <- rbind(ci_pub, ci_worst, ci_best)
ci_all$label <- paste(ci_all$arm, ci_all$scenario, sep = " — ")

p_overlay <- ggplot(ci_all, aes(x = day, y = ci, color = arm, linetype = scenario)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = c("Rivaroxaban" = col_rivaro, "Apixaban" = col_apix)) +
  scale_linetype_manual(values = c("Published" = "solid", "Worst-case Apix" = "dashed", "Best-case Apix" = "dotted")) +
  scale_x_continuous(breaks = seq(0, 90, 10)) +
  labs(title = "KM Overlay: Published vs Best/Worst-Case Scenarios",
       x = "Days", y = "Cumulative Incidence (%)",
       color = "Arm", linetype = "Scenario") +
  theme_classic(base_size = 13) +
  theme(legend.position = "bottom")

ggsave("COBRA/cobra_best_overlay.pdf", plot = p_overlay, width = 10, height = 6)
cat("\nOverlay plot saved: COBRA/cobra_best_overlay.pdf\n")

cat("\n====== ANALYSIS COMPLETE ======\n")
