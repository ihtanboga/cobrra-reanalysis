# ============================================================
# COBRA: Dose-Phase Cox + Time-Varying HR
# ============================================================

library(dplyr)
library(survival)
library(splines)
library(ggplot2)
library(ggpubr)

allIPD <- read.csv("COBRA/cobra_allIPD.csv", stringsAsFactors = FALSE)
allIPD$arm <- factor(allIPD$arm, levels = c("Rivaroxaban", "Apixaban"))
allIPD$id <- seq_len(nrow(allIPD))
allIPD$tx <- as.numeric(allIPD$arm == "Apixaban")

col_rivaro <- "#7B2D8E"
col_apix   <- "#2166AC"

# ============================================================
# PART 1: PIECEWISE COX — 4 Dose-Phase Arms
# ============================================================
cat("====== PART 1: Dose-Phase Piecewise Cox ======\n\n")

cuts_dose <- c(7, 21)
split_dose <- survSplit(Surv(time, status) ~ ., data = allIPD, cut = cuts_dose, episode = "period")

split_dose$dose <- with(split_dose, case_when(
  arm == "Apixaban"    & period == 1 ~ "A10",
  arm == "Apixaban"    & period == 2 ~ "A5",
  arm == "Apixaban"    & period == 3 ~ "A5",
  arm == "Rivaroxaban" & period == 1 ~ "R30",
  arm == "Rivaroxaban" & period == 2 ~ "R30",
  arm == "Rivaroxaban" & period == 3 ~ "R20"
))
split_dose$dose <- factor(split_dose$dose, levels = c("R30", "R20", "A10", "A5"))

cat("Events by dose-phase:\n")
ev_table <- split_dose %>%
  group_by(dose) %>%
  summarise(
    Events = sum(status),
    Person_days = round(sum(time - tstart), 1),
    Rate_per1000pd = round(Events / sum(time - tstart) * 1000, 2),
    .groups = "drop"
  )
print(as.data.frame(ev_table))

cx_dose <- coxph(Surv(tstart, time, status) ~ dose, data = split_dose)
cat("\nHazard Ratios (vs R30):\n")
hr_df <- data.frame(
  Dose = names(coef(cx_dose)),
  HR = round(exp(coef(cx_dose)), 3),
  LCL = round(exp(confint(cx_dose))[, 1], 3),
  UCL = round(exp(confint(cx_dose))[, 2], 3),
  P = round(summary(cx_dose)$coefficients[, "Pr(>|z|)"], 4)
)
print(hr_df)

# Contrasts
coefs <- coef(cx_dose); vcm <- vcov(cx_dose)
ct <- c(-1, 0, 1)
loghr <- sum(ct * coefs)
se_lr <- sqrt(as.numeric(t(ct) %*% vcm %*% ct))
cat(sprintf("\n--- Key Pairwise ---\n"))
cat(sprintf("A10 vs R30 (loading vs loading): HR = %.3f (%.3f-%.3f), P = %.4f\n",
            hr_df$HR[2], hr_df$LCL[2], hr_df$UCL[2], hr_df$P[2]))
cat(sprintf("A5  vs R20 (maint. vs maint.):   HR = %.3f (%.3f-%.3f), P = %.4f\n",
            exp(loghr), exp(loghr - 1.96*se_lr), exp(loghr + 1.96*se_lr),
            2*pnorm(-abs(loghr/se_lr))))
cat(sprintf("R20 vs R30 (rivaro dose drop):   HR = %.3f (%.3f-%.3f), P = %.4f\n",
            hr_df$HR[1], hr_df$LCL[1], hr_df$UCL[1], hr_df$P[1]))

# ============================================================
# PART 2: Time-Varying HR via coxph(tt()) + natural splines
# ============================================================
cat("\n\n====== PART 2: Time-Varying HR (coxph tt()) ======\n\n")

# tt() function: treatment effect varies with time via natural spline
# Use 2 df (3 knots equivalent) for smooth HR(t)
fit_tt <- coxph(
  Surv(time, status) ~ tt(tx),
  data = allIPD,
  tt = function(x, t, ...) {
    # x = tx (0/1), t = event time
    # Create basis: x, x*ns(t, df=2)
    sp <- ns(t, df = 2, Boundary.knots = c(0, 90))
    cbind(apix = x, apix_t1 = x * sp[, 1], apix_t2 = x * sp[, 2])
  }
)

cat("Time-varying Cox (tt + ns(t, df=2)):\n")
print(summary(fit_tt))

# Compute HR(t) over grid
tgrid <- seq(1, 88, by = 0.5)
b <- coef(fit_tt)
sp_grid <- ns(tgrid, df = 2, Boundary.knots = c(0, 90))
loghr_grid <- b[1] + b[2] * sp_grid[, 1] + b[3] * sp_grid[, 2]

# SE via delta method
V <- vcov(fit_tt)
se_grid <- numeric(length(tgrid))
for (i in seq_along(tgrid)) {
  g <- c(1, sp_grid[i, 1], sp_grid[i, 2])
  se_grid[i] <- sqrt(as.numeric(t(g) %*% V %*% g))
}

hr_time <- data.frame(
  day = tgrid,
  HR = exp(loghr_grid),
  LCL = exp(loghr_grid - 1.96 * se_grid),
  UCL = exp(loghr_grid + 1.96 * se_grid)
)

cat("\nHR(t) at key timepoints:\n")
for (d in c(3, 7, 14, 21, 30, 45, 60, 80)) {
  row <- hr_time[which.min(abs(hr_time$day - d)), ]
  cat(sprintf("  Day %2d: HR = %.3f (%.3f-%.3f)\n", d, row$HR, row$LCL, row$UCL))
}

# ---- HR(t) Plot ----
p_hr <- ggplot(hr_time, aes(x = day, y = HR)) +
  geom_ribbon(aes(ymin = pmax(LCL, 0), ymax = pmin(UCL, 3)),
              fill = "#2166AC", alpha = 0.15) +
  geom_line(color = "#2166AC", linewidth = 1.1) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_vline(xintercept = c(7, 21), linetype = "dotted", color = "grey40") +
  annotate("text", x = 8.5, y = 1.8,
           label = "Day 7\nApix 10\u21925mg", hjust = 0, size = 3.2, color = "grey30") +
  annotate("text", x = 22.5, y = 1.8,
           label = "Day 21\nRivaro 15\u219220mg", hjust = 0, size = 3.2, color = "grey30") +
  scale_x_continuous(breaks = seq(0, 90, 10)) +
  coord_cartesian(ylim = c(0, 2.5)) +
  labs(
    title = "Time-Varying Hazard Ratio (Apixaban vs Rivaroxaban)",
    subtitle = "Cox PH with tt(treatment) + ns(time, df=2)",
    x = "Days since Randomization",
    y = "Hazard Ratio (Apixaban vs Rivaroxaban)"
  ) +
  theme_classic(base_size = 13)

ggsave("COBRA/cobra_tv_hr.pdf", plot = p_hr, width = 10, height = 6)
cat("\nSaved: COBRA/cobra_tv_hr.pdf\n")

# ---- Predicted log-hazard by arm (Predict-style) ----
# Rivaroxaban = baseline (logHR = 0), Apixaban = logHR(t) above
df_pred <- data.frame(
  day = rep(tgrid, 2),
  arm = rep(c("Rivaroxaban", "Apixaban"), each = length(tgrid)),
  loghaz = c(rep(0, length(tgrid)), loghr_grid),
  lower = c(rep(0, length(tgrid)), loghr_grid - 1.96 * se_grid),
  upper = c(rep(0, length(tgrid)), loghr_grid + 1.96 * se_grid)
)

p_pred <- ggplot(df_pred, aes(x = day, y = loghaz, color = arm, fill = arm)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.12, color = NA) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = c(7, 21), linetype = "dotted", color = "grey40") +
  annotate("text", x = 8, y = max(df_pred$upper) * 0.7, label = "D7", size = 3.5) +
  annotate("text", x = 22, y = max(df_pred$upper) * 0.7, label = "D21", size = 3.5) +
  scale_color_manual(values = c("Rivaroxaban" = col_rivaro, "Apixaban" = col_apix)) +
  scale_fill_manual(values = c("Rivaroxaban" = col_rivaro, "Apixaban" = col_apix)) +
  scale_x_continuous(breaks = seq(0, 90, 10)) +
  labs(
    title = "Predicted Log-Hazard by Treatment Over Time",
    subtitle = "Rivaroxaban = reference (0), Apixaban = log(HR(t))",
    x = "Days since Randomization",
    y = "Log Relative Hazard",
    color = "Treatment", fill = "Treatment"
  ) +
  theme_classic(base_size = 13) +
  theme(legend.position = c(0.85, 0.85))

ggsave("COBRA/cobra_predict_rms.pdf", plot = p_pred, width = 10, height = 6)
cat("Saved: COBRA/cobra_predict_rms.pdf\n")

# ---- Test: is HR time-varying? ----
# Compare tt() model (3 df) vs constant-HR model (1 df) via LR
fit_const <- coxph(Surv(time, status) ~ tx, data = allIPD)
lr_diff <- fit_tt$loglik[2] - fit_const$loglik[2]
df_diff <- length(coef(fit_tt)) - length(coef(fit_const))
p_lr <- pchisq(2 * lr_diff, df = df_diff, lower.tail = FALSE)
cat(sprintf("\nLR test (constant vs time-varying HR): chi2 = %.3f, df = %d, P = %.4f\n",
            2 * lr_diff, df_diff, p_lr))

cat("\n====== DONE ======\n")
