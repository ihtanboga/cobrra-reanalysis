# ============================================================
# COBRA: Dose-Phase Cox + rms Time-Varying Effect
# A10 (Apix 10mg, day 0-7), A5 (Apix 5mg, day 7+)
# R30 (Rivaro 15mg, day 0-21), R20 (Rivaro 20mg, day 21+)
# ============================================================

library(dplyr)
library(survival)
library(rms)
library(ggplot2)

allIPD <- read.csv("COBRA/cobra_allIPD.csv", stringsAsFactors = FALSE)
allIPD$arm <- factor(allIPD$arm, levels = c("Rivaroxaban", "Apixaban"))
allIPD$id <- seq_len(nrow(allIPD))

# ============================================================
# PART 1: PIECEWISE COX — 4 Dose-Phase Arms
# ============================================================
cat("====== PART 1: Dose-Phase Piecewise Cox ======\n\n")

# Split at day 7 and day 21
cuts <- c(7, 21)
split_data <- survSplit(Surv(time, status) ~ ., data = allIPD, cut = cuts, episode = "period")

# Assign dose-phase label
split_data$dose <- with(split_data, case_when(
  arm == "Apixaban"    & period == 1 ~ "A10",   # Apix 10mg BID, day 0-7
  arm == "Apixaban"    & period == 2 ~ "A5",    # Apix 5mg BID, day 7-21
  arm == "Apixaban"    & period == 3 ~ "A5",    # Apix 5mg BID, day 21+
  arm == "Rivaroxaban" & period == 1 ~ "R30",   # Rivaro 15mg BID, day 0-7
  arm == "Rivaroxaban" & period == 2 ~ "R30",   # Rivaro 15mg BID, day 7-21
  arm == "Rivaroxaban" & period == 3 ~ "R20"    # Rivaro 20mg QD, day 21+
))
split_data$dose <- factor(split_data$dose, levels = c("R30", "R20", "A10", "A5"))

# Event counts per dose-phase
cat("Events by dose-phase:\n")
ev_table <- split_data %>%
  group_by(dose) %>%
  summarise(
    N_intervals = n(),
    Events = sum(status),
    Person_days = sum(time - tstart),
    Rate_per1000pd = round(Events / Person_days * 1000, 2),
    .groups = "drop"
  )
print(as.data.frame(ev_table))

# Cox model: dose-phase
cx_dose <- coxph(Surv(tstart, time, status) ~ dose, data = split_data)
cat("\nCox model — Dose-phase (ref = R30):\n")
print(summary(cx_dose))

# Extract HRs
hr_table <- data.frame(
  Comparison = names(coef(cx_dose)),
  HR = round(exp(coef(cx_dose)), 3),
  LCL = round(exp(confint(cx_dose))[, 1], 3),
  UCL = round(exp(confint(cx_dose))[, 2], 3),
  P = round(summary(cx_dose)$coefficients[, "Pr(>|z|)"], 4)
)
cat("\nHazard Ratios (vs R30 reference):\n")
print(hr_table)

# Pairwise: A10 vs R30, A5 vs R20
cat("\n--- Key Pairwise Comparisons ---\n")

# A10 vs R30: already in model (doseA10 vs ref)
cat(sprintf("A10 vs R30: HR = %.3f (%.3f-%.3f), P = %.4f\n",
            hr_table$HR[hr_table$Comparison == "doseA10"],
            hr_table$LCL[hr_table$Comparison == "doseA10"],
            hr_table$UCL[hr_table$Comparison == "doseA10"],
            hr_table$P[hr_table$Comparison == "doseA10"]))

# A5 vs R20: need contrast
# log(HR_A5/R20) = log(HR_A5/R30) - log(HR_R20/R30)
coefs <- coef(cx_dose)
vcov_m <- vcov(cx_dose)
# A5 vs R20 contrast: c(0, -1, 0, 1) for (R20, A10, A5) → depends on order
# dose levels: R30(ref), R20, A10, A5
# coefs: doseR20, doseA10, doseA5
contrast_A5vR20 <- c(-1, 0, 1)  # A5 - R20
loghr <- sum(contrast_A5vR20 * coefs)
se_loghr <- sqrt(t(contrast_A5vR20) %*% vcov_m %*% contrast_A5vR20)
hr_a5r20 <- exp(loghr)
ci_a5r20 <- exp(loghr + c(-1.96, 1.96) * as.numeric(se_loghr))
p_a5r20 <- 2 * pnorm(-abs(loghr / as.numeric(se_loghr)))
cat(sprintf("A5 vs R20:  HR = %.3f (%.3f-%.3f), P = %.4f\n",
            hr_a5r20, ci_a5r20[1], ci_a5r20[2], p_a5r20))

# A10 vs R30 (loading vs loading)
# A5 vs R30
cat(sprintf("A5 vs R30:  HR = %.3f (%.3f-%.3f), P = %.4f\n",
            hr_table$HR[hr_table$Comparison == "doseA5"],
            hr_table$LCL[hr_table$Comparison == "doseA5"],
            hr_table$UCL[hr_table$Comparison == "doseA5"],
            hr_table$P[hr_table$Comparison == "doseA5"]))

# R20 vs R30
cat(sprintf("R20 vs R30: HR = %.3f (%.3f-%.3f), P = %.4f\n",
            hr_table$HR[hr_table$Comparison == "doseR20"],
            hr_table$LCL[hr_table$Comparison == "doseR20"],
            hr_table$UCL[hr_table$Comparison == "doseR20"],
            hr_table$P[hr_table$Comparison == "doseR20"]))

# ============================================================
# PART 2: rms — Time-Varying Treatment Effect
# tx * rcs(time, 4) with Frank Harrell's cph
# ============================================================
cat("\n\n====== PART 2: rms Time-Varying Effect ======\n\n")

# rms needs tstart/time format with time-varying covariate
# Use the split dataset but create midpoint time for the interaction
split_data$tmid <- (split_data$tstart + split_data$time) / 2
split_data$tx <- ifelse(split_data$arm == "Apixaban", 1, 0)
split_data$tx <- factor(split_data$tx, levels = c(0, 1), labels = c("Rivaroxaban", "Apixaban"))

dd <- datadist(split_data)
options(datadist = "dd")

# cph with time-varying effect: tx * rcs(tmid, 4)
fit_rms <- cph(Surv(tstart, time, status) ~ tx * rcs(tmid, 4),
               data = split_data, x = TRUE, y = TRUE)

cat("rms model summary:\n")
print(fit_rms)
cat("\nANOVA:\n")
print(anova(fit_rms))

# Predict HR over time grid
tgrid <- seq(1, 85, by = 1)

# Use contrast() for HR at each time point
ct <- contrast(fit_rms,
               list(tx = "Apixaban",    tmid = tgrid),
               list(tx = "Rivaroxaban", tmid = tgrid))

hr_df <- data.frame(
  time = tgrid,
  logHR = ct$Contrast,
  se = ct$SE
)
hr_df$HR <- exp(hr_df$logHR)
hr_df$LCL <- exp(hr_df$logHR - 1.96 * hr_df$se)
hr_df$UCL <- exp(hr_df$logHR + 1.96 * hr_df$se)

# ggplot: HR(t) curve
p_hr <- ggplot(hr_df, aes(x = time, y = HR)) +
  geom_ribbon(aes(ymin = LCL, ymax = UCL), fill = "#2166AC", alpha = 0.18) +
  geom_line(color = "#2166AC", linewidth = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(7, 21), linetype = "dotted", color = "grey50") +
  annotate("text", x = 7, y = max(hr_df$UCL, na.rm = TRUE) * 0.95,
           label = "Day 7\n(Apix dose change)", hjust = -0.05, size = 3) +
  annotate("text", x = 21, y = max(hr_df$UCL, na.rm = TRUE) * 0.95,
           label = "Day 21\n(Rivaro dose change)", hjust = -0.05, size = 3) +
  scale_x_continuous(breaks = seq(0, 90, 10)) +
  scale_y_continuous(breaks = seq(0, 3, 0.5)) +
  coord_cartesian(ylim = c(0, max(hr_df$UCL, na.rm = TRUE) * 1.05)) +
  labs(
    title = "Time-Varying HR (Apixaban vs Rivaroxaban)",
    subtitle = "rms: cph(Surv ~ tx * rcs(time, 4))",
    x = "Days since Randomization",
    y = "Hazard Ratio (Apixaban vs Rivaroxaban)"
  ) +
  theme_classic(base_size = 13)

print(p_hr)

# Save PDF
ggsave("COBRA/cobra_tv_hr.pdf", plot = p_hr, width = 10, height = 6)
cat("\nPDF saved: COBRA/cobra_tv_hr.pdf\n")

# Also try Predict() approach
cat("\n--- Predict() output ---\n")
pred <- Predict(fit_rms, tmid = tgrid, tx)
p_pred <- ggplot(pred) +
  geom_vline(xintercept = c(7, 21), linetype = "dotted", color = "grey50") +
  labs(
    title = "Predicted Log-Hazard by Treatment Over Time",
    subtitle = "rms: cph(Surv ~ tx * rcs(time, 4))",
    x = "Days", y = "Log Relative Hazard"
  ) +
  theme_classic(base_size = 13)

ggsave("COBRA/cobra_predict_rms.pdf", plot = p_pred, width = 10, height = 6)
cat("PDF saved: COBRA/cobra_predict_rms.pdf\n")

cat("\n====== DONE ======\n")
