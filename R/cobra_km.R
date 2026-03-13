# ============================================================
# COBRA: Reconstructed Kaplan-Meier — Apixaban vs Rivaroxaban
# ============================================================

library(readr)
library(dplyr)
library(reconstructKM)
library(survival)
library(survminer)

# ===== 1) Digitized clicks dosyalarını oku =====
# Format: noktalı virgül ayraç, virgül ondalık (Avrupa)
rivaro_raw <- read_delim(
  "COBRA/Rivaro_Dataset.csv",
  delim = ";", col_names = c("time", "cum_inc_pct"),
  locale = locale(decimal_mark = ","),
  show_col_types = FALSE, trim_ws = TRUE
)

apix_raw <- read_delim(
  "COBRA/Apix_Dataset.csv",
  delim = ";", col_names = c("time", "cum_inc_pct"),
  locale = locale(decimal_mark = ","),
  show_col_types = FALSE, trim_ws = TRUE
)

# ===== 2) Kümülatif insidans (%) → Survival (0-1) =====
# Orijinal grafik y ekseni %0-10, digitize değerler doğrudan %
rivaro_clicks <- rivaro_raw %>%
  mutate(survival = 1 - cum_inc_pct / 100) %>%
  select(time, survival)

apix_clicks <- apix_raw %>%
  mutate(survival = 1 - cum_inc_pct / 100) %>%
  select(time, survival)

# (0, 1) satırını ekle (yoksa)
if (!any(abs(rivaro_clicks$time) < 1e-6 & abs(rivaro_clicks$survival - 1) < 1e-6)) {
  rivaro_clicks <- bind_rows(data.frame(time = 0, survival = 1), rivaro_clicks)
}
if (!any(abs(apix_clicks$time) < 1e-6 & abs(apix_clicks$survival - 1) < 1e-6)) {
  apix_clicks <- bind_rows(data.frame(time = 0, survival = 1), apix_clicks)
}

# Sırala, aynı time'da min(survival), monoton azalan garanti
collapse_corners <- function(df) {
  df <- df %>%
    group_by(time) %>%
    summarise(survival = min(survival), .groups = "drop") %>%
    arrange(time)
  df$survival <- cummin(pmin(df$survival, 1))
  df$survival[df$survival < 0] <- 0
  df
}

rivaro_clicks <- collapse_corners(rivaro_clicks)
apix_clicks   <- collapse_corners(apix_clicks)

# ===== 3) Number at Risk (NAR) tabloları =====
# Ekran görüntüsünden okunan değerler (0, 10, 20, ... 90 ay)
rivaro_NAR <- data.frame(
  time = seq(0, 90, by = 10),
  NAR  = c(1355, 1319, 1299, 1292, 1287, 1281, 1278, 1270, 1265, 1261)
)

apix_NAR <- data.frame(
  time = seq(0, 90, by = 10),
  NAR  = c(1345, 1330, 1323, 1320, 1315, 1313, 1313, 1309, 1307, 1301)
)

# Clicks'i NAR max zamanına kırp (digitize taşmaları önlemek için)
max_nar_time <- max(rivaro_NAR$time)
rivaro_clicks <- rivaro_clicks %>% filter(time <= max_nar_time)
apix_clicks   <- apix_clicks %>% filter(time <= max_nar_time)

# ===== 4) reconstructKM: format_raw_tabs =====
rivaro_aug <- format_raw_tabs(raw_NAR = rivaro_NAR, raw_surv = rivaro_clicks)
apix_aug   <- format_raw_tabs(raw_NAR = apix_NAR,   raw_surv = apix_clicks)

# Sütun adı düzeltme (gerekirse) + sıralama garantisi
fix_aug_surv <- function(aug) {
  as <- aug$aug_surv
  if (!"surv" %in% names(as) && "survival" %in% names(as))
    names(as)[names(as) == "survival"] <- "surv"
  as <- as[order(as$time), ]
  aug$aug_surv <- as
  aug
}
rivaro_aug <- fix_aug_surv(rivaro_aug)
apix_aug   <- fix_aug_surv(apix_aug)

# ===== 5) IPD Rekonstrüksiyonu =====
rivaro_recon <- KM_reconstruct(aug_NAR = rivaro_aug$aug_NAR, aug_surv = rivaro_aug$aug_surv)
apix_recon   <- KM_reconstruct(aug_NAR = apix_aug$aug_NAR,   aug_surv = apix_aug$aug_surv)

# ===== 6) Arm'ları birleştir =====
rivaro_IPD <- data.frame(
  arm    = "Rivaroxaban",
  time   = rivaro_recon$IPD_time,
  status = rivaro_recon$IPD_event
)
apix_IPD <- data.frame(
  arm    = "Apixaban",
  time   = apix_recon$IPD_time,
  status = apix_recon$IPD_event
)

allIPD <- rbind(rivaro_IPD, apix_IPD)
allIPD$arm <- factor(allIPD$arm, levels = c("Rivaroxaban", "Apixaban"))

# Kontrol
cat("Toplam IPD satır:", nrow(allIPD), "\n")
cat("  Apixaban   :", sum(allIPD$arm == "Apixaban"), "\n")
cat("  Rivaroxaban:", sum(allIPD$arm == "Rivaroxaban"), "\n")
cat("Max event time (gün): Apixaban =",
    round(max(allIPD$time[allIPD$arm == "Apixaban" & allIPD$status == 1], na.rm = TRUE), 1),
    ", Rivaroxaban =",
    round(max(allIPD$time[allIPD$arm == "Rivaroxaban" & allIPD$status == 1], na.rm = TRUE), 1), "\n")

# ===== 7) Kaplan-Meier Çizimi =====
col_apix   <- "#2166AC"   # mavi
col_rivaro <- "#7B2D8E"   # mor (orijinal grafikteki renk)

fit <- survfit(Surv(time, status) ~ arm, data = allIPD)

p <- ggsurvplot(
  fit,
  data        = allIPD,
  fun         = "event",                    # kümülatif insidans (1 - S)
  conf.int    = FALSE,
  censor      = FALSE,
  xlim        = c(0, 90),
  ylim        = c(0, 0.10),                # %0-10 → 0-0.10
  break.time.by = 10,
  palette     = c(col_rivaro, col_apix),
  legend.title = "",
  legend.labs  = c("Rivaroxaban", "Apixaban"),
  legend       = c(0.25, 0.85),
  xlab        = "Days",
  ylab        = "Cumulative Incidence",
  risk.table       = TRUE,
  risk.table.title = "No. at Risk",
  ggtheme     = theme_classic(base_size = 13)
)

# Y eksenini % olarak göster
p$plot <- p$plot +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 0.10),
    breaks = seq(0, 0.10, by = 0.02)
  )

# ===== 8) Cox Regresyonu (HR, 95% CI, p) =====
coxfit <- coxph(Surv(time, status) ~ arm, data = allIPD)
sum_cx <- summary(coxfit)
hr  <- exp(coef(coxfit))[1]
ci  <- exp(confint(coxfit))[1, ]
pvl <- sum_cx$coefficients[1, "Pr(>|z|)"]
fmt_p <- ifelse(pvl < 0.001, "<0.001", sprintf("%.3f", pvl))

cat("\n=== Cox Regression ===\n")
cat(sprintf("HR (Apixaban vs Rivaroxaban): %.2f (95%% CI %.2f–%.2f), P = %s\n",
            hr, ci[1], ci[2], fmt_p))

# Grafiğe HR anotasyonu
hr_lab <- sprintf("HR %.2f (95%% CI %.2f\u2013%.2f)\nP = %s", hr, ci[1], ci[2], fmt_p)
p$plot <- p$plot +
  annotate("text", x = 45, y = 0.09, hjust = 0.5,
           label = hr_lab, size = 4)

print(p)

# ===== 9) PH Varsayımı Testi =====
ph_test <- cox.zph(coxfit)
cat("\n=== Proportional Hazards Test ===\n")
print(ph_test)

# ===== 10) IPD'yi kaydet =====
write.csv(allIPD, "COBRA/cobra_allIPD.csv", row.names = FALSE)
cat("\nIPD kaydedildi: COBRA/cobra_allIPD.csv\n")

# ===== 11) PDF olarak kaydet =====
# ggsurvplot + risk.table → ggarrange ile tek sayfada birleştir
library(ggpubr)
combined <- ggarrange(
  p$plot, p$table,
  ncol = 1, nrow = 2,
  heights = c(3, 1)
)
ggsave("COBRA/cobra_km_plot.pdf", plot = combined, width = 9, height = 7)
cat("PDF kaydedildi: COBRA/cobra_km_plot.pdf\n")

# ============================================================
# ===== 12) LANDMARK ANALİZİ: 21. günden sonra =====
# ============================================================
# X ekseni GÜN cinsinden (0-90 gün), landmark = 21. gün
landmark <- 21

cat("\n=== Landmark Analizi (Day", landmark, ") ===\n")

# 21. günde hâlâ risk altında olanları tut
lm_IPD <- allIPD %>% filter(time > landmark)

# Time zero = Day 21 → zamanı kaydır
lm_IPD$time <- lm_IPD$time - landmark

cat("Landmark IPD satır:", nrow(lm_IPD), "\n")
cat("  Apixaban   :", sum(lm_IPD$arm == "Apixaban"), "\n")
cat("  Rivaroxaban:", sum(lm_IPD$arm == "Rivaroxaban"), "\n")
cat("  Çıkarılan (ilk 21 gün):", nrow(allIPD) - nrow(lm_IPD), "\n")

# KM fit
fit_lm <- survfit(Surv(time, status) ~ arm, data = lm_IPD)

# Çizim
p_lm <- ggsurvplot(
  fit_lm,
  data        = lm_IPD,
  fun         = "event",
  conf.int    = FALSE,
  censor      = FALSE,
  xlim        = c(0, 70),
  ylim        = c(0, 0.06),
  break.time.by = 10,
  palette     = c(col_rivaro, col_apix),
  legend.title = "",
  legend.labs  = c("Rivaroxaban", "Apixaban"),
  legend       = c(0.25, 0.85),
  xlab        = "Days after Landmark (Day 21)",
  ylab        = "Cumulative Incidence",
  risk.table       = TRUE,
  risk.table.title = "No. at Risk",
  ggtheme     = theme_classic(base_size = 13)
)

p_lm$plot <- p_lm$plot +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 0.06),
    breaks = seq(0, 0.06, by = 0.01)
  )

# Cox regresyon (landmark)
coxfit_lm <- coxph(Surv(time, status) ~ arm, data = lm_IPD)
sum_lm    <- summary(coxfit_lm)
hr_lm     <- exp(coef(coxfit_lm))[1]
ci_lm     <- exp(confint(coxfit_lm))[1, ]
pvl_lm    <- sum_lm$coefficients[1, "Pr(>|z|)"]
fmt_p_lm  <- ifelse(pvl_lm < 0.001, "<0.001", sprintf("%.3f", pvl_lm))

cat(sprintf("\nHR (Apixaban vs Rivaroxaban, landmark Day 21): %.2f (95%% CI %.2f–%.2f), P = %s\n",
            hr_lm, ci_lm[1], ci_lm[2], fmt_p_lm))

hr_lab_lm <- sprintf("HR %.2f (95%% CI %.2f\u2013%.2f)\nP = %s", hr_lm, ci_lm[1], ci_lm[2], fmt_p_lm)
p_lm$plot <- p_lm$plot +
  annotate("text", x = 35, y = 0.055, hjust = 0.5,
           label = hr_lab_lm, size = 4)

print(p_lm)

# PH testi (landmark)
ph_lm <- cox.zph(coxfit_lm)
cat("\n=== PH Test (Landmark) ===\n")
print(ph_lm)

# PDF kaydet
combined_lm <- ggarrange(
  p_lm$plot, p_lm$table,
  ncol = 1, nrow = 2,
  heights = c(3, 1)
)
ggsave("COBRA/cobra_landmark21_km.pdf", plot = combined_lm, width = 9, height = 7)
cat("Landmark PDF kaydedildi: COBRA/cobra_landmark21_km.pdf\n")
