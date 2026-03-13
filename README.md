# COBRRA Trial — Independent Methodological Reanalysis

> **"The COBRRA Trial: A Methodological Commentary on Apixaban versus Rivaroxaban for Bleeding Risk in Acute Venous Thromboembolism"**
>
> *Ibrahim Halil Tanboga, MD, PhD*
> Nisantasi University, Medical School — Department of Cardiology and Biostatistics

---

## 📖 Overview

This repository contains all data and R analysis code supporting an **independent methodological reanalysis** of the **COBRRA trial** (Castellucci et al., *New England Journal of Medicine*, March 2026; 394:1051–60).

COBRRA is the first large-scale, head-to-head randomized controlled trial comparing **apixaban** versus **rivaroxaban** for bleeding safety in acute venous thromboembolism (VTE). The trial randomized 2,760 patients and found clinically relevant bleeding in **3.3% (apixaban) vs. 7.1% (rivaroxaban)** — a 54% relative risk reduction (RR 0.46; 95% CI 0.33–0.65; P<0.001).

This reanalysis goes beyond the published results by reconstructing pseudo-individual patient data (pseudo-IPD) from published Kaplan–Meier curves and applying a suite of advanced survival analyses to rigorously test the robustness and temporal dynamics of the reported treatment effect.

---

## 🗂️ Repository Structure

```
├── data/
│   ├── Apix_Dataset.csv        # Digitized Kaplan–Meier data — Apixaban arm
│   ├── Rivaro_Dataset.csv      # Digitized Kaplan–Meier data — Rivaroxaban arm
│   └── cobra_allIPD.csv        # Reconstructed pseudo-IPD (Guyot algorithm output)
│
├── R/
│   ├── cobra_km.R              # Kaplan–Meier reconstruction & primary Cox model
│   ├── cobra_landmark_full.R   # Landmark analyses (day 21 and day 30)
│   ├── cobra_dose.R            # Phase-specific piecewise Cox model (dose-phase coding)
│   ├── cobra_dose2.R           # Time-varying hazard ratio model (spline interaction)
│   └── cobra_best.R            # Best-case / worst-case sensitivity analysis
│
└── README.md
```

---

## 🔬 What Was Done

### 1. Kaplan–Meier Digitization & Pseudo-IPD Reconstruction

Published Kaplan–Meier survival curves from the COBRRA paper were digitized using **[WebPlotDigitizer v5.2](https://automeris.io/wpd/?v=5_2)**, an open-source tool for extracting numerical data from graphical figures. The digitized coordinate pairs, combined with the number-at-risk tables reported at days 0, 30, 60, and 90, served as input to the **Guyot algorithm** — a validated method for reconstructing approximate individual patient-level event times from aggregate KM data.

The reconstruction yielded **pseudo-IPD** for 1,345 apixaban and 1,355 rivaroxaban patients, closely matching the published 90-day cumulative incidences (3.27% vs. 6.94% reconstructed; 3.3% vs. 7.1% published). The reconstructed Cox model returned **HR 0.46 (95% CI 0.32–0.66)**, virtually identical to the published relative risk of 0.46.

### 2. Time-Varying Absolute Risk Difference Analysis

Using the reconstructed KM curves, the cumulative absolute risk difference (rivaroxaban − apixaban) was computed at each time point across the 90-day follow-up. Key finding: **~70% of the total 3.67 percentage-point difference had accumulated by day 21**, coinciding with rivaroxaban's intensified 21-day loading regimen (15 mg BID → 30 mg/day).

| Day | Risk Difference | % of Total |
|-----|-----------------|------------|
| 7   | 1.40%           | 38%        |
| 21  | 2.57%           | 70%        |
| 90  | 3.67%           | 100%       |

### 3. Landmark Analyses

To test whether the treatment effect persists *after* the differential loading period, landmark KM analyses were performed starting at **day 21** (when both drugs transition to maintenance doses) and **day 30**.

- **Day-21 landmark**: HR 0.58 (95% CI 0.34–0.98; P = 0.042) — attenuated but statistically significant
- **Day-30 landmark**: HR 0.60 (95% CI 0.35–1.04; P = 0.076) — marginally non-significant

This refutes the claim that the result is *exclusively* driven by the loading phase.

### 4. Phase-Specific Piecewise Cox Model

Treatment was re-coded as a **time-varying exposure** reflecting four clinically meaningful dose-phases:

| Phase | Regimen | Duration |
|-------|---------|----------|
| A10   | Apixaban 10 mg BID (loading) | Days 0–7 |
| A5    | Apixaban 5 mg BID (maintenance) | Days 7–90 |
| R30   | Rivaroxaban 15 mg BID (loading) | Days 0–21 |
| R20   | Rivaroxaban 20 mg QD (maintenance) | Days 21–90 |

**Key pairwise results:**

| Comparison | HR | 95% CI | P |
|------------|----|--------|---|
| A10 vs. R30 (loading vs. loading) | 0.39 | 0.20–0.76 | 0.005 |
| A5 vs. R20 (maintenance vs. maintenance) | 0.58 | 0.34–0.98 | **0.043** |

Apixaban's advantage is present in **both** the loading phase and the maintenance phase — not exclusively during the period of differential dosing.

### 5. Time-Varying Hazard Ratio Model (Spline Interaction)

A Cox model with a **treatment × natural cubic spline(time)** interaction term was fitted to characterize how the hazard ratio evolves continuously over 90 days.

| Day | HR (Apix vs. Rivar) | 95% CI |
|-----|----------------------|--------|
| 7   | 0.42                 | 0.25–0.69 |
| 21  | 0.37                 | 0.22–0.61 |
| 60  | 0.46                 | 0.24–0.86 |
| 80  | 0.74                 | 0.35–1.57 |

The likelihood ratio test for time-varying HR did not reject proportional hazards (P = 0.370), indicating the constant-HR model cannot be statistically distinguished from a time-varying one given the available power (138 total events).

### 6. Best-Case / Worst-Case Sensitivity Analysis

60 patients (25 apixaban, 35 rivaroxaban) were excluded post-randomization from the published "ITT" analysis. We systematically tested all plausible and extreme assumptions for these patients.

| Scenario | Risk Ratio | P |
|----------|------------|---|
| Published (mITT) | 0.46 | <0.001 |
| Equal-risk missing | 0.47 | <0.001 |
| Worst-case apixaban (mathematical extreme) | **0.73** | **0.044** |
| Best-case apixaban | 0.34 | <0.001 |

Even the most extreme mathematically possible scenario (all 25 missing apixaban patients have events, none of the 35 rivaroxaban patients do) yields a result that remains **statistically significant and directionally consistent**. The tipping point cannot be reached within the constraints of the available data.

---

## ✅ Key Conclusions

1. **Kaplan–Meier reconstruction** closely replicates the published COBRRA curves (HR 0.46), validating the pseudo-IPD for further analysis.

2. **~70% of the absolute risk difference** is established within the first 21 days, coinciding with rivaroxaban's more intensive loading regimen.

3. **The loading phase does not fully explain the result.** Apixaban retains a statistically significant advantage in the maintenance-phase comparison (A5 vs. R20: HR 0.58, P = 0.043) and in the day-21 landmark analysis (HR 0.58, P = 0.042).

4. **The primary finding is robust** to all tested missing-data assumptions — it cannot be reversed even under the most extreme worst-case scenario.

5. **COBRRA's result is directionally consistent** with a decade of observational, propensity-matched, and meta-analytic data suggesting a bleeding advantage for apixaban, and is pharmacologically coherent given apixaban's lower total loading dose and shorter loading duration.

---

## ⚙️ Requirements

All analyses were performed in **R** (≥ 4.2.0). Key packages used:

```r
# Survival analysis
library(survival)
library(survminer)

# KM reconstruction
library(IPDfromKM)   # or custom Guyot algorithm implementation

# Restricted cubic splines / time-varying models
library(rms)
library(splines)

# Data wrangling & visualization
library(tidyverse)
library(ggplot2)
library(patchwork)
```

---

## 📊 Data Sources

The CSV datasets in this repository were generated by:

1. **Digitizing** the published Kaplan–Meier figure from the COBRRA paper using **[WebPlotDigitizer v5.2](https://automeris.io/wpd/?v=5_2)** — an open-source, browser-based tool for extracting data from scientific figures.
2. **Reconstructing pseudo-IPD** from the digitized coordinates and number-at-risk tables using the Guyot algorithm.

> ⚠️ These are **reconstructed pseudo-IPD**, not original patient-level data. Analyses should be interpreted accordingly and cannot substitute for access to the original trial dataset.

---

## 📄 Reference

> Castellucci LA, Chen VM, Kovacs MJ, et al. Bleeding Risk with Apixaban vs. Rivaroxaban in Acute Venous Thromboembolism. *N Engl J Med* 2026;394:1051–60.

---

## 👤 Author

**Ibrahim Halil Tanboga, MD, PhD**
Department of Cardiology and Biostatistics
Nisantasi University Medical School

---

*This repository supports an independent methodological commentary and is not affiliated with the COBRRA trial investigators or study sponsors.*
