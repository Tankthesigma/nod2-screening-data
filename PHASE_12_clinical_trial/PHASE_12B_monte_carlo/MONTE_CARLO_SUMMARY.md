# Monte Carlo Clinical Trial Simulation Results

> **DISCLAIMER:** This is a computational simulation for research/educational purposes. Simulated patients are generated using assumed response rates. Results predict statistical outcomes and do not guarantee clinical efficacy.

---

## Overview

| Parameter | Value |
|-----------|-------|
| **Simulations** | 1,000 Phase II trials |
| **Patients per trial** | 210 (70 per arm) |
| **Primary comparison** | Placebo vs Febuxostat 80mg |
| **Primary analysis** | Logistic regression (adjusted) |
| **Random seed** | 42 |

---

## Assumed Response Rates

| Arm | Assumed Response Rate |
|-----|----------------------|
| Placebo | 25% |
| Febuxostat 40mg | 40% |
| Febuxostat 80mg | 50% |

*Response rates derived from historical Crohn's trial placebo rates and computational binding predictions.*

---

## Key Results

### Empirical Power

| Metric | Value |
|--------|-------|
| **Empirical Power (p < 0.05)** | **88.1%** |
| Strong Evidence (p < 0.01) | 69.4% |
| Model Convergence Rate | 100% |

### Mean Response Rates (95% CI)

| Arm | Mean | 95% CI |
|-----|------|--------|
| Placebo | 26.0% | 15.7% - 37.1% |
| Febuxostat 40mg | 41.0% | 30.0% - 52.9% |
| Febuxostat 80mg | 51.4% | 40.0% - 62.9% |

### Effect Size

| Metric | Value |
|--------|-------|
| **Mean Risk Difference** | **25.4 percentage points** |
| 95% CI | 10.0 - 41.4 pp |
| Mean Odds Ratio | 3.42 |
| Odds Ratio 95% CI | 1.52 - 7.05 |

---

## Power Curve Summary

| True Effect Size | Empirical Power |
|------------------|-----------------|
| 10% | 30.0% |
| 15% | 45.2% |
| 20% | 74.2% |
| **25%** | **89.0%** |
| 30% | 96.4% |
| 35% | 98.6% |

**Minimum Detectable Effect (MDE) at 80% power:** ~21%

---

## Interpretation

Based on 1,000 simulated trials with assumed response rates derived from our computational binding data:

- A Phase II trial of **210 patients** has approximately **88% probability** of detecting a statistically significant benefit for Febuxostat 80mg over placebo (p < 0.05)
- The trial is adequately powered (>80%) to detect our assumed 25% effect size
- Even if the true effect is smaller (20%), the trial still has 74% power

---

## ISEF Talking Point

> "Monte Carlo simulations of 1,000 virtual clinical trials predict an **88% probability** of successful Phase II outcome (p < 0.05), with a mean effect size of 25 percentage points (NNT = 4). This supports advancement of Febuxostat for genotype-selected Crohn's disease."

---

## Output Files

```
PHASE_12B_monte_carlo/
├── MONTE_CARLO_SIMULATION.py      # Main simulation code
├── MONTE_CARLO_VISUALIZATIONS.py  # Figure generation code
├── MONTE_CARLO_SUMMARY.md         # This file
└── results/
    ├── simulation_summary.txt     # Text summary
    ├── all_trials_data.csv        # Raw data (1000 trials)
    ├── patient_example.csv        # Example patient cohort
    ├── power_curve_data.csv       # Power curve data
    ├── response_rates_barchart.png
    ├── pvalue_distribution.png
    ├── effect_size_distribution.png
    ├── power_curve.png
    └── forest_plot.png
```

---

## Reproducibility

- **Random seed:** 42
- **Code:** `MONTE_CARLO_SIMULATION.py`, `MONTE_CARLO_VISUALIZATIONS.py`
- **Dependencies:** numpy, pandas, scipy, statsmodels, matplotlib

To reproduce:
```bash
python MONTE_CARLO_SIMULATION.py
python MONTE_CARLO_VISUALIZATIONS.py
```

---

## Limitations

1. **Assumed response rates:** Based on computational predictions, not clinical data
2. **Simplified patient model:** Real patients have more heterogeneity
3. **No dropout modeling:** Assumes complete follow-up after dropout adjustment
4. **Single endpoint:** CDAI-100 only; real trials have multiple endpoints

---

## References

1. FDA Guidance for Industry: Adaptive Designs for Clinical Trials (2019)
2. Lachin JM. Introduction to sample size determination. Controlled Clinical Trials (1981)
3. Julious SA. Sample sizes for clinical trials. CRC Press (2010)
