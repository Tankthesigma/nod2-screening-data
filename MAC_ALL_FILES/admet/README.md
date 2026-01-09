# NOD2-Scout Phase 5: ADMET Profiling

## Directory Structure
```
admet/
├── inputs/           # Top 100 FDA + 100 Natural products
├── outputs/          # Scored, tiered, filtered compounds
├── figures/          # 6 publication figures
├── logs/             # Failure and summary logs
├── ADMET_Report.pdf  # Complete PDF report
└── README.md
```

## Methods
- **Descriptors:** RDKit (Crippen LogP, ESOL solubility, PAINS/Brenk filters)
- **Absorption:** BOILED-Egg proxy (GI/BBB), bioavailability rules
- **Safety:** Structural alerts, CYP risk flags
- **Tiering:** Fixed cutoffs (ADMET > 0.7, Safety > 0.8 for Tier 1)

## Limitations
- ML toxicity models (hERG, AMES, hepatotox) not available offline
- Flagged as NA and marked for experimental validation

## Reproduction
```bash
python scripts/10_admet_profiling.py
```

## Citation
NOD2-Scout | ISEF 2026
