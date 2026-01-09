# PHASE A1 ANALYSIS REPORT
## NOD2 Mutant MD Simulations
**Generated:** 2026-01-09 15:16

---
## Executive Summary

### Pocket Occupancy Results

| Mutant | Ligand | Avg Occupancy | Verdict |
|--------|--------|---------------|---------|
| G908R | febuxostat | 81.9% | ~ PARTIAL BINDING |
| G908R | natural | 100.0% | BINDING PRESERVED |
| R702W | febuxostat | 100.0% | BINDING PRESERVED |
| R702W | natural | 100.0% | BINDING PRESERVED |

---
## Methods

### Simulation Parameters
- **Mutations:** R702W, G908R (Crohn's disease-associated)
- **Ligands:** Febuxostat (drug), Natural (MDP control)
- **Replicates:** 3 per condition
- **Duration:** 20 ns each (240 ns total)
- **Force field:** AMBER14 (protein) + OpenFF-2.1.0 (ligand)
- **Temperature:** 310.15 K

### Analysis Definitions
- **Pocket residues:** GLU1008, ASP1011, ARG1037
- **Pocket occupancy cutoff:** 5.0 Å
- **Contact cutoff:** 4.0 Å
- **H-bond criteria:** 3.5 A, >=150 deg

---
## Detailed Results

### Table 1: Pocket Occupancy

| ID | Mutant | Ligand | Rep | Occupancy | Status |
|----|--------|--------|-----|-----------|--------|
| 1 | G908R | febuxostat | 1 | 45.8% | POOR |
| 2 | G908R | febuxostat | 2 | 100.0% | EXCELLENT |
| 3 | G908R | febuxostat | 3 | 100.0% | EXCELLENT |
| 4 | G908R | natural | 1 | 100.0% | EXCELLENT |
| 5 | G908R | natural | 2 | 100.0% | EXCELLENT |
| 6 | G908R | natural | 3 | 100.0% | EXCELLENT |
| 7 | R702W | febuxostat | 1 | 100.0% | EXCELLENT |
| 8 | R702W | febuxostat | 2 | 100.0% | EXCELLENT |
| 9 | R702W | febuxostat | 3 | 100.0% | EXCELLENT |
| 10 | R702W | natural | 1 | 100.0% | EXCELLENT |
| 11 | R702W | natural | 2 | 100.0% | EXCELLENT |
| 12 | R702W | natural | 3 | 100.0% | EXCELLENT |

### Table 2: Ligand RMSD

| ID | Mutant | Ligand | Rep | Mean (Å) | SD (Å) | Max (Å) | Status |
|----|--------|--------|-----|----------|--------|---------|--------|
| 1 | G908R | febuxostat | 1 | 181.46 | 5.43 | 188.11 | UNSTABLE |
| 2 | G908R | febuxostat | 2 | 174.25 | 21.77 | 188.45 | UNSTABLE |
| 3 | G908R | febuxostat | 3 | 112.22 | 55.22 | 189.38 | UNSTABLE |
| 4 | G908R | natural | 1 | 155.99 | 29.27 | 192.22 | UNSTABLE |
| 5 | G908R | natural | 2 | 178.49 | 19.72 | 191.50 | UNSTABLE |
| 6 | G908R | natural | 3 | 131.13 | 12.72 | 186.06 | UNSTABLE |
| 7 | R702W | febuxostat | 1 | 132.27 | 13.67 | 186.58 | UNSTABLE |
| 8 | R702W | febuxostat | 2 | 134.46 | 14.53 | 197.00 | UNSTABLE |
| 9 | R702W | febuxostat | 3 | 172.67 | 23.08 | 190.45 | UNSTABLE |
| 10 | R702W | natural | 1 | 181.62 | 6.55 | 186.75 | UNSTABLE |
| 11 | R702W | natural | 2 | 179.38 | 8.28 | 185.62 | UNSTABLE |
| 12 | R702W | natural | 3 | 110.07 | 47.25 | 183.46 | UNSTABLE |

### Table 3: Key Residue Contact Frequencies

| Mutant | Ligand | Rep | GLU1008 | ASP1011 | ARG1037 |
|--------|--------|-----|---------|---------|---------|
| G908R | febuxostat | 1 | 35.6% | 1.1% | 14.5% |
| G908R | febuxostat | 2 | 100.0% | 92.8% | 100.0% |
| G908R | febuxostat | 3 | 99.7% | 100.0% | 100.0% |
| G908R | natural | 1 | 4.4% | 98.9% | 100.0% |
| G908R | natural | 2 | 58.5% | 96.8% | 99.8% |
| G908R | natural | 3 | 100.0% | 11.5% | 85.5% |
| R702W | febuxostat | 1 | 99.9% | 100.0% | 100.0% |
| R702W | febuxostat | 2 | 42.5% | 92.0% | 100.0% |
| R702W | febuxostat | 3 | 100.0% | 100.0% | 100.0% |
| R702W | natural | 1 | 100.0% | 64.5% | 98.2% |
| R702W | natural | 2 | 100.0% | 68.8% | 100.0% |
| R702W | natural | 3 | 100.0% | 34.0% | 100.0% |

### Table 4: Binding Site RMSD

| ID | Mutant | Ligand | Rep | Mean (Å) | Status |
|----|--------|--------|-----|----------|--------|
| 1 | G908R | febuxostat | 1 | 54.81 | UNSTABLE |
| 2 | G908R | febuxostat | 2 | 51.70 | UNSTABLE |
| 3 | G908R | febuxostat | 3 | 55.52 | UNSTABLE |
| 4 | G908R | natural | 1 | 63.48 | UNSTABLE |
| 5 | G908R | natural | 2 | 50.74 | UNSTABLE |
| 6 | G908R | natural | 3 | 62.31 | UNSTABLE |
| 7 | R702W | febuxostat | 1 | 64.47 | UNSTABLE |
| 8 | R702W | febuxostat | 2 | 66.36 | UNSTABLE |
| 9 | R702W | febuxostat | 3 | 50.90 | UNSTABLE |
| 10 | R702W | natural | 1 | 48.65 | UNSTABLE |
| 11 | R702W | natural | 2 | 51.21 | UNSTABLE |
| 12 | R702W | natural | 3 | 64.32 | UNSTABLE |

### Table 5: Hydrogen Bond Summary

| ID | Mutant | Ligand | Rep | Total Events | Unique Pairs | Persistent (>25%) | Status |
|----|--------|--------|-----|--------------|--------------|-------------------|--------|
| 1 | G908R | febuxostat | 1 | 0 | 0 | 0 | OK |
| 2 | G908R | febuxostat | 2 | 0 | 0 | 0 | OK |
| 3 | G908R | febuxostat | 3 | 0 | 0 | 0 | OK |
| 4 | G908R | natural | 1 | 0 | 0 | 0 | OK |
| 5 | G908R | natural | 2 | 0 | 0 | 0 | OK |
| 6 | G908R | natural | 3 | 0 | 0 | 0 | OK |
| 7 | R702W | febuxostat | 1 | 0 | 0 | 0 | OK |
| 8 | R702W | febuxostat | 2 | 0 | 0 | 0 | OK |
| 9 | R702W | febuxostat | 3 | 0 | 0 | 0 | OK |
| 10 | R702W | natural | 1 | 0 | 0 | 0 | OK |
| 11 | R702W | natural | 2 | 0 | 0 | 0 | OK |
| 12 | R702W | natural | 3 | 0 | 0 | 0 | OK |

### Table 6: Ligand COM Distance from Pocket

| ID | Mutant | Ligand | Rep | Mean (Å) | SD (Å) | Drift (Å) | Status |
|----|--------|--------|-----|----------|--------|-----------|--------|
| 1 | G908R | febuxostat | 1 | 31.20 | 7.19 | +11.88 | SIGNIFICANT_DRIFT |
| 2 | G908R | febuxostat | 2 | 39.39 | 15.53 | +35.36 | SIGNIFICANT_DRIFT |
| 3 | G908R | febuxostat | 3 | 21.47 | 13.76 | +19.21 | SIGNIFICANT_DRIFT |
| 4 | G908R | natural | 1 | 35.36 | 8.76 | +16.13 | SIGNIFICANT_DRIFT |
| 5 | G908R | natural | 2 | 29.05 | 10.47 | +8.82 | SIGNIFICANT_DRIFT |
| 6 | G908R | natural | 3 | 34.92 | 11.84 | +28.73 | SIGNIFICANT_DRIFT |
| 7 | R702W | febuxostat | 1 | 30.12 | 5.45 | +7.83 | SIGNIFICANT_DRIFT |
| 8 | R702W | febuxostat | 2 | 24.35 | 6.22 | +11.93 | SIGNIFICANT_DRIFT |
| 9 | R702W | febuxostat | 3 | 31.25 | 10.62 | +26.18 | SIGNIFICANT_DRIFT |
| 10 | R702W | natural | 1 | 29.07 | 10.01 | +22.15 | SIGNIFICANT_DRIFT |
| 11 | R702W | natural | 2 | 25.76 | 7.05 | +14.92 | SIGNIFICANT_DRIFT |
| 12 | R702W | natural | 3 | 28.80 | 10.64 | +8.25 | SIGNIFICANT_DRIFT |

### Table 7: Temperature Stability

| ID | Mutant | Ligand | Rep | Temp (K) | Status |
|----|--------|--------|-----|----------|--------|
| 1 | G908R | febuxostat | 1 | 311.1±0.6 | STABLE |
| 2 | G908R | febuxostat | 2 | 311.2±0.7 | STABLE |
| 3 | G908R | febuxostat | 3 | 311.1±0.7 | STABLE |
| 4 | G908R | natural | 1 | 311.1±0.7 | STABLE |
| 5 | G908R | natural | 2 | 311.2±0.7 | STABLE |
| 6 | G908R | natural | 3 | 311.1±0.7 | STABLE |
| 7 | R702W | febuxostat | 1 | 311.1±0.7 | STABLE |
| 8 | R702W | febuxostat | 2 | 311.2±0.7 | STABLE |
| 9 | R702W | febuxostat | 3 | 311.2±0.7 | STABLE |
| 10 | R702W | natural | 1 | 311.1±0.7 | STABLE |
| 11 | R702W | natural | 2 | 311.1±0.7 | STABLE |
| 12 | R702W | natural | 3 | 311.1±0.7 | STABLE |

---
## Conclusions

### Key Findings

1. **R702W Mutant:**
   - Febuxostat binding: 100.0% occupancy
   - Natural ligand binding: 100.0% occupancy

2. **G908R Mutant:**
   - Febuxostat binding: 81.9% occupancy
   - Natural ligand binding: 100.0% occupancy

### Precision Medicine Implications

Based on pocket occupancy results:
- **Occupancy >=95%:** Febuxostat binding preserved, patient likely responsive
- **Occupancy 80-95%:** Partial binding, monitor treatment response
- **Occupancy <80%:** Binding impaired, consider alternative therapy

---
## References

- Wild-type NOD2 baseline: 99-100% pocket occupancy (Phase 6/7)
- R702W mutation: ~10% of Crohn's patients
- G908R mutation: ~4.6% of Crohn's patients
