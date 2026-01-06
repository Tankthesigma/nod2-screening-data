# Phase 6 MD Simulation Parameters

> Reference parameters from wild-type NOD2 LRR simulations (Phase 6)
> These parameters are replicated exactly for Phase A1 mutant simulations.

---

## System Setup

| Parameter | Value |
|-----------|-------|
| Force Field (Protein) | AMBER14 (`amber14-all.xml`) |
| Water Model | TIP3P-FB (`amber14/tip3pfb.xml`) |
| Ligand Force Field | OpenFF-2.1.0 (GAFF fallback) |
| Box Padding | 1.0 nm |
| Ionic Strength | 0.15 M NaCl |

---

## Simulation Parameters

| Parameter | Value |
|-----------|-------|
| Temperature | 310.15 K (37C) |
| Pressure | 1.0 atm |
| Timestep | 4.0 fs (with HMR) |
| Integrator | Langevin Middle |
| Friction | 1.0 ps^-1 |
| Constraint Tolerance | 1e-6 |

---

## Hydrogen Mass Repartitioning (HMR)

| Parameter | Value |
|-----------|-------|
| Enabled | Yes |
| Hydrogen Mass | 3.0 amu |
| Constraints | HBonds |
| Rigid Water | Yes |

---

## Electrostatics

| Parameter | Value |
|-----------|-------|
| Method | PME (Particle Mesh Ewald) |
| Cutoff | 1.0 nm |
| Dispersion Correction | Enabled |

---

## Protocol

### 1. Energy Minimization
- **Steps:** 5,000

### 2. NVT Equilibration (100 ps)
Gradual heating protocol:
| Stage | Temperature | Duration |
|-------|-------------|----------|
| 1 | 50 K | 5 ps |
| 2 | 100 K | 10 ps |
| 3 | 150 K | 10 ps |
| 4 | 200 K | 15 ps |
| 5 | 250 K | 20 ps |
| 6 | 310 K | 40 ps |

### 3. NPT Equilibration
- **Duration:** 500 ps
- **Barostat:** Monte Carlo (100 step interval)

### 4. Production MD
- **Duration:** 20 ns per replicate
- **Replicates:** 3

---

## Output Settings

| Parameter | Value |
|-----------|-------|
| Report Interval | 2,500 steps (10 ps) |
| Trajectory Interval | 2,500 steps (10 ps) |
| Checkpoint Interval | 250,000 steps (1 ns) |
| Trajectory Format | DCD |

---

## Phase A1 Simulation Matrix

| Mutant | Ligand | Replicates | Total Time |
|--------|--------|------------|------------|
| R702W | Febuxostat | 3 | 60 ns |
| R702W | 20a-Dihydrocortisol | 3 | 60 ns |
| G908R | Febuxostat | 3 | 60 ns |
| G908R | 20a-Dihydrocortisol | 3 | 60 ns |
| **Total** | | **12** | **240 ns** |

---

## Expected Runtime

Based on Phase 6 performance on RTX 4060 Ti (~200-250 ns/day with HMR):
- **Per simulation (20 ns):** ~2-2.5 hours
- **Total (240 ns):** ~24-30 hours (sequential)
- **With parallelization:** Scales linearly with GPU count

---

## Notes

1. Mutant structures have backbone atoms preserved from WT; sidechains rebuilt by PDBFixer
2. Ligand poses transferred directly from Phase 6 WT docking (same coordinates)
3. Different random seeds per replicate ensure independent sampling
4. NaN energy checks at each stage for early failure detection
