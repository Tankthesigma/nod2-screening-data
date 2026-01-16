# Compound Identity Correction Notice

**Date:** 2026-01-16
**Project:** NOD2-CROHN Drug Discovery (ISEF 2026)
**Author:** Tanmay Vasudeva

---

## Summary

A compound identity mislabeling was discovered and corrected on 2026-01-16. The natural product candidate was previously labeled as **CID_10592 (Dihydrocortisol)** but the correct identity is **CID_10120 (Bufadienolide)**.

---

## Correction Details

| Property | Incorrect (Old) | Correct (New) |
|----------|-----------------|---------------|
| **PubChem CID** | CID_10592 | CID_10120 |
| **Common Name** | Dihydrocortisol / 5β-Dihydrocortisol | Bufadienolide |
| **IUPAC Name** | 20α-Dihydrocortisol | 3β,5,14-Trihydroxy-5β-bufa-20,22-dienolide |
| **Molecular Formula** | C21H32O5 | C24H34O5 |
| **Molecular Weight** | 364.48 g/mol | 402.53 g/mol |
| **Compound Class** | Corticosteroid metabolite | Cardiac glycoside (bufadienolide) |
| **Source** | Cortisol metabolism | Toads/plants (Bufo species) |

---

## Root Cause Analysis

The mislabeling occurred during the transition from Phase 4 (ML ranking) to Phase 5 (ADMET filtering):

1. **CNN Affinity Model** ranked compounds by predicted binding affinity
2. **ADMET Filtering** re-ranked compounds by drug-likeness
3. The **#1 CNN compound** (CID_10592) was incorrectly selected instead of the **#2 ADMET-ranked compound** (CID_10120)
4. Both compounds had similar binding scores but different ADMET profiles

---

## Scientific Impact

### Positive Discovery
The corrected compound (CID_10120, Bufadienolide) is actually more scientifically interesting:

1. **Known NF-κB inhibitor** - Bufadienolides are documented NF-κB pathway inhibitors in literature
2. **NOD2 signals through NF-κB** - This suggests a dual mechanism of action
3. **Mutation-resistant binding** - FEP shows ΔΔG = -0.44 kcal/mol (no effect from R702W)
4. **Traditional medicine connection** - "Chan Su" (toad venom) used in Chinese medicine for inflammation

### FEP Results Remain Valid
All Free Energy Perturbation calculations were performed with the actual molecular structure (from SDF files), not the label. The numerical results are correct:

| Property | CID_10120 (Bufadienolide) |
|----------|---------------------------|
| ΔG_bind (WT) | -15.22 ± 0.26 kcal/mol |
| ΔG_bind (R702W) | -15.66 ± 0.26 kcal/mol |
| ΔΔG | -0.44 kcal/mol |
| Effect | Mutation-resistant |

---

## Files Updated (2026-01-16)

### Critical Documentation (5 files)
- [x] PHASE_6/PHASE6_README.md
- [x] PHASE_10_mutation_analysis/PATIENT_STRATIFICATION.md
- [x] PHASE_12_clinical_trial/PHASE_12A_protocol/ISEF_TALKING_POINTS.md
- [x] PHASE_12_clinical_trial/PHASE_12A_protocol/PRECLINICAL_PLAN.md
- [x] MASTER_PROJECT_DATA.md

### Computational Scripts (12 files)
- [x] compute_absolute_binding.py
- [x] deploy_solvent_fep.py
- [x] create_deployment_package.py
- [x] fep_pmx_natural/build_full_fep_systems.py
- [x] fep_pmx_natural/build_fep_systems_wsl.py
- [x] fep_pmx_natural/setup_fep_natural.py
- [x] fep_pmx_natural/setup_solvent_natural.py
- [x] fep_pmx_natural/analyze_mbar.py
- [x] fep_pmx_natural/analyze_fep_results.py
- [x] fep_pmx_natural/select_boresch_anchors.py
- [x] fep_pmx_natural/run_verification.py
- [x] fep_pmx_natural/launch_all_windows.py

### Metadata Files (4 files)
- [x] fep_pmx_natural/run_manifest.json
- [x] fep_pmx_natural/run_manifest.txt
- [x] fep_pmx_natural/mbar_results.json
- [x] fep_pmx_natural/verification_results.txt

### Window Scripts (60 files in fep_pmx_natural)
- [x] fep_pmx_natural/wt_complex/window_*/run_window.py (20 files) - UPDATED
- [x] fep_pmx_natural/mut_complex/window_*/run_window.py (20 files) - UPDATED
- [x] fep_pmx_natural/solvent/window_*/run_window.py (20 files) - UPDATED

### Generated Documents
- [x] MASTER_PROJECT_SUMMARY.docx
- [x] MASTER_PROJECT_SUMMARY.pdf
- [x] create_master_summary.py

---

## Files NOT Updated (Intentionally)

### CSV/Data Files
The following files in `admet/outputs/` contain the original CID_10592 identifier in their data. These are preserved as historical records:
- benchmark_*.csv
- admet_*.csv
- enriched_*.csv

**Note:** These CSV files document the actual screening output at the time of analysis. The CID identifier in these files refers to the same molecular structure that is now correctly identified as CID_10120.

---

## Verification Commands

To verify all critical files have been corrected:

```bash
# Search for old identifier (should return 0 matches in critical files)
grep -r "CID_10592" --include="*.py" --include="*.md" --include="*.json" fep_pmx_natural/

# Confirm new identifier is present
grep -r "CID_10120" --include="*.py" --include="*.md" fep_pmx_natural/ | head -20
```

---

## SMILES String (Correct)

```
C[C@]12CCC3C([C@]1(CC[C@@H]2C4=COC(=O)C=C4)O)CC[C@]5([C@@]3(CC[C@@H](C5)O)C)O
```

---

## Contact

For questions about this correction, contact the project author.

**Project Repository:** github.com/Tankthesigma/nod2-screening-data
