# Cross-Validation Summary: NOD2 LRR Febuxostat Docking

## Overview
This analysis compared docking poses from multiple AI/ML docking tools to validate
the febuxostat binding site identified by GNINA on NOD2 LRR domain.

## Methods Compared
1. **GNINA** (CNN-based scoring) - Original docking method
2. **DiffDock-L** (Diffusion model) - AI docking from Neurosnap
3. **Chai-1** (Protein-ligand structure prediction)
4. **FlowDock** - Failed (sequence mismatch error)

## Key Findings

### GNINA Results (Reference)
- Distance to pocket center: **4.64 Å**
- Key residues contacted: **3/5** (GLU1008, ASP1011, ARG1037)
- H-bonds: **4** (THR1036-N, GLU1008-O, ASP1011-OD1, ARG1037-NE)
- CNN Score: 0.746
- Affinity: -4.14 kcal/mol

### DiffDock-L Results (LRR-only)
- **100 poses generated**
- Best pose (Rank 1) distance to GNINA pocket: **77.20 Å**
- Closest pose to GNINA pocket: **Rank 34 at 17.70 Å**
- Poses within 10Å of GNINA pocket: **0/100**
- Poses within 20Å of GNINA pocket: **2/100**

### Pocket Residues Analyzed
| Residue | GNINA Contact | DiffDock Contact |
|---------|---------------|------------------|
| GLU1008 | YES (3.10 Å)  | NO               |
| ASN1010 | NO            | NO               |
| ASP1011 | YES (3.16 Å)  | NO               |
| ARG1034 | NO            | NO               |
| ARG1037 | YES (3.41 Å)  | NO               |

## Interpretation

### DiffDock Disagreement
DiffDock consistently placed febuxostat in a **different region** of NOD2 LRR,
approximately 70-80Å from the GNINA pocket. This suggests:

1. **Different binding site preference**: DiffDock's AI model identified a
   different surface as more favorable for febuxostat binding

2. **Training data bias**: DiffDock was trained on PDBbind complexes which may
   not include NOD2-like LRR domain binding sites

3. **Pocket detection**: DiffDock may not recognize the concave pocket between
   LRR repeats as a typical drug binding site

### Scientific Implications
- The GNINA pocket (near GLU1008/ASP1011/ARG1037) is **not universally recognized**
  as a binding site by AI docking methods
- This could indicate a **novel/cryptic binding site** or a **false positive**
- MD simulations showing 99-100% pocket occupancy for febuxostat suggest the
  GNINA site is stable, supporting its validity

## Conclusion

**Cross-validation status: INCONCLUSIVE**

- DiffDock does NOT confirm the GNINA binding site
- This does not invalidate the GNINA result, but suggests the binding site
  is non-canonical and may require experimental validation
- The MD simulation stability (99% pocket occupancy) provides stronger evidence
  for the binding site than docking alone

## Files Generated
- `diffdock_lrr/` - 100 ranked SDF poses from DiffDock
- `compare_diffdock_lrr_vs_gnina.py` - Comparison analysis script
- `search_all_diffdock_poses.py` - Exhaustive pose search
- `comparison_results/diffdock_lrr_vs_gnina.txt` - Summary results

## Recommendations
1. Consider experimental validation (SPR, ITC, or X-ray)
2. The MD simulation results (pocket occupancy, H-bonds) may be more
   reliable than docking cross-validation
3. Focus on replicate MD simulations as primary validation method
