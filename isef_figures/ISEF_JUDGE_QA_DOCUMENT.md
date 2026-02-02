# ISEF 2026 Judge Q&A Document
## Deep Learning Guided Discovery and FEP Validation of NOD2 Binders for the R702W Crohn's Variant
**Author:** Tanmay Vasudeva | **School:** Texas Virtual Academy at Hallsville

---

## ONE-SENTENCE PROJECT SUMMARY
"I built a screening pipeline for NOD2 binders and used FEP to quantify how the R702W variant changes binding for specific ligands, with apo MD providing supporting context."

---

# THE 3 HEADLINE RESULTS

### 1. FEP: Febuxostat is Mutation-Sensitive
- **ΔΔG = +2.34 kcal/mol**
- **50x weaker binding in R702W** (derived from ΔΔG)
- **8σ significance, p<0.001**

### 2. FEP: Bufadienolide Shows No Significant Difference
- **ΔΔG = -0.44 kcal/mol**
- **Not significant (1.8σ)**

### 3. Apo MD: R702W Has Reduced Sampling vs WT
- WT % frames >5Å = **12.6%**
- R702W % frames >5Å = **0%**

**Story:** Ligand-dependent binding change + mechanistic context

---

# PART 1: VIRTUAL SCREENING / MOLECULAR DOCKING

## Q: What docking software did you use?
**A:** GNINA, a CNN-based molecular docking program that uses deep learning to score protein-ligand interactions.

## Q: Why GNINA instead of AutoDock Vina?
**A:** GNINA uses convolutional neural networks trained on experimental binding data, providing more accurate affinity predictions than traditional scoring functions. It also includes SMINA rescoring for additional validation.

## Q: What were your grid box parameters?
| Parameter | Value |
|-----------|-------|
| Grid Box Size | 30 Å × 30 Å × 30 Å (cubic) |
| Box Center | Geometric centroid of LRR domain CA atoms (residues 750-1040) |
| LRR Domain Residues | 750-1040 (291 residues) |

## Q: What search parameters did you use?
| Parameter | Value |
|-----------|-------|
| Exhaustiveness | 32 |
| Number of Modes | 9 |
| Energy Range | GNINA default |

## Q: How does CNN scoring work?
- **CNN Model:** Default GNINA ensemble
- **Score Type:** CNNaffinity (POSITIVE scale - higher = better binding)
- **Score Range:** 2.47 - 7.34 kcal/mol across library
- **Score Mean:** 4.14 kcal/mol

## Q: How did you prepare the receptor?
| Step | Method |
|------|--------|
| Input Format | PDBQT |
| Structure Source | AlphaFold (AF-Q9HC29) |
| Hydrogens | Retained in PDBQT |
| Charges | Gasteiger-Marsili |
| Total Residues | 1040 |

## Q: Did you validate your docking with other methods?
**A:** Yes, I performed cross-validation with DiffDock, Chai-1, and FlowDock. However, the results were **INCONCLUSIVE**:

| Method | COM Distance to GNINA | Key Contacts | Confidence |
|--------|----------------------|--------------|------------|
| GNINA | 0.0 Å (reference) | 3 | CNN=0.746 |
| DiffDock | 82.11 Å | 0 | 100 poses |
| Chai-1 | 53.40 Å | 0 | pTM=0.235, ipTM=0.105 (LOW) |
| FlowDock | FAILED | - | Sequence error |

**IMPORTANT:** The methods disagree, and Chai-1 had low confidence (pTM=0.235, want >0.6). **MD validation (70-80% pocket occupancy) is the primary evidence for binding site**, not docking cross-validation.

---

# PART 2: COMPOUND LIBRARY

## Q: Where did your compounds come from?
| Source | Count | Type |
|--------|-------|------|
| Natural Products (COCONUT database) | 7,414 | Various scaffolds |
| FDA-Approved Drugs (eDrug3D) | 2,152 | Repurposing candidates |
| **TOTAL** | **9,566** | Combined library |

## Q: Why both natural products and FDA drugs?
**A:** FDA-approved drugs offer faster translation to clinic (known safety profiles), while natural products provide chemical diversity that may not exist in synthetic libraries.

---

# PART 3: MACHINE LEARNING (NOD2-Scout)

## Q: What ML algorithm did you use?
**A:** XGBoost (eXtreme Gradient Boosting) classifier

## Q: What were your hyperparameters?
```python
params = {
    'n_estimators': 300,
    'max_depth': 5,
    'learning_rate': 0.05,
    'subsample': 0.8,
    'colsample_bytree': 0.8,
    'min_child_weight': 3,
    'random_state': 42
}
```

## Q: What features did you use?
| Type | Details | Count |
|------|---------|-------|
| Morgan Fingerprints | Radius 2, 2048 bits | 2,048 |
| Molecular Descriptors | logP, TPSA, HBD, HBA, aromatic rings, fraction_csp3, QED, PAINS alerts, permeability score, PgP risk | 10 |
| **TOTAL FEATURES** | - | **2,058** |

## Q: How did you label your training data?
**A:** Pseudo-labeling based on docking scores:
- **Positive class (1):** Top 20% CNN affinity per molecular weight bin
- **Negative class (0):** Bottom 20% CNN affinity per molecular weight bin
- **MW bins:** 5 bins to prevent size bias
- **Total labeled:** 856 compounds (728 train, 128 test)
- **Class balance:** 50/50

## Q: How did you validate the model?
**A:** 5-Fold Scaffold-Split Cross-Validation using Murcko scaffolds (RDKit). This ensures chemically similar compounds don't appear in both training and test sets, preventing data leakage.

## Q: What was your model performance?
| Metric | Value |
|--------|-------|
| **Scaffold CV AUC (RANGE)** | **0.85-0.93** |
| Scaffold CV Mean ± Std | 0.9177 ± 0.0326 |
| Test AUROC | 0.9005 |
| Random CV AUROC | 0.9417 |
| Logistic Regression Baseline | 0.9091 |
| **Shuffled Label Control** | **0.4776** (confirms real signal) |

## Q: What does the shuffled label control show?
**A:** Training on randomly shuffled labels produces AUC ~0.48 (random chance), confirming the model learns genuine binding patterns rather than spurious correlations.

## Q: How did you combine multiple scores?
**Composite Score Weights:**
| Component | Weight |
|-----------|--------|
| ML Score | 35% |
| CNN Affinity | 30% |
| QED (drug-likeness) | 15% |
| Permeability | 10% |
| Safety | 10% |

---

# PART 4: ADMET FILTERING

## Q: What ADMET filters did you apply?
| Filter | Cutoff |
|--------|--------|
| Molecular Weight | ≤500 Da |
| LogP | ≤5 |
| H-bond Donors | ≤5 |
| H-bond Acceptors | ≤10 |
| Rotatable Bonds | ≤10 |
| TPSA | ≤140 Ų |
| PAINS Alerts | 0 |

## Q: What is Febuxostat's ADMET profile?
| Property | Value |
|----------|-------|
| MW | 315.37 Da |
| LogP | 2.39 |
| TPSA | 86.04 Ų |
| ADMET Score | 0.926 |
| Safety Score | 0.954 |

---

# PART 5: MOLECULAR DYNAMICS

## Q: What MD software and force fields did you use?
| Parameter | Value |
|-----------|-------|
| MD Engine | OpenMM |
| Protein Force Field | AMBER14 |
| Water Model | TIP3P-FB |
| Ligand Force Field | OpenFF-2.1.0 |

## Q: How did you set up your system?
| Parameter | Value |
|-----------|-------|
| Box Padding | 1.0 nm |
| Ionic Strength | 0.15 M NaCl |
| Minimization Steps | 5,000 |
| Nonbonded Cutoff | 1.0 nm |

## Q: What was your equilibration protocol?
| Phase | Duration |
|-------|----------|
| NVT | 100 ps (heating 50K → 310K) |
| NPT | 500 ps |

## Q: What were your production parameters?
| Parameter | Value |
|-----------|-------|
| Temperature | 310.15 K |
| Pressure | 1.0 atm |
| Timestep | 4.0 fs (with hydrogen mass repartitioning) |
| Production Length | 20 ns per system |
| Replicates | 3 per condition |

## Q: What was your total simulation time?
| System | Total |
|--------|-------|
| WT Complexes | 280 ns |
| Mutant Complexes | 240 ns |
| **TOTAL** | **520 ns** |

## Q: How did you validate binding?
**Pocket Occupancy (5Å threshold):**
| Compound | Mean Occupancy |
|----------|----------------|
| Febuxostat | 70.3% |
| Bufadienolide | 79.7% |
| Decoy (negative control) | 0.0% |

## Q: What binding residues did you identify?
| Residue | Contact Frequency |
|---------|-------------------|
| GLU1008 | 50-71% |
| ARG1037 | 71-86% |
| ASP1011 | 67-93% |
| ARG1034 | 82-88% |
| LEU1007 | 61-65% |

---

# PART 6: FREE ENERGY PERTURBATION (FEP)

## Q: What is FEP and why is it the "gold standard"?
**A:** Free Energy Perturbation calculates the thermodynamic binding free energy by gradually transforming a ligand to "nothing" and measuring the free energy change. It's considered the gold standard because it:
1. Provides quantitative ΔG values (not just rankings)
2. Includes entropic contributions
3. Accounts for protein flexibility
4. Has been validated against experimental data

## Q: What FEP parameters did you use?
| Parameter | Value |
|-----------|-------|
| Software | OpenMM + openmmtools |
| Estimator | MBAR (PyMBAR) |
| Temperature | 310.0 K |
| Timestep | 2.0 fs |
| Equilibration | 100 ps per window |
| Production | 1 ns per window |

## Q: How many lambda windows did you use?
| System | Windows |
|--------|---------|
| WT Complex | 20 |
| MUT Complex | 20 |
| Solvent | 20 per compound |
| **TOTAL** | **120 windows** |

## Q: What were the FEP results for Febuxostat?
| Metric | Value |
|--------|-------|
| ΔG_bind (WT) | -10.36 ± 0.18 kcal/mol |
| ΔG_bind (R702W) | -8.02 ± 0.19 kcal/mol |
| **ΔΔG** | **+2.34 kcal/mol** |
| Fold Change | **50x weaker** |
| Significance | **8σ (p<0.001)** |

## Q: What were the FEP results for Bufadienolide?
| Metric | Value |
|--------|-------|
| ΔG_bind (WT) | -15.22 ± 0.26 kcal/mol |
| ΔG_bind (R702W) | -15.66 ± 0.26 kcal/mol |
| **ΔΔG** | **-0.44 kcal/mol** |
| Fold Change | Mutation-resistant (trend) |
| Significance | **1.8σ (Not Significant)** |

## Q: How do you calculate fold change from ΔΔG?
**A:** Using the equation: Fold change = e^(ΔΔG/RT)
- At 310K, RT ≈ 0.616 kcal/mol
- For ΔΔG = +2.34: e^(2.34/0.616) ≈ 45-50x weaker

---

# PART 7: STRUCTURE INFORMATION

## Q: What protein structure did you use?
| Property | Value |
|----------|-------|
| UniProt ID | Q9HC29 |
| Source | AlphaFold |
| Total Residues | 1040 |
| LRR Domain | 744-1040 |
| HD2 Domain | 629-743 |
| R702W Position | Residue 702 (HD2 domain) |
| Distance R702 → Binding Pocket | **79.4 Å** |

## Q: Why is the mutation distance important?
**A:** The R702W mutation is located 79.4 Å away from the binding pocket, indicating that any effect on binding must be **allosteric** (transmitted through the protein structure) rather than direct contact with the ligand.

---

# PART 8: SELECTION FUNNEL

## Q: How many compounds made it through each stage?
| Stage | Count |
|-------|-------|
| Starting Library | 9,566 |
| After ADMET | 144 |
| Tier 1 Candidates | 8 |
| MD Validated | 4 |
| FEP Validated | 2 |

---

# PART 9: CONTROLS & LIMITATIONS

## Q: What controls did you include?
| Control Type | Result | Interpretation |
|--------------|--------|----------------|
| Decoy (negative) | 0% pocket occupancy | Confirms specificity |
| Shuffled labels | AUC 0.48 | Confirms ML learns real signal |
| Budesonide (positive) | 48.7% occupancy | Known NOD2 binder |

## Q: What are the limitations of your study?
1. **AlphaFold structure uncertainty:** 17.96 Å RMSD vs experimental PDB 5IRM
2. **Docking cross-validation inconclusive:** GNINA, DiffDock, and Chai-1 predict different sites; MD validation is primary evidence
3. **Bufadienolide result not statistically significant:** 1.8σ is below the 2σ threshold
4. **No experimental validation:** All results are computational predictions
5. **ML pseudo-labels:** Model trained on docking scores, not experimental binding data
6. **FEP on only 2 ligands:** Due to computational cost, only 2 compounds were validated with FEP

---

# INTERVIEW DEFENSE - KEY QUESTIONS

## Q: "Why docking and ML if you have FEP?"
**A:** "Docking and ML prioritize candidates across 9,566 compounds efficiently. FEP is used as a high-accuracy validation step for selected ligands because it's computationally expensive (~120 windows per comparison)."

## Q: "Why only two ligands in FEP?"
**A:** "FEP is computationally expensive - each compound requires 120 lambda windows with 1 ns production per window. I used it as a validation layer to quantify mutation effects precisely for representative ligands from different chemical classes."

## Q: "Does reduced sampling mean anything biologically?"
**A:** "I'm not claiming function or therapy. I'm presenting an apo dynamics difference that is **consistent with** the ligand-dependent binding differences observed in FEP. The reduced conformational heterogeneity in R702W (0% vs 12.6% frames >5Å) provides mechanistic context."

## Q: "How do you address the AlphaFold structure limitation?"
**A:** "The AlphaFold structure differs from the experimental PDB 5IRM by 17.96 Å RMSD. However, for the LRR binding domain specifically, AlphaFold provides good confidence (pLDDT >90). The MD validation with 70-80% pocket occupancy confirms stable binding regardless of global structure differences."

## Q: "Your docking cross-validation failed. How can you trust the binding site?"
**A:** "You're right that the docking methods disagree. However, the binding site is validated by MD simulations showing 70-80% pocket occupancy over 520 ns total simulation time, with the decoy showing 0% occupancy. The MD validation is the primary evidence, not docking cross-validation."

---

# SAFE LITERATURE PHRASING

When discussing prior literature:
- **SAY:** "Some prior literature describes R702W as 'destabilizing' NOD2."
- **SAY:** "Here I evaluated conformational sampling in apo MD and found reduced sampling in R702W under this setup."
- **DON'T SAY:** "R702W is unstable" (making claims about stability)
- **DON'T SAY:** "This proves the literature is wrong" (overreaching)

---

# WORDS TO NEVER USE

- "Precision Medicine", "N=1"
- "rescue", "drug rescue"
- "stabilize", "stabilization"
- "therapeutic", "treatment"
- "activation", "unlocking"
- Any claims about therapeutic effect

**Instead use:** "binder discovery," "binding validation," "computational screening"

---

# QUICK REFERENCE - VERIFIED KEY NUMBERS

| Metric | Value |
|--------|-------|
| Compounds screened | **9,566** |
| ML AUC (scaffold CV) | **0.85-0.93** |
| Total MD time | **520 ns** |
| FEP windows | **120** |
| Febuxostat ΔΔG | **+2.34 kcal/mol (50x weaker, 8σ, p<0.001)** |
| Bufadienolide ΔΔG | **-0.44 kcal/mol (1.8σ, NS)** |
| R702 distance to pocket | **79.4 Å** |
| WT % frames >5Å | **12.6%** |
| R702W % frames >5Å | **0%** |
| WT apo RMSD | **4.05 Å** |
| R702W apo RMSD | **3.27 Å** |
| Febuxostat pocket occupancy | **70%** |
| Bufadienolide pocket occupancy | **80%** |
| Decoy pocket occupancy | **0%** |

---

*Document generated for ISEF 2026 presentation preparation*
*All values verified from codebase at C:\Users\vasud\nod2-screening-data*
