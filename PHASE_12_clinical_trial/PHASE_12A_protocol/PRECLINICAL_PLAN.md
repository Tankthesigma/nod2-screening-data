# Preclinical Development Plan: CID_10120 Bufadienolide

> **DISCLAIMER:** This is a computationally designed preclinical plan for research/educational purposes and is not medical advice.

> **COMPOUND IDENTITY CORRECTION (2026-01-16):** This compound was previously mislabeled as "CID_10592 / CID_10120 Bufadienolide". The correct identity is CID_10120, a bufadienolide (cardiac glycoside). All compound information below has been corrected.

---

## Compound Information

| Property | Value |
|----------|-------|
| **Name** | 3β,5,14-Trihydroxy-5β-bufa-20,22-dienolide |
| **Common Name** | Bufadienolide |
| **Alias** | natural_top (internal project name) |
| **PubChem CID** | **10120** (previously mislabeled as 10592) |
| **SMILES** | `C[C@]12CCC3C([C@]1(CC[C@@H]2C4=COC(=O)C=C4)O)CC[C@]5([C@@]3(CC[C@@H](C5)O)C)O` |
| **Formula** | **C24H34O5** |
| **Molecular Weight** | **402.53 g/mol** |
| **Source** | Cardiac glycoside (bufadienolide class, found in toads/plants) |
| **Status** | Novel therapeutic candidate, NOT FDA-approved |
| **Note** | Known NF-κB inhibitor; NOD2 signals through NF-κB pathway |

---

## Development Goal

**Objective:** IND (Investigational New Drug) application to FDA for Phase I clinical trial

**Indication:** NOD2-associated Crohn's disease (genotype-stratified)

**Mechanism:** Direct stabilization of NOD2 LRR domain (validated by MD simulation: 100% pocket occupancy in 2/3 replicates)

---

## Preclinical Development Timeline

### Phase 1: In Vitro Studies (Months 1-6)

| Study | Objective | Duration |
|-------|-----------|----------|
| NOD2 binding assay | Confirm direct binding to LRR domain | 2 months |
| Surface plasmon resonance (SPR) | Measure binding affinity (Kd) | 2 months |
| Cellular NF-kB reporter assay | Demonstrate functional NOD2 modulation | 3 months |
| Cytokine profiling | IL-1beta, TNF-alpha, IL-6 in stimulated PBMCs | 2 months |
| Selectivity panel | Screen against related receptors (NOD1, TLRs) | 2 months |

**Go/No-Go Criteria:**
- Kd < 10 uM for NOD2 LRR binding
- >50% reduction in NF-kB activation at 10 uM
- Selectivity >10x vs NOD1

---

### Phase 2: ADMET Profiling (Months 4-9)

| Study | Objective | Duration |
|-------|-----------|----------|
| Metabolic stability | Human/mouse/rat liver microsomes | 1 month |
| CYP inhibition panel | CYP1A2, 2C9, 2C19, 2D6, 3A4 | 1 month |
| Plasma protein binding | Human, mouse, rat | 1 month |
| Permeability (Caco-2) | Intestinal absorption prediction | 1 month |
| P-gp efflux | MDR1 substrate/inhibitor assessment | 1 month |
| hERG channel | Cardiac safety (IC50) | 1 month |
| Ames test | Mutagenicity screening | 2 months |
| Micronucleus assay | Genotoxicity | 2 months |

**Go/No-Go Criteria:**
- T1/2 > 30 min in human liver microsomes
- No significant CYP3A4 inhibition (IC50 > 10 uM)
- hERG IC50 > 30 uM
- Negative Ames test

---

### Phase 3: Animal Efficacy & Safety (Months 7-18)

#### Efficacy Studies

| Model | Species | Objective | Duration |
|-------|---------|-----------|----------|
| DSS colitis | Mouse | Acute inflammation model | 3 months |
| TNBS colitis | Rat | Chronic inflammation model | 3 months |
| NOD2-knockout rescue | Mouse | Target engagement proof | 4 months |
| Pharmacokinetics | Mouse/Rat/Dog | PK parameters (Cmax, AUC, T1/2) | 3 months |

#### Safety Studies (GLP)

| Study | Species | Duration |
|-------|---------|----------|
| Single-dose toxicity | Rat, Dog | 1 month |
| 14-day repeat dose | Rat, Dog | 2 months |
| 28-day repeat dose | Rat, Dog | 3 months |
| Reproductive toxicity (Seg I) | Rat | 3 months |
| Cardiovascular safety | Dog (telemetry) | 2 months |

**Go/No-Go Criteria:**
- >30% improvement in colitis scores vs vehicle
- NOAEL identification with adequate safety margin (>10x)
- No significant cardiovascular findings

---

### Phase 4: IND-Enabling Studies (Months 16-24)

| Activity | Description | Duration |
|----------|-------------|----------|
| GLP 28-day toxicity (pivotal) | Rat and dog, full histopathology | 4 months |
| CMC development | Synthesis scale-up, stability, formulation | 6 months |
| Bioanalytical method validation | GLP-compliant LC-MS/MS assay | 2 months |
| IND writing | Module 2-5 preparation | 3 months |
| Pre-IND meeting | FDA Type B meeting | 1 month |

---

## Regulatory Pathway

```
Month 0-6:   In vitro studies
Month 4-9:   ADMET profiling
Month 7-18:  Animal efficacy + safety
Month 16-24: IND-enabling studies
Month 24-27: IND submission and FDA review
Month 28+:   Phase I clinical trial (if approved)
```

**Total timeline to Phase I:** ~30 months

---

## Budget Estimate (Rough Order of Magnitude)

| Phase | Estimated Cost |
|-------|----------------|
| In vitro studies | $200,000 - $400,000 |
| ADMET profiling | $150,000 - $300,000 |
| Animal efficacy | $300,000 - $500,000 |
| GLP toxicology | $800,000 - $1,500,000 |
| CMC/formulation | $500,000 - $1,000,000 |
| Regulatory/IND | $200,000 - $400,000 |
| **Total** | **$2.2M - $4.1M** |

*Costs vary significantly by CRO selection, study scope, and geographic location.*

---

## Risk Assessment

| Risk | Likelihood | Mitigation |
|------|------------|------------|
| Poor oral bioavailability | Medium | Prodrug strategy, formulation optimization |
| Off-target steroid effects | Medium | Selectivity profiling, dose optimization |
| Metabolic instability | Low-Medium | Structural analogs if needed |
| Toxicity findings | Low | Conservative dose selection |
| Competition from Febuxostat | High | Position as next-gen therapy |

---

## Comparison to Febuxostat Path

| Factor | Febuxostat | CID_10120 Bufadienolide |
|--------|------------|-------------------------|
| Regulatory status | FDA-approved (gout) | Novel compound |
| Path to Crohn's trial | Phase II (repurposing) | Phase I (new IND) |
| Timeline | 2-3 years | 5-7 years |
| Cost | $5-10M | $50-100M+ |
| Risk | Lower | Higher |

**Recommendation:** Prioritize Febuxostat Phase II trial while developing CID_10120 Bufadienolide as backup/next-generation candidate.

---

## References

1. FDA Guidance for Industry: IND Applications (2015)
2. ICH M3(R2): Nonclinical Safety Studies
3. ICH S7A: Safety Pharmacology Studies
4. 21 CFR Part 312: Investigational New Drug Application
