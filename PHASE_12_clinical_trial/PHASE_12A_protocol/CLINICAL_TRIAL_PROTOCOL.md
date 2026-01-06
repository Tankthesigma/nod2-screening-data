# Phase II Clinical Trial Protocol: Febuxostat in Genotype-Stratified Crohn's Disease

> **DISCLAIMER:** This is a computationally designed Phase II protocol for research/educational purposes and is not medical advice.

---

## A. Study Overview

**Title:** Phase II Randomized, Double-Blind, Placebo-Controlled Trial of Febuxostat in Genotype-Stratified Crohn's Disease

**Study Type:** Randomized, double-blind, placebo-controlled

**Duration:** 12-week treatment + 4-week follow-up (16 weeks total)

**Population:** Crohn's disease patients who are wild-type (WT) or carriers of specified NOD2 missense variants (genotype-stratified)

**Sponsor:** [To be determined]

**ClinicalTrials.gov ID:** [To be assigned]

---

## B. Patient Population

### Inclusion Criteria

1. Age 18-65 years
2. Confirmed Crohn's disease diagnosis (endoscopic + histologic confirmation)
3. Moderate-to-severe disease activity (CDAI 220-450)
4. NOD2 genotype confirmed via clinically validated assay (CLIA-certified laboratory or equivalent)
   - Variant calling must include: R702W, G908R, N852S, M863V, and L1007fs
5. Eligible genotypes: WT, R702W, G908R, N852S, or M863V
6. Stable background medications allowed:
   - 5-ASA: stable dose for >=4 weeks prior to screening
   - Corticosteroids: <=20mg prednisone equivalent, stable or tapering dose for >=2 weeks
   - Immunomodulators (azathioprine, 6-MP, methotrexate): stable dose for >=8 weeks

### Exclusion Criteria

1. **L1007fs carriers** (validated pocket residues GLU1008, ASP1011, ARG1037 are deleted by truncation; binding to validated LRR pocket not possible)
2. Prior biologic exposure requires washout per drug label and investigator judgment
   - Patients will be stratified by prior biologic exposure (yes/no)
3. Severe hepatic impairment (ALT or AST >3x upper limit of normal)
4. Renal impairment (eGFR <30 mL/min/1.73m2)
5. Pregnancy or breastfeeding
6. Active infection or intra-abdominal abscess
7. Known hypersensitivity to febuxostat or xanthine oxidase inhibitors
8. Current use of azathioprine or 6-mercaptopurine at doses requiring XO inhibitor caution

---

## C. Treatment Arms (1:1:1 Randomization)

| Arm | Treatment | Dose | Route | Duration |
|-----|-----------|------|-------|----------|
| 1 | Placebo | -- | Oral, once daily | 12 weeks |
| 2 | Febuxostat | 40 mg | Oral, once daily | 12 weeks |
| 3 | Febuxostat | 80 mg | Oral, once daily | 12 weeks |

**Note:** Febuxostat 40mg arm is exploratory. Primary comparison is Placebo vs Febuxostat 80mg.

---

## D. Stratified Randomization

Randomization will be stratified by:

1. **Genotype group:**
   - WT (no NOD2 mutation)
   - Any eligible missense mutation (R702W, G908R, N852S, or M863V)

2. **Baseline CDAI:**
   - 220-300 (moderate)
   - 301-450 (severe)

3. **Prior biologic exposure:**
   - Yes
   - No

Randomization will use permuted blocks within each stratum (block sizes 3, 6).

---

## E. Rescue Therapy Policy

**Criteria for rescue therapy initiation:**
- CDAI increases >=100 points from baseline at any scheduled visit, OR
- CDAI >450 at any scheduled visit

**Rescue therapy options:**
- Corticosteroids (prednisone or equivalent)
- Biologic therapy (per investigator judgment)

**Study continuation:**
- Patients receiving rescue therapy continue in study
- Classified as "treatment failure" for primary endpoint analysis
- All rescue therapy use documented

---

## F. Endpoints

### Primary Endpoint

**CDAI-100 Response at Week 12**
- Definition: >=100 point decrease in CDAI from baseline
- Outcome: Binary (yes/no)
- Primary comparison: Placebo vs Febuxostat 80mg
- Analysis: Logistic regression, adjusted for baseline CDAI, genotype stratum, and study site

### Secondary Endpoints (Exploratory)

*Report effect sizes with 95% confidence intervals. No p-value inflation.*

| Endpoint | Timepoint | Definition |
|----------|-----------|------------|
| Clinical remission | Week 12 | CDAI <150 |
| CRP change | Weeks 4, 8, 12 | Percent change from baseline |
| Fecal calprotectin change | Weeks 4, 8, 12 | Percent change from baseline |
| Endoscopic improvement | Week 12 | SES-CD score reduction >=50% |
| Quality of life | Week 12 | IBDQ score change from baseline |
| Febuxostat 40mg response | Week 12 | CDAI-100 response (exploratory) |

### Safety Endpoints

- Adverse events (graded by CTCAE v5.0)
- Liver function tests (ALT, AST, total bilirubin)
- Serum uric acid levels
- Complete blood count (CBC)
- Cardiovascular events

---

## G. Analysis Populations

| Population | Definition | Use |
|------------|------------|-----|
| **ITT (Intent-to-Treat)** | All randomized patients | PRIMARY analysis |
| **Per-Protocol** | Patients completing >=80% of doses without major protocol deviations | Supportive analysis |
| **Safety** | All patients receiving >=1 dose of study drug | Safety analysis |

### Missing Data Handling

- Primary analysis: Non-responder imputation (NRI) for missing Week 12 CDAI
- Sensitivity analyses: Multiple imputation, tipping point analysis

---

## H. Safety Monitoring

### Data Safety Monitoring Board (DSMB)

- Independent DSMB reviews safety data at 50% enrollment
- **No formal efficacy interim analysis** (avoids alpha spending complexity)
- DSMB membership: 2 gastroenterologists, 1 biostatistician, 1 clinical pharmacologist

### Stopping Rules

| Trigger | Action |
|---------|--------|
| >15% Grade 3+ hepatic AEs in any arm | Pause enrollment, DSMB review |
| Any treatment-related death | Halt trial, emergency DSMB review |
| Clear evidence of harm (DSMB judgment) | Halt trial |

### Laboratory Monitoring Schedule

| Assessment | Screening | Baseline | Week 4 | Week 8 | Week 12 | Week 16 (Follow-up) |
|------------|-----------|----------|--------|--------|---------|---------------------|
| ALT, AST, Bilirubin | X | X | X | X | X | X |
| Serum uric acid | X | X | X | X | X | X |
| CBC | X | X | X | X | X | X |
| CRP | | X | X | X | X | |
| Fecal calprotectin | | X | X | X | X | |
| CDAI | | X | X | X | X | X |
| NOD2 genotyping | X | | | | | |

---

## I. Study Timeline

| Phase | Duration | Activities |
|-------|----------|------------|
| Start-up | Months 1-3 | Site selection, IRB approval, drug supply |
| Enrollment | Months 4-15 | Patient screening, genotyping, randomization |
| Treatment | Rolling | 12-week treatment per patient |
| Follow-up | Rolling | 4-week post-treatment follow-up |
| Analysis | Months 16-18 | Database lock, statistical analysis |
| Reporting | Months 19-21 | Manuscript preparation, regulatory submission |

**Total estimated duration:** 21 months

---

## J. Regulatory Considerations

- IND application to FDA (Febuxostat new indication)
- IRB/Ethics committee approval at all sites
- Informed consent with genetic testing disclosure
- Data privacy compliance (HIPAA, GDPR as applicable)

---

## K. References

1. Hugot JP et al. (2001) Nature 411:599-603
2. Ogura Y et al. (2001) Nature 411:603-606
3. Febuxostat prescribing information (Takeda)
4. CDAI scoring: Best WR et al. Gastroenterology 1976;70:439-444
5. SES-CD: Daperno M et al. Gastrointest Endosc 2004;60:505-512
