# ISEF 2026 Talking Points: NOD2-CROHN Drug Discovery

> **For judges and presentation**

---

## Four Key Innovation Points

### 1. Precision Medicine: Genotype-Guided Patient Selection

**The Problem:** Current Crohn's therapies treat all patients the same, ignoring genetic differences.

**Our Solution:** We use NOD2 genotyping to select patients most likely to benefit.

- **INCLUDE:** Patients with WT, R702W, G908R, N852S, or M863V genotypes (binding pocket intact)
- **EXCLUDE:** L1007fs carriers (frameshift mutation deletes pocket residues GLU1008, ASP1011, ARG1037)

**Why it matters:** This is the **first genotype-stratified therapy proposal** for Crohn's disease, moving from "one-size-fits-all" to personalized medicine.

---

### 2. Drug Repurposing: 10+ Years of Safety Data

**The Problem:** New drug development takes 10-15 years and costs $1-2 billion with high failure rates.

**Our Solution:** Repurpose Febuxostat, an FDA-approved gout medication.

- **Already proven safe:** 10+ years of post-marketing safety data
- **Known pharmacology:** Well-characterized ADMET profile
- **Faster path:** Phase II trial (2-3 years) vs new drug (10+ years)
- **Lower cost:** ~$10M vs ~$1B for de novo development

**Why it matters:** Repurposing dramatically reduces risk, cost, and time to patients.

---

### 3. Cost Impact: Making Treatment Accessible

**The Problem:** Current biologic therapies cost $20,000-$30,000+ per year, limiting access worldwide.

**Our Solution:** Febuxostat is available as a generic oral tablet.

| Therapy Type | Approximate Annual Cost |
|--------------|------------------------|
| Biologics (Humira, Remicade, etc.) | $20,000 - $30,000+ |
| **Febuxostat (generic)** | **<$1,000** |

**Potential savings:** >95% cost reduction

**Why it matters:** Affordable oral therapy could transform Crohn's treatment globally, especially in resource-limited settings.

---

### 4. Scientific Innovation: First NOD2-Targeted Small Molecule

**The Problem:** All current Crohn's therapies target downstream inflammation (TNF-alpha, interleukins), not the root cause.

**Our Solution:** Directly stabilize the NOD2 protein where the genetic defect occurs.

**Computational Evidence:**
- Virtual screening of 11,000+ compounds identified Febuxostat
- Molecular dynamics simulation: 99-100% pocket occupancy (2/3 replicates)
- Distance analysis confirms missense mutations don't affect binding site

**Why it matters:** This is a **paradigm shift** from treating symptoms to addressing the underlying molecular defect.

---

## One-Sentence Summary

> "We designed the first genotype-stratified clinical trial for Crohn's disease, repurposing a safe, affordable generic drug to target the NOD2 protein based on computational drug discovery and molecular dynamics validation."

---

## Anticipated Judge Questions & Answers

**Q: Why exclude L1007fs patients?**
A: The L1007fs frameshift truncates the protein at residue 1007, deleting the binding pocket residues (GLU1008, ASP1011, ARG1037). Without these residues, the drug has nowhere to bind. Our distance analysis confirmed this computationally.

**Q: How confident are you in the computational predictions?**
A: We validated our docking with 60 ns molecular dynamics simulations showing stable binding (99-100% pocket occupancy in 2/3 replicates). However, computational predictions must be confirmed experimentally - that's why we designed a rigorous Phase II trial.

**Q: What if the trial fails?**
A: Drug development has inherent uncertainty. We estimate 40-50% probability of a positive Phase II result. Even a negative result would be scientifically valuable, ruling out this mechanism and informing future research.

**Q: How is this different from existing research?**
A: To our knowledge, no one has proposed a genotype-stratified, NOD2-targeted small molecule therapy for Crohn's disease. We combined computational screening, MD validation, and clinical trial design into a complete translational pipeline.

**Q: What about the natural compound (CID_10120 Bufadienolide)?**
A: Our virtual screening also identified CID_10120, a bufadienolide (cardiac glycoside class), as a strong binder. **NOTE:** This compound was previously mislabeled as "CID_10592 / Dihydrocortisol" - the correct identity is CID_10120. Bufadienolides are known NF-κB inhibitors, and since NOD2 signals through NF-κB, this suggests a dual mechanism. However, unlike Febuxostat, this is a novel therapeutic candidate requiring full preclinical development. We prioritize Febuxostat because its FDA-approved status enables a direct Phase II trial.

---

## Visual Aids to Prepare

1. **NOD2 Structure:** 3D image showing binding pocket + mutation locations
2. **Patient Flow:** Screening -> Genotyping -> Eligibility -> Randomization
3. **Cost Comparison:** Bar chart of biologic vs febuxostat costs
4. **MD Simulation:** Pocket occupancy over time (stable binding)
5. **Precision Medicine Diagram:** Genotype decision tree

---

## Contact for Collaboration

*[Your information here]*

**Project Repository:** github.com/Tankthesigma/nod2-screening-data
