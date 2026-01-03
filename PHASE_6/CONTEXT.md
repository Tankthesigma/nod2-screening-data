# NOD2-SCOUT PHASE 6 - MD SIMULATION CONTEXT

## WHAT THIS IS
ISEF 2026 project - drug screening for Crohn's disease targeting NOD2 protein.
Phase 6 = Molecular Dynamics simulations to validate docked drug candidates.

## CURRENT STATUS: DEBUGGING NVT CRASH

### The Setup (2-GPU split)
- **vast_5090.py**: Runs on rented RTX 5090 (10 sims, 200ns) - auto-shutdown when done
- **local_4060.py**: Runs on local RTX 4060 Ti (4 sims, 80ns) - no shutdown

### WHAT'S BEEN FIXED SO FAR
1. ✅ HMR (Hydrogen Mass Repartitioning) - 3.0 amu for 4fs timestep
2. ✅ Constraints changed from AllBonds → HBonds
3. ✅ Position restraints use minimized coords (not original)
4. ✅ Dispersion correction added for NPT density
5. ✅ Barostat frequency 25 → 100 steps
6. ✅ SDF loading handles multi-conformer files (returns list)
7. ✅ PDBFixer integration for missing atoms
8. ✅ Fresh AlphaFold receptor download (corrupted PDBs replaced)
9. ✅ OpenFF SystemGenerator fix - nonbondedMethod in forcefield_kwargs
10. ✅ Gradual heating NVT: 100K → 200K → 310K
11. ✅ Lower restraint force: 1000 → 100 kJ/mol/nm²
12. ✅ Re-minimize after adding restraints

### CURRENT ERROR (if still crashing)
```
Particle coordinate is NaN
```
This happens during NVT equilibration. Means atoms are flying off to infinity.

### POSSIBLE REMAINING ISSUES
1. **Ligand not in binding site** - receptor downloaded fresh from AlphaFold, ligand SDF has docked coordinates. They might not align!
2. **Ligand clashing with receptor** - need to check if docked pose is compatible with AlphaFold structure
3. **Missing ligand in system** - receptor-only system being created, ligand not added

### KEY FILES
```
PHASE_6/
├── gpu_scripts/split_2gpu/
│   ├── vast_5090.py      # Vast.ai script (10 sims)
│   ├── local_4060.py     # Local 4060 Ti script (4 sims)
│   ├── vast_launch.sh    # Vast.ai launcher
│   └── local_launch.sh   # Local launcher
├── download_receptor.py   # Downloads fresh NOD2 from AlphaFold
├── fix_receptor.py        # Fixes corrupted PDBs (deprecated)
└── structures/
    ├── NOD2_LRR_clean.pdb        # Clean receptor (from AlphaFold)
    ├── febuxostat_docked.sdf     # Docked ligand poses
    ├── ursodiol_docked.sdf
    ├── budesonide_docked.sdf
    ├── natural_top_docked.sdf
    └── decoy_docked.sdf
```

### SIMULATION PROTOCOL
1. Load receptor PDB (NOD2_LRR_clean.pdb)
2. Load ligand SDF
3. PDBFixer cleans receptor
4. SystemGenerator parameterizes (OpenFF → GAFF fallback)
5. Solvate with TIP3P water + 0.15M ions
6. Minimize (5000 steps)
7. NVT with restraints + gradual heating (100ps)
8. NPT unrestrained (500ps)
9. Production (20ns)

### TO DEBUG NEXT
If still crashing:
1. Check if ligand coordinates in SDF match receptor binding site
2. Try running APO (no ligand) first - if that works, problem is ligand
3. Visualize receptor + ligand in PyMOL to check alignment
4. May need to re-dock ligands to the AlphaFold structure

### GITHUB
https://github.com/Tankthesigma/nod2-screening-data

### COMMANDS TO RUN

**On Vast.ai:**
```bash
export GITHUB_TOKEN=<your_token_here>
git clone https://github.com/Tankthesigma/nod2-screening-data.git
cd nod2-screening-data/PHASE_6
bash gpu_scripts/split_2gpu/vast_launch.sh
```

**On Local 4060 Ti:**
```bash
git pull
cd PHASE_6
bash gpu_scripts/split_2gpu/local_launch.sh
```

---
Last updated: 2026-01-03
