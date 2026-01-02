# PHASE 6: MOLECULAR DYNAMICS SIMULATION
## Easy Step-by-Step Guide for NOD2-Scout (ISEF 2026)

---

## WHAT IS THIS?

Phase 6 validates that your #1 drug candidate (Febuxostat) actually **sticks** to the NOD2 protein. We simulate the protein-drug complex for 20 nanoseconds to see if the drug stays bound or floats away.

**Think of it like this:** Docking (Phase 4) was a snapshot. MD simulation is a movie showing if the drug stays in place.

---

## WHAT YOU NEED

| Item | Cost | Where to Get |
|------|------|--------------|
| GPU Instance | $10-20 | Vast.ai (cheapest) |
| ~2 hours | Free | Your time |

---

## THE FILES EXPLAINED

### 📁 `structures/` - The Molecules (UPLOAD THIS)

| File | What It Is |
|------|-----------|
| `complex_febuxostat.pdb` | NOD2 protein + Febuxostat (YOUR #1 CANDIDATE) |
| `complex_ursodiol.pdb` | NOD2 protein + Ursodiol (backup candidate) |
| `complex_budesonide.pdb` | NOD2 protein + Budesonide (KNOWN CROHN'S DRUG - positive control) |
| `complex_decoy.pdb` | NOD2 protein + bad compound (should fail - negative control) |
| `complex_apo.pdb` | NOD2 protein alone (no drug - baseline) |
| `receptor.pdb` | Just the NOD2 protein |
| `ligand_*.pdb` | Individual drug molecules |

### 📁 `gpu_scripts/` - The Code (UPLOAD THIS)

| File | What It Does |
|------|-------------|
| `master_md_pipeline.sh` | **RUN THIS** - Does everything automatically |
| `03_run_md_febuxostat.py` | Simulates Febuxostat binding (20ns x 3 runs) |
| `03_run_md_ursodiol.py` | Simulates Ursodiol binding |
| `03_run_md_budesonide.py` | Simulates Budesonide (control) |
| `03_run_md_apo.py` | Simulates empty protein |
| `04_analyze_*.py` | Analyzes if drug stayed bound |
| `05_mmgbsa_*.py` | Calculates binding strength |

### 📁 `figures/` - The Results (ALREADY MADE - Placeholders)

8 publication-ready figures will be regenerated with real data after simulation.

---

## STEP-BY-STEP INSTRUCTIONS

### STEP 1: Rent a GPU (5 minutes)

1. Go to **https://vast.ai**
2. Create account, add $20 credits
3. Click "Search" and filter:
   - GPU: RTX 4090
   - Disk: 50GB+
   - Image: nvidia/cuda or pytorch/pytorch
4. Click "RENT" on cheapest option (~$0.40/hour)
5. Wait for "Running" status, then click "Connect"

### STEP 2: Set Up the GPU (10 minutes)

Copy-paste these commands ONE BY ONE:

```bash
# 1. Update system
sudo apt update && sudo apt install -y wget git

# 2. Install Miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda3
source $HOME/miniconda3/bin/activate

# 3. Create environment with OpenMM
conda create -n nod2md python=3.10 -y
conda activate nod2md
conda install -c conda-forge openmm cudatoolkit=11.8 -y
pip install mdtraj numpy pandas matplotlib seaborn

# 4. Verify GPU works
python -c "from openmm import Platform; print('GPU:', Platform.getPlatformByName('CUDA').getName())"
```

You should see: `GPU: CUDA`

### STEP 3: Download Your Files (2 minutes)

```bash
# Clone your repo
git clone https://github.com/Tankthesigma/nod2-screening-data.git
cd nod2-screening-data/md_simulation

# Check files are there
ls structures/
ls gpu_scripts/
```

You should see all the `.pdb` and `.py` files.

### STEP 4: Run the Simulation (12-24 hours)

**Option A: Run Everything Automatically**
```bash
chmod +x gpu_scripts/master_md_pipeline.sh
./gpu_scripts/master_md_pipeline.sh
```

**Option B: Run Just Febuxostat First (faster test)**
```bash
python gpu_scripts/03_run_md_febuxostat.py --pdb structures/complex_febuxostat.pdb --replicate 0
```

### STEP 5: Download Results

When done, package and download:

```bash
# On GPU - package results
tar -czvf md_results.tar.gz trajectories/ analysis/ *.log

# On YOUR computer - download
scp user@gpu-ip:~/nod2-screening-data/md_simulation/md_results.tar.gz ~/Downloads/
```

---

## WHAT THE SIMULATION DOES

```
1. MINIMIZATION (5 min)
   └── Relaxes the structure, removes bad contacts

2. HEATING (10 min)
   └── Slowly heats from 0K to 310K (body temperature)

3. EQUILIBRATION (30 min)
   └── Lets the system stabilize at body temperature

4. PRODUCTION (2-3 hours per compound)
   └── Records 20 nanoseconds of motion
   └── Saves snapshots every 10 picoseconds
   └── Runs 3 times for statistics
```

---

## WHAT SUCCESS LOOKS LIKE

### Good Result (Drug Stays Bound):
- RMSD < 2.5 Å (low movement)
- H-bonds: 5+ persistent
- Distance to LEU1007: < 5 Å

### Bad Result (Drug Falls Off):
- RMSD > 4 Å (high movement)
- H-bonds: < 2
- Distance to LEU1007: > 10 Å

---

## THE 8 FIGURES EXPLAINED

| Figure | What It Shows |
|--------|--------------|
| Fig 1: RMSD | Does the drug wiggle too much? Lower = better |
| Fig 2: RMSF | Which parts of protein move most? |
| Fig 3: H-bonds | How many hydrogen bonds form? More = stronger binding |
| Fig 4: 1007fs Distance | Is drug near the Crohn's mutation site? |
| Fig 5: MM-GBSA | Binding energy in kcal/mol. More negative = stronger |
| Fig 6: Key Residues | Distance to important amino acids |
| Fig 7: Heatmap | All metrics compared at once |
| Fig 8: Dashboard | Final summary for your poster |

---

## TROUBLESHOOTING

### "CUDA not found"
```bash
conda install -c conda-forge cudatoolkit=11.8 -y
```

### "Out of memory"
- Close other processes: `nvidia-smi` then `kill <PID>`
- Or rent a GPU with more VRAM (A100 has 80GB)

### "Module not found"
```bash
pip install mdtraj openmm pdbfixer
```

### Simulation crashes
- Check log files: `cat logs/*.log`
- Usually means bad structure - email me

---

## TIME & COST ESTIMATE

| Compound | Time (RTX 4090) |
|----------|-----------------|
| Febuxostat | ~6 hours (3 x 2h) |
| Ursodiol | ~6 hours |
| Budesonide | ~6 hours |
| Apo | ~4 hours |
| **TOTAL** | **~22 hours** |

**Cost on Vast.ai:** $0.40/hr × 24h = **$10-15**

---

## QUICK REFERENCE COMMANDS

```bash
# Activate environment
conda activate nod2md

# Run all simulations
./gpu_scripts/master_md_pipeline.sh

# Run single compound
python gpu_scripts/03_run_md_febuxostat.py --pdb structures/complex_febuxostat.pdb --replicate 0

# Check progress
tail -f logs/febuxostat_rep1.log

# Analyze results
python gpu_scripts/04_analyze_febuxostat.py --traj trajectories/febuxostat_rep1.dcd --top structures/complex_febuxostat.pdb
```

---

## FILE STRUCTURE AFTER COMPLETION

```
md_simulation/
├── structures/          # Input files (already done)
├── gpu_scripts/         # Scripts (already done)
├── trajectories/        # OUTPUT: .dcd movie files
│   ├── febuxostat_rep1.dcd
│   ├── febuxostat_rep2.dcd
│   └── ...
├── analysis/            # OUTPUT: .json results
│   ├── febuxostat_analysis.json
│   └── ...
├── logs/                # OUTPUT: progress logs
│   ├── febuxostat_rep1.log
│   └── ...
└── figures/             # OUTPUT: regenerated figures
    ├── fig1_rmsd_comparison.png
    └── ...
```

---

## NEED HELP?

1. Check the log files first: `cat logs/*.log`
2. Google the error message
3. Open an issue on GitHub

---

**Created:** 2026-01-02
**Author:** NOD2-Scout Pipeline
**Project:** ISEF 2026 - Crohn's Disease Drug Discovery
