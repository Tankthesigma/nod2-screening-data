# GPU RENTAL INSTRUCTIONS FOR NOD2-SCOUT MD SIMULATIONS

## Overview
- **Total simulation time**: ~220 ns
- **Estimated GPU time**: 6-12 hours on RTX 4090
- **Recommended providers**: Vast.ai, Lambda Labs, RunPod

## Step 1: Choose a Provider

### Option A: Vast.ai (Recommended - Cheapest)
1. Go to https://vast.ai/
2. Create account and add credits ($10-20 should be enough)
3. Search for: RTX 4090, 24GB VRAM, Ubuntu 22.04
4. Look for instances with good network speed (>100 Mbps)
5. Rent for 12-24 hours (safety margin)

**Estimated cost**: $0.40-0.80/hour × 12h = $5-10

### Option B: Lambda Labs
1. Go to https://lambdalabs.com/
2. Select "1x A10" or "1x A100" instance
3. Choose PyTorch image

**Estimated cost**: $0.75/hour × 12h = $9

### Option C: RunPod
1. Go to https://runpod.io/
2. Select GPU pod with RTX 4090 or A100
3. Choose "Community Cloud" for lower prices

**Estimated cost**: $0.50-1.00/hour × 12h = $6-12

---

## Step 2: Set Up Environment

Once connected to your GPU instance, run:

```bash
# Update system
sudo apt update && sudo apt upgrade -y

# Install conda (if not present)
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda3
source $HOME/miniconda3/bin/activate

# Create environment
conda create -n nod2md python=3.10 -y
conda activate nod2md

# Install OpenMM with CUDA
conda install -c conda-forge openmm cudatoolkit=11.8 -y

# Install other dependencies
pip install mdtraj pdbfixer numpy pandas matplotlib seaborn scipy

# Verify CUDA
python -c "from openmm import Platform; print(Platform.getPlatformByName('CUDA').getName())"
```

---

## Step 3: Upload Files

Upload the following from your local machine:

```bash
# From local machine
scp -r md_simulation/structures/ user@gpu-instance:~/nod2/
scp -r md_simulation/gpu_scripts/ user@gpu-instance:~/nod2/
```

Or use the provider's file upload interface.

---

## Step 4: Run Simulations

```bash
cd ~/nod2

# 1. Prepare structures (if not done)
python gpu_scripts/01_prepare_receptor.sh

# 2. Run MD for each compound
for compound in febuxostat ursodiol budesonide; do
    for rep in 0 1 2; do
        echo "Running $compound replicate $rep"
        python gpu_scripts/03_run_md_$compound.py \
            --pdb structures/complex_$compound.pdb \
            --replicate $rep
    done
done

# 3. Run analysis
for compound in febuxostat ursodiol budesonide; do
    python gpu_scripts/04_analyze_$compound.py \
        --traj ${compound}_rep1.dcd \
        --top structures/complex_$compound.pdb
done
```

---

## Step 5: Download Results

```bash
# From GPU instance
tar -czvf md_results.tar.gz *.dcd *.log *_analysis.json

# From local machine
scp user@gpu-instance:~/nod2/md_results.tar.gz ./
```

---

## Troubleshooting

### "CUDA out of memory"
- Reduce system size or use implicit solvent
- Try mixed precision: `CudaPrecision: 'mixed'`

### "Platform not found"
- Reinstall: `conda install -c conda-forge openmm cudatoolkit=11.8`
- Check CUDA: `nvidia-smi`

### Slow performance
- Verify GPU is being used: Check for "CUDA" in output
- Close other processes: `nvidia-smi` then `kill <PID>`

---

## Time Estimates (RTX 4090)

| System | Time/replicate | Total (3 rep) |
|--------|---------------|---------------|
| Febuxostat | ~2 hours | ~6 hours |
| Ursodiol | ~2 hours | ~6 hours |
| Budesonide | ~2 hours | ~6 hours |
| Apo (control) | ~1.5 hours | ~4.5 hours |
| **TOTAL** | - | **~22.5 hours** |

*Note: Times may vary based on system size and GPU model*

---

## Budget Summary

| Provider | Cost/hour | 24h rental |
|----------|-----------|------------|
| Vast.ai | $0.40-0.80 | $10-20 |
| Lambda | $0.75 | $18 |
| RunPod | $0.50-1.00 | $12-24 |

**Recommendation**: Vast.ai with RTX 4090 for best price/performance.

---

Generated: 2026-01-02 15:19:09
NOD2-Scout Phase 6 MD Simulation Pipeline
