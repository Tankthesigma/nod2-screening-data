#!/bin/bash
# Vast.ai setup script for FEP windows
# Run this on each Vast.ai instance after copying files

# Install miniconda if not present
if ! command -v conda &> /dev/null; then
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
    bash miniconda.sh -b -p $HOME/miniconda3
    eval "$($HOME/miniconda3/bin/conda shell.bash hook)"
    conda init
fi

# Create fep environment
conda create -n fep python=3.10 -y
conda activate fep

# Install OpenMM and dependencies
conda install -c conda-forge openmm cudatoolkit=11.8 -y
pip install openmmtools pymbar numpy

echo "Setup complete. Run windows with:"
echo "  CUDA_VISIBLE_DEVICES=0 python run_failed_windows_fresh.py --system <sys> --window <idx>"
