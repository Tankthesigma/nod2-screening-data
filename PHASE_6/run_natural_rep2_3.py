#!/usr/bin/env python3
"""Run natural_top rep2 and rep3 on other GPU"""
import os
os.chdir('/workspace/nod2-screening-data/PHASE_6')
import sys
sys.path.insert(0, 'gpu_scripts/split_2gpu')
from vast_other import run_simulation, setup_directories

setup_directories()
print('Running natural_top rep2 and rep3')
run_simulation('natural_top', 2, 'NOD2_LRR_clean.pdb', 'natural_top_docked.sdf', 1, 2)
run_simulation('natural_top', 3, 'NOD2_LRR_clean.pdb', 'natural_top_docked.sdf', 2, 2)
print('DONE')
