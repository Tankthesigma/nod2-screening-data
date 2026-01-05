#!/usr/bin/env python3
"""Run ONLY ursodiol_rep3 on 5090"""
import os
os.chdir('/workspace/nod2-screening-data/PHASE_6')

# Quick inline simulation
from openmm import *
from openmm.app import *
from openmm.unit import *
import sys
sys.path.insert(0, 'gpu_scripts/split_2gpu')
from vast_5090 import run_simulation, setup_directories, preflight_check, SIMULATIONS

# Override to just ursodiol_rep3
SIMULATIONS = [("ursodiol", 3, "NOD2_LRR_clean.pdb", "ursodiol_docked.sdf")]

setup_directories()
print('Running ursodiol_rep3 ONLY')
run_simulation('ursodiol', 3, 'NOD2_LRR_clean.pdb', 'ursodiol_docked.sdf', 1, 1)
print('DONE')
