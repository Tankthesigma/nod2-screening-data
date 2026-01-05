import subprocess
import sys
sys.path.insert(0, 'gpu_scripts/split_2gpu')

# Patch the script to only run ursodiol_rep3
import vast_5090
vast_5090.SIMULATIONS = [("ursodiol", 3, "NOD2_LRR_clean.pdb", "ursodiol_docked.sdf")]
vast_5090.main()
