===============================================
FEP VAST.AI DEPLOYMENT - UPDATED
===============================================

TOTAL: 47 windows remaining
- Your local GPU: 5 windows (wt_complex 00, 14-17)
- Vast.ai (8 GPUs): 42 windows

-----------------------------------------------
STEP 1: PACKAGE DATA FOR UPLOAD
-----------------------------------------------
In PowerShell:

  cd C:\Users\vasud\nod2-screening-data
  tar -cvzf fep_data.tar.gz fep_pmx vast_deploy

Creates ~500MB file.

-----------------------------------------------
STEP 2: RENT VAST.AI INSTANCE
-----------------------------------------------
- Rent 8x RTX 4070S ($0.65/hr) or similar
- Select "PyTorch" template (has CUDA)
- Note the SSH command they give you

-----------------------------------------------
STEP 3: UPLOAD TO VAST.AI
-----------------------------------------------
From PowerShell (replace PORT and IP):

  scp -P <PORT> fep_data.tar.gz root@<VAST_IP>:/workspace/

-----------------------------------------------
STEP 4: SETUP ON VAST.AI (SSH in first)
-----------------------------------------------
  cd /workspace
  tar -xvzf fep_data.tar.gz
  mkdir -p fep/data
  mv fep_pmx fep/data/
  mv vast_deploy/* fep/
  cd fep
  chmod +x *.sh

  # Install dependencies ONCE (not per GPU!)
  ./setup_once.sh

-----------------------------------------------
STEP 5: LAUNCH ALL 8 GPUS
-----------------------------------------------
  ./run_all_gpus.sh

Monitor:
  tail -f gpu_0.log           # Watch one GPU
  grep 'DONE' gpu_*.log       # See completions
  grep 'COMPLETE' gpu_*.log   # See finished GPUs

-----------------------------------------------
STEP 6: RUN LOCALLY (SAME TIME)
-----------------------------------------------
Double-click: run_local.bat

Or:
  cd C:\Users\vasud\nod2-screening-data\vast_deploy
  run_local.bat

-----------------------------------------------
STEP 7: DOWNLOAD RESULTS
-----------------------------------------------
After all GPUs show COMPLETE (~2 hours):

On Vast.ai:
  cd /workspace/fep/data
  find fep_pmx -name "u_nk.npy" | tar -cvzf results.tar.gz -T -

Download to local:
  scp -P <PORT> root@<VAST_IP>:/workspace/fep/data/results.tar.gz .

-----------------------------------------------
STEP 8: STOP INSTANCE!
-----------------------------------------------
Don't forget or you keep paying!

===============================================
WINDOW DISTRIBUTION
===============================================

Your GPU (5 windows):
  wt_complex: 00, 14, 15, 16, 17

Vast.ai GPU 0 (5 windows):
  wt_complex: 18, 19
  mut_complex: 00, 01, 02

Vast.ai GPU 1-3 (5 each):
  mut_complex: 03-17

Vast.ai GPU 4 (5 windows):
  mut_complex: 18, 19
  solvent: 00, 01, 02

Vast.ai GPU 5-6 (5 each):
  solvent: 03-12

Vast.ai GPU 7 (7 windows):
  solvent: 13-19

===============================================
ESTIMATED TIME & COST
===============================================
- Vast.ai: ~2 hours x $0.65 = ~$1.30
- Your GPU: ~1.5 hours (5 windows)
- Total: < $2

===============================================
