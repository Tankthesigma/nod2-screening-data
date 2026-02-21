import os
import sys
import glob
import pickle
import shutil
import time
import numpy as np

# Force unbuffered output
sys.stdout.reconfigure(line_buffering=True)
sys.stderr.reconfigure(line_buffering=True)
os.environ["PYTHONUNBUFFERED"] = "1"

print("="*60)
print("FEP SETUP SCRIPT (FIXED)")
print("="*60)

print()
print("Working directory:", os.getcwd())
print("[WARNING] Run from project root!")

try:
    from rdkit import Chem
    print("[PASS] RDKit")
except ImportError:
    print("[FAIL] RDKit! Install: conda install -c conda-forge rdkit")
    sys.exit(1)

import mdtraj as md
from openmm import XmlSerializer, LangevinMiddleIntegrator, Platform
from openmm.app import Simulation, PDBFile, Topology
import openmm.unit as unit
import openmm

print()
print("Versions:")
print("  OpenMM:", openmm.__version__)

try:
    import openmmtools
    print("  OpenMMTools:", openmmtools.__version__)
except:
    try:
        import openmmtools
        print("  OpenMMTools: (version unavailable)")
    except:
        print("  OpenMMTools: MISSING")
        sys.exit(1)

try:
    import perses
    print("  Perses:", perses.__version__)
except:
    try:
        import perses
        print("  Perses: (version unavailable)")
    except:
        print("  Perses: MISSING")
        sys.exit(1)

# Try both old and new perses import paths
try:
    from perses.app.relative_point_mutation_setup import PointMutationExecutor
    print("  Perses API: relative_point_mutation_setup")
except ImportError:
    try:
        from perses.app.setup_relative_calculation import PointMutationExecutor
        print("  Perses API: setup_relative_calculation")
    except ImportError:
        print("[FAIL] Cannot import PointMutationExecutor from perses!")
        sys.exit(1)

print()
print("="*60)
print("PART 1: FINDING FILES")
print("="*60)

# FIXED: Added PHASE_A1_mutant_MD to PDB patterns for consistency
dcd_patterns = ["**/WT*[Ff]ebuxostat*rep2*.dcd", "**/WT*[Ff]ebuxostat*.dcd",
                "**/febuxostat*rep2*.dcd", "PHASE_6/**/*.dcd", "PHASE_A1_mutant_MD/**/*.dcd"]
pdb_patterns = ["**/WT*[Ff]ebuxostat*rep2*solvated*.pdb", "**/WT*[Ff]ebuxostat*rep2*.pdb",
                "**/WT*[Ff]ebuxostat*solvated*.pdb", "PHASE_6/**/*solvated*.pdb",
                "PHASE_A1_mutant_MD/**/*solvated*.pdb", "PHASE_A1_mutant_MD/**/*.pdb"]
sdf_patterns = ["**/febuxostat.sdf", "**/[Ff]ebuxostat*.sdf", "**/ligand*.sdf", "PHASE_6/**/*.sdf"]

def find_files(patterns):
    found = []
    for p in patterns:
        found.extend(glob.glob(p, recursive=True))
    # FIXED: Sort for deterministic file selection
    return sorted(set(found))

dcd_files = find_files(dcd_patterns)
pdb_files = find_files(pdb_patterns)
sdf_files = find_files(sdf_patterns)

print("Found:", len(dcd_files), "DCD,", len(pdb_files), "PDB,", len(sdf_files), "SDF")

if not dcd_files:
    print("[FAIL] No DCD!")
    sys.exit(1)

# PRIORITY: febuxostat WT trajectories (exclude mutants R702W, G908R)
def is_wt_febuxostat(f):
    fl = f.lower()
    return "febuxostat" in fl and "r702w" not in fl and "g908r" not in fl and "mutant" not in fl

wt_febuxostat = [f for f in dcd_files if is_wt_febuxostat(f)]
print(f"[INFO] Found {len(wt_febuxostat)} WT Febuxostat trajectories")

if wt_febuxostat:
    # Prefer rep2, then largest file
    rep2_feb = [f for f in wt_febuxostat if "rep2" in f.lower()]
    if rep2_feb:
        dcd_path = max(rep2_feb, key=os.path.getsize)
    else:
        dcd_path = max(wt_febuxostat, key=os.path.getsize)
else:
    print("[FAIL] No WT Febuxostat trajectories found!")
    print("  Available DCDs:")
    for f in dcd_files[:10]:
        print(f"    {os.path.basename(f)}")
    sys.exit(1)

print("[SELECTED] DCD:", os.path.basename(dcd_path))
if "mutant" in dcd_path.lower() or "r702w" in dcd_path.lower():
    print("[FAIL] Selected DCD appears to be mutant trajectory! Cannot do R702W FEP on R702W structure.")
    sys.exit(1)

# FIXED: More robust PDB matching - collect all viable PDBs and select best
viable_pdbs = []
for pdb in pdb_files:
    try:
        test_frame = md.load_frame(dcd_path, 0, top=pdb)
        if test_frame.n_atoms > 1000:
            viable_pdbs.append({"path": pdb, "n_atoms": test_frame.n_atoms})
    except:
        pass

if not viable_pdbs:
    print("[FAIL] No PDB matches DCD!")
    sys.exit(1)

if len(viable_pdbs) > 1:
    print(f"[WARNING] {len(viable_pdbs)} PDBs match DCD:")
    for p in viable_pdbs:
        print(f"  {os.path.basename(p['path'])}: {p['n_atoms']} atoms")
    # Prefer WT febuxostat, then largest atom count
    def pdb_score(p):
        name = p["path"].lower()
        wt_bonus = 100 if "wt" in name and "febuxostat" in name else 0
        return wt_bonus + p["n_atoms"]
    viable_pdbs.sort(key=pdb_score, reverse=True)
    print("  Using:", os.path.basename(viable_pdbs[0]["path"]))

pdb_path = viable_pdbs[0]["path"]
print("[SELECTED] PDB:", os.path.basename(pdb_path), "(", viable_pdbs[0]["n_atoms"], "atoms)")

if not sdf_files:
    print("[FAIL] No SDF!")
    sys.exit(1)

print()
print("="*60)
print("PART 2: PREPARING STRUCTURES")
print("="*60)

for d in ["fep_study", "fep_study/inputs", "fep_study/complex", "fep_study/protein",
          "fep_study/results", "fep_study/figures"]:
    os.makedirs(d, exist_ok=True)

# FIXED: Use load_frame to avoid OOM on large trajectories
print("[INFO] Loading last frame only (memory-safe)...")
try:
    n_frames = len(md.open(dcd_path))
    print("  Total frames in trajectory:", n_frames)
    if n_frames == 0:
        print("[FAIL] DCD has 0 frames!")
        sys.exit(1)
    last_frame = md.load_frame(dcd_path, n_frames - 1, top=pdb_path)
    print("  Last frame atoms:", last_frame.n_atoms)
except Exception as e:
    print("[FAIL] Could not load DCD:", str(e)[:100])
    print("  DCD:", dcd_path)
    print("  PDB:", pdb_path)
    sys.exit(1)

# FIXED: Image/center the trajectory to handle PBC
print("[INFO] Imaging trajectory to handle PBC...")
try:
    last_frame.image_molecules(inplace=True)
    print("  [PASS] Molecules imaged")
except Exception as e:
    print("  [WARNING] Imaging failed:", str(e)[:50], "- proceeding anyway")

def count_heavy(residue):
    count = 0
    for a in residue.atoms:
        if a.element is not None:
            if a.element.symbol != "H":
                count += 1
        elif not a.name.startswith("H"):
            count += 1
    return count

EXCLUDE = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE",
           "PRO","SER","THR","TRP","TYR","VAL","HIE","HID","HIP","HSD","HSE","HSP","CYX","CYM",
           "ACE","NME","NH2","HOH","WAT","TIP3","TIP4","SOL","H2O","OH2","TIP","NA","CL","K",
           "MG","CA","ZN","FE","CU","MN","NA+","CL-","K+","SOD","CLA","POT","HEM","NAG","GOL",
           "PEG","EDO","DMS","ACT","BME"}
LOW_PRI = {"SO4","PO4","CIT","TRS"}

res_counts = {}
for r in last_frame.topology.residues:
    n = r.name.upper()
    res_counts[n] = res_counts.get(n, 0) + 1

candidates = []
for r in last_frame.topology.residues:
    n = r.name.upper()
    if n not in EXCLUDE and res_counts[n] == 1:
        pri = 0 if n in LOW_PRI else 1
        candidates.append({"res": r, "name": r.name, "heavy": count_heavy(r), "pri": pri, "idx": r.index})

if not candidates:
    print("[FAIL] No ligand!")
    sys.exit(1)

candidates.sort(key=lambda x: (-x["pri"], -x["heavy"]))
ligand_res = candidates[0]["res"]
ligand_name = candidates[0]["name"]
ligand_heavy = candidates[0]["heavy"]
ligand_res_idx = candidates[0]["idx"]

print("[SELECTED] Ligand:", ligand_name, "(", ligand_heavy, "heavy)")

# Validate ligand name looks like febuxostat (or similar drug)
EXPECTED_LIGANDS = {"FBX", "FEB", "FEBUXOSTAT", "UNK", "UNL", "LIG"}
if ligand_name.upper() not in EXPECTED_LIGANDS:
    print(f"[WARNING] Ligand name '{ligand_name}' not in expected list: {EXPECTED_LIGANDS}")
    print("  Verify this is the correct ligand for FEP.")

# FIXED: Prioritize febuxostat SDF by name, then heavy atom count
febuxostat_sdfs = [s for s in sdf_files if "febuxostat" in s.lower()]
if febuxostat_sdfs:
    print(f"[INFO] Found {len(febuxostat_sdfs)} Febuxostat SDFs")
    # Use the first febuxostat SDF that matches heavy atom count
    best_sdf, best_diff = None, 999
    for sdf in febuxostat_sdfs:
        try:
            supplier = Chem.SDMolSupplier(sdf)
            for mol in supplier:
                if mol:
                    diff = abs(mol.GetNumHeavyAtoms() - ligand_heavy)
                    if diff < best_diff:
                        best_diff, best_sdf = diff, sdf
        except:
            pass
else:
    print("[WARNING] No Febuxostat SDFs found, using heavy atom matching")
    best_sdf, best_diff = None, 999
    for sdf in sdf_files:
        try:
            supplier = Chem.SDMolSupplier(sdf)
            for mol in supplier:
                if mol:
                    diff = abs(mol.GetNumHeavyAtoms() - ligand_heavy)
                    if diff < best_diff:
                        best_diff, best_sdf = diff, sdf
        except:
            pass

if not best_sdf:
    print("[FAIL] No valid SDF!")
    sys.exit(1)

print("[SELECTED] SDF:", os.path.basename(best_sdf), "(diff=", best_diff, ")")

# Warn if SDF has multiple molecules (perses only reads the first)
try:
    supplier = Chem.SDMolSupplier(best_sdf)
    mol_count = sum(1 for m in supplier if m)
    if mol_count > 1:
        print(f"[WARNING] SDF contains {mol_count} molecules! Perses will only use the first.")
        print("  If this SDF has multiple conformers/tautomers, verify the first is correct.")
except:
    pass

# FIXED: Handle multiple resSeq=702 matches, filter to protein only
PROTEIN_RES = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
               "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
               "HIE","HID","HIP","HSD","HSE","HSP","CYX","CYM","ARN","AR0"}

res702_matches = []
for chain in last_frame.topology.chains:
    for r in chain.residues:
        if r.resSeq == 702:
            res702_matches.append((chain, r))

# Filter to protein residues only (exclude water, ions)
protein_res702 = [(c, r) for c, r in res702_matches if r.name.upper() in PROTEIN_RES]

if len(res702_matches) > 1:
    print(f"[INFO] Found {len(res702_matches)} residues with resSeq=702:")
    for c, r in res702_matches:
        marker = " <- PROTEIN" if (c, r) in protein_res702 else " (non-protein, ignored)"
        print(f"  Chain {c.chain_id or c.index}: {r.name}{marker}")

if not protein_res702:
    print("[FAIL] No protein residue 702 found!")
    for c in last_frame.topology.chains:
        prot = [r for r in c.residues if r.name.upper() in PROTEIN_RES]
        if prot:
            print("  Chain", c.chain_id or c.index, ":", min(r.resSeq for r in prot), "-", max(r.resSeq for r in prot))
    sys.exit(1)

if len(protein_res702) > 1:
    print(f"[INFO] Multiple protein residues with resSeq=702 ({len(protein_res702)} found).")
    print("  Selecting chain closest to ligand...")
    
    # Calculate distance from each res702 to ligand centroid
    ligand_atoms = [a for a in last_frame.topology.atoms if a.residue.index == ligand_res_idx]
    ligand_coords = last_frame.xyz[0, [a.index for a in ligand_atoms], :] * 10  # nm to A
    ligand_centroid = np.mean(ligand_coords, axis=0)
    
    best_dist = float('inf')
    best_match = None
    for chain, res in protein_res702:
        res_atoms = [a for a in res.atoms]
        res_coords = last_frame.xyz[0, [a.index for a in res_atoms], :] * 10
        res_centroid = np.mean(res_coords, axis=0)
        dist = np.linalg.norm(res_centroid - ligand_centroid)
        print(f"    Chain {chain.chain_id or chain.index}, res {res.name}: {dist:.1f}A from ligand")
        if dist < best_dist:
            best_dist = dist
            best_match = (chain, res)
    
    target_chain, target_res = best_match
    print(f"  [SELECTED] Chain {target_chain.chain_id or target_chain.index} (closest at {best_dist:.1f}A)")
else:
    target_chain, target_res = protein_res702[0]

# FIXED: Check if already TRP - cannot do TRP->TRP mutation
if target_res.name.upper() in ["TRP"]:
    print("[FAIL] Residue 702 is already TRP! Cannot perform R702W mutation on mutant structure.")
    print("  This trajectory appears to be from a mutant simulation, not WT.")
    sys.exit(1)

raw_chain_id = target_chain.chain_id if target_chain.chain_id else None
index_str = str(target_chain.index)

chain_options_raw = []
if raw_chain_id:
    chain_options_raw.append(raw_chain_id)
if "A" not in chain_options_raw:
    chain_options_raw.append("A")
if index_str not in chain_options_raw:
    chain_options_raw.append(index_str)

chain_options = list(dict.fromkeys(chain_options_raw))

print("[INFO] Chain options to try:", chain_options)
print("[PASS] Residue 702:", target_res.name)

if target_res.name.upper() not in ["ARG","ARN","AR0"]:
    print("[FAIL] Expected ARG at 702, found", target_res.name, "- cannot proceed with R702W mutation")
    sys.exit(1)

print()
print("[WARNING] R->W: +1->0 charge change.")
print("  Perses should auto-add counterion for charge neutrality.")
print("  Will verify system charge after setup.")

# Include structural ions with protein (Mg, Zn, Ca, etc.)
STRUCTURAL_IONS = {"MG", "ZN", "CA", "MN", "FE", "CU", "MG2", "ZN2", "CA2", "NA", "CL"}
protein_idx = list(last_frame.topology.select("protein"))
print(f"[INFO] Protein atoms from MDTraj select: {len(protein_idx)}")

# Get ligand atom indices first so we can exclude them
ligand_idx = [a.index for a in last_frame.topology.atoms if a.residue.index == ligand_res_idx]
ligand_idx_set = set(ligand_idx)
print(f"[INFO] Ligand atoms to exclude: {len(ligand_idx)}")

# CRITICAL: Remove ligand atoms from protein selection (MDTraj may include non-standard residues)
protein_idx = [idx for idx in protein_idx if idx not in ligand_idx_set]
print(f"[INFO] Protein atoms after ligand removal: {len(protein_idx)}")

ion_idx = []
for r in last_frame.topology.residues:
    if r.name.upper() in STRUCTURAL_IONS:
        # Check if ion is close to protein (within 5A of any protein atom)
        ion_atoms = [a.index for a in r.atoms]
        ion_coords = last_frame.xyz[0, ion_atoms, :] * 10  # nm to A
        prot_coords = last_frame.xyz[0, protein_idx[:1000], :] * 10  # Sample protein atoms
        min_dist = np.min([np.linalg.norm(ion_coords[0] - pc) for pc in prot_coords])
        if min_dist < 5.0:
            ion_idx.extend(ion_atoms)
            print(f"[INFO] Including structural ion {r.name} (dist to protein: {min_dist:.1f}A)")

protein_with_ions_idx = sorted(protein_idx + ion_idx)

# Verify no ligand atoms in the final list
overlap = set(protein_with_ions_idx) & ligand_idx_set
if overlap:
    print(f"[FAIL] Ligand atoms still in protein selection: {len(overlap)} atoms")
    sys.exit(1)

# FIXED: For complex leg, save PROTEIN + IONS PDB (perses adds ligand from SDF)
protein_only_for_complex = last_frame.atom_slice(protein_with_ions_idx)
protein_only_for_complex.save("fep_study/inputs/WT_protein_for_complex.pdb")
print("[PASS] WT_protein_for_complex.pdb (", protein_only_for_complex.n_atoms, "atoms) - PROTEIN ONLY (no ligand)")

# FIXED: Post-slice verification - check BOTH resSeq and residue name match original
sliced_res702 = None
sliced_res702_candidates = []
for r in protein_only_for_complex.topology.residues:
    if r.resSeq == 702:
        sliced_res702_candidates.append(r)
    # Also check by original residue index in case of renumbering
    if r.name.upper() == target_res.name.upper():
        # Could be a match if nearby in sequence
        pass

# Verify we found the right residue
if not sliced_res702_candidates:
    print("[FAIL] No residue with resSeq=702 in sliced protein! MDTraj may have renumbered.")
    print("  Original residue:", target_res.name, "at resSeq", target_res.resSeq)
    # Show what residues we do have
    all_res = list(protein_only_for_complex.topology.residues)
    if all_res:
        print("  Sliced protein resSeq range:", min(r.resSeq for r in all_res), "-", max(r.resSeq for r in all_res))
    sys.exit(1)

# Check the residue is still ARG
sliced_res702 = sliced_res702_candidates[0]
if sliced_res702.name.upper() not in ["ARG","ARN","AR0"]:
    print("[FAIL] Post-slice: Residue 702 is", sliced_res702.name, "not ARG!")
    print("  This suggests MDTraj renumbered residues incorrectly.")
    sys.exit(1)

print("[PASS] Post-slice: Residue 702 verified as", sliced_res702.name)

# Also save protein+ions for protein leg
protein_only = last_frame.atom_slice(protein_with_ions_idx)
protein_only.save("fep_study/inputs/WT_protein_stripped.pdb")
print("[PASS] WT_protein_stripped.pdb (", protein_only.n_atoms, "atoms, includes structural ions)")

# Save ligand coordinates for reference (not used by perses, just for checking)
ligand_only = last_frame.atom_slice(ligand_idx)
ligand_only.save("fep_study/inputs/ligand_from_md.pdb")
print("[PASS] ligand_from_md.pdb (", ligand_only.n_atoms, "atoms) - reference only")

# STEP 1: Translate SDF to MD pose centroid
print()
print("[INFO] Translating SDF to MD-bound pose...")
aligned_sdf_path = "fep_study/inputs/febuxostat.sdf"
try:
    from rdkit.Chem import AllChem
    
    # Load SDF
    sdf_mol = Chem.SDMolSupplier(best_sdf)[0]
    if sdf_mol is None:
        raise ValueError("Could not read SDF")
    
    # Load MD ligand to get target centroid
    md_lig = md.load("fep_study/inputs/ligand_from_md.pdb")
    md_coords = md_lig.xyz[0] * 10  # nm to Angstrom
    md_centroid = np.mean(md_coords, axis=0)
    
    # Get SDF centroid
    conf = sdf_mol.GetConformer()
    sdf_coords = []
    for i in range(sdf_mol.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        sdf_coords.append([pos.x, pos.y, pos.z])
    sdf_coords = np.array(sdf_coords)
    sdf_centroid = np.mean(sdf_coords, axis=0)
    
    # Calculate translation vector
    translation = md_centroid - sdf_centroid
    print(f"  Translation vector: ({translation[0]:.1f}, {translation[1]:.1f}, {translation[2]:.1f}) A")
    
    # Apply translation to all atoms
    for i in range(sdf_mol.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, (pos.x + translation[0], 
                                  pos.y + translation[1], 
                                  pos.z + translation[2]))
    
    # Verify new centroid
    new_coords = []
    for i in range(sdf_mol.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        new_coords.append([pos.x, pos.y, pos.z])
    new_centroid = np.mean(new_coords, axis=0)
    print(f"  New SDF centroid: ({new_centroid[0]:.1f}, {new_centroid[1]:.1f}, {new_centroid[2]:.1f})")
    print(f"  MD centroid:      ({md_centroid[0]:.1f}, {md_centroid[1]:.1f}, {md_centroid[2]:.1f})")
    
    # Save translated SDF
    writer = Chem.SDWriter(aligned_sdf_path)
    writer.write(sdf_mol)
    writer.close()
    print("  [PASS] Saved translated SDF")
except Exception as e:
    print(f"  [WARNING] SDF translation error: {e}")
    print("  Using original SDF (perses will handle positioning)")
    shutil.copy(best_sdf, aligned_sdf_path)

# STEP 2: Validate the ALIGNED SDF
# CRITICAL: Validate ligand pose - SDF must match MD-bound position
print()
print("[INFO] Validating ligand pose (SDF vs MD-bound)...")
def validate_ligand_pose(sdf_path, md_ligand_pdb):
    """Check that SDF coordinates match MD-bound pose."""
    try:
        # Load SDF
        supplier = Chem.SDMolSupplier(sdf_path)
        sdf_mol = None
        for mol in supplier:
            if mol:
                sdf_mol = mol
                break
        if not sdf_mol:
            print("  [WARNING] Could not read SDF molecule")
            return False, "SDF unreadable"
        
        # Check if SDF has valid 3D coordinates
        conf = sdf_mol.GetConformer()
        sdf_coords = []
        z_coords = []
        for i in range(sdf_mol.GetNumAtoms()):
            pos = conf.GetAtomPosition(i)
            sdf_coords.append((pos.x, pos.y, pos.z))
            z_coords.append(pos.z)
        
        # Check for 2D coordinates (all Z near zero)
        z_range = max(z_coords) - min(z_coords)
        if z_range < 0.1:
            print(f"  [FAIL] SDF appears to be 2D (Z-range={z_range:.3f}A)")
            return False, "SDF is 2D, not 3D"
        
        # Check for all-zero coordinates
        coord_range = max(max(c) for c in sdf_coords) - min(min(c) for c in sdf_coords)
        if coord_range < 1.0:
            print(f"  [FAIL] SDF has invalid coordinates (range={coord_range:.3f}A)")
            return False, "SDF coordinates invalid"
        
        # Load MD ligand
        md_lig = md.load(md_ligand_pdb)
        md_coords = md_lig.xyz[0] * 10  # nm to Angstrom
        
        # Compare heavy atom counts
        sdf_heavy = sdf_mol.GetNumHeavyAtoms()
        md_heavy = 0
        for a in md_lig.topology.atoms:
            if a.element is not None and a.element.symbol != "H":
                md_heavy += 1
            elif a.element is None and not a.name.startswith("H"):
                md_heavy += 1  # Fallback for missing element
        
        if sdf_heavy != md_heavy:
            print(f"  [FAIL] Heavy atom mismatch: SDF={sdf_heavy}, MD={md_heavy}")
            print("  This indicates wrong ligand or different protonation state.")
            return False, f"Heavy atom mismatch: SDF={sdf_heavy} vs MD={md_heavy}"
        
        # Calculate centroid distance
        sdf_centroid = np.mean(sdf_coords, axis=0)
        md_centroid = np.mean(md_coords, axis=0)
        centroid_dist = np.linalg.norm(sdf_centroid - md_centroid)
        
        print(f"  SDF centroid: ({sdf_centroid[0]:.1f}, {sdf_centroid[1]:.1f}, {sdf_centroid[2]:.1f})")
        print(f"  MD centroid:  ({md_centroid[0]:.1f}, {md_centroid[1]:.1f}, {md_centroid[2]:.1f})")
        print(f"  Centroid distance: {centroid_dist:.2f} A")
        
        # Also calculate coordinate spread to detect rotated poses
        sdf_spread = np.std(sdf_coords, axis=0)
        md_spread = np.std(md_coords, axis=0)
        spread_diff = np.linalg.norm(sdf_spread - md_spread)
        print(f"  Coordinate spread difference: {spread_diff:.2f} A")
        
        if centroid_dist > 5.0:
            print("  [FAIL] Ligand poses differ by >5A! SDF may not be in bound pose.")
            print("  The SDF coordinates must match the MD simulation pose.")
            return False, f"Centroid distance {centroid_dist:.1f}A > 5A"
        elif spread_diff > 3.0:
            print("  [FAIL] Ligand appears rotated/different conformation (spread diff >3A).")
            return False, f"Conformation mismatch (spread diff {spread_diff:.1f}A)"
        elif centroid_dist > 2.0:
            print("  [WARNING] Ligand poses differ by >2A. Verify SDF is correct pose.")
            return True, f"Warning: {centroid_dist:.1f}A difference"
        else:
            print("  [PASS] Ligand poses match (<2A centroid, <3A spread)")
            return True, "OK"
            
    except Exception as e:
        print(f"  [FAIL] Pose validation error: {e}")
        return False, f"Validation error: {e}"

pose_ok, pose_msg = validate_ligand_pose(aligned_sdf_path, "fep_study/inputs/ligand_from_md.pdb")
if not pose_ok:
    print("[FAIL] Ligand pose validation failed:", pose_msg)
    print("  Ensure the SDF contains the docked/bound pose, not a 2D or random conformation.")
    sys.exit(1)

print("[PASS] Using aligned febuxostat.sdf")

print()
print("="*60, flush=True)
print("PART 3A: COMPLEX HYBRID", flush=True)
print("="*60, flush=True)
print("[INFO] This step can take 5-15 minutes for large systems...", flush=True)
print("[INFO] Perses will: parameterize ligand, solvate system, create hybrid topology", flush=True)

# Charge verification function
def verify_system_charge(system, name):
    """Verify system is approximately neutral (for PME)."""
    from openmm import NonbondedForce
    total_charge = 0.0
    nb_count = 0
    counted_particles = set()
    
    for i in range(system.getNumForces()):
        force = system.getForce(i)
        if isinstance(force, NonbondedForce):
            nb_count += 1
            for j in range(force.getNumParticles()):
                if j not in counted_particles:
                    charge, sigma, epsilon = force.getParticleParameters(j)
                    total_charge += charge.value_in_unit(unit.elementary_charge)
                    counted_particles.add(j)
    
    if nb_count == 0:
        print(f"[WARNING] {name}: No NonbondedForce found! Cannot verify charge.")
        return None
    if nb_count > 1:
        print(f"[INFO] {name}: Found {nb_count} NonbondedForces (counted particles once)")
    
    print(f"  [{name}] Total charge: {total_charge:.2f} e")
    if abs(total_charge) > 1.5:
        print(f"[WARNING] System {name} has significant net charge! May cause PME artifacts.")
    else:
        print(f"[PASS] {name} is approximately neutral")
    return total_charge

# FIXED: Explicit forcefield specification for reproducibility
# Using amber14 to match typical MD setups
FORCEFIELD_FILES = ['amber14/protein.ff14SB.xml', 'amber14/tip3p.xml']
print("[INFO] Forcefields:", FORCEFIELD_FILES)

def try_create_executor(pdb_file, chain_opts, ligand_sdf=None):
    """Create PointMutationExecutor with multiple chain fallbacks."""
    last_error = None
    for cid in chain_opts:
        chain_error = None  # Reset per chain
        # Try multiple configurations: OpenFF -> GAFF, with ligand_input -> ligand_file fallback
        configs = [
            {"use_openff": True, "use_ff_files": True, "ligand_param": "ligand_input"},
            {"use_openff": True, "use_ff_files": False, "ligand_param": "ligand_input"},
            {"use_openff": False, "use_ff_files": True, "ligand_param": "ligand_input"},
            {"use_openff": False, "use_ff_files": False, "ligand_param": "ligand_input"},
            # Fallback to ligand_file for older perses versions
            {"use_openff": True, "use_ff_files": True, "ligand_param": "ligand_file"},
            {"use_openff": False, "use_ff_files": True, "ligand_param": "ligand_file"},
            {"use_openff": False, "use_ff_files": False, "ligand_param": "ligand_file"},
        ]
        for cfg_idx, cfg in enumerate(configs):
            try:
                ff_info = "OpenFF" if cfg["use_openff"] else "GAFF"
                lig_param = cfg["ligand_param"]
                print(f"[DEBUG] Trying config {cfg_idx+1}/{len(configs)}: chain={cid}, {ff_info}, {lig_param}", flush=True)

                kwargs = {
                    "protein_filename": pdb_file,
                    "mutation_chain_id": cid,
                    "mutation_residue_id": "702",
                    "proposed_residue": "TRP",
                    "conduct_endstate_validation": False,
                    "ionic_strength": 0.15 * unit.molar,
                }
                if cfg["use_openff"]:
                    kwargs["small_molecule_forcefields"] = "openff-2.1.0"
                if cfg["use_ff_files"]:
                    kwargs["forcefield_files"] = FORCEFIELD_FILES
                if ligand_sdf:
                    kwargs[cfg["ligand_param"]] = ligand_sdf

                print(f"[DEBUG] Calling PointMutationExecutor... (this may take several minutes)", flush=True)
                start_time = time.time()
                executor = PointMutationExecutor(**kwargs)
                elapsed = time.time() - start_time
                print(f"[DEBUG] PointMutationExecutor completed in {elapsed:.1f}s", flush=True)
                ff_info = "OpenFF" if cfg["use_openff"] else "GAFF"
                lig_param = cfg["ligand_param"]
                print(f"[PASS] Created with chain={cid} ({ff_info}, {lig_param})")
                return executor
            except TypeError as te:
                # Parameter not accepted, try next config
                print(f"[DEBUG] TypeError (param not accepted): {str(te)[:60]}", flush=True)
                chain_error = chain_error or te  # Keep first error
                continue
            except Exception as exc:
                print(f"[DEBUG] Exception: {type(exc).__name__}: {str(exc)[:80]}", flush=True)
                chain_error = exc
                # Continue trying other configs (don't break!)
                continue
        
        # All configs failed for this chain
        if chain_error:
            last_error = chain_error
            print(f"[INFO] chain={cid} failed: {str(chain_error)[:80]}")

    raise RuntimeError("All chain options failed. Last: " + str(last_error))

def try_create_protein_executor(pdb_file, chain_opts):
    """Create PointMutationExecutor for protein-only leg."""
    last_error = None
    for cid in chain_opts:
        try:
            kwargs = {
                "protein_filename": pdb_file,
                "mutation_chain_id": cid,
                "mutation_residue_id": "702",
                "proposed_residue": "TRP",
                "conduct_endstate_validation": False,
                "ionic_strength": 0.15 * unit.molar,
            }

            # Try with forcefield specification
            try:
                kwargs["forcefield_files"] = FORCEFIELD_FILES
                executor = PointMutationExecutor(**kwargs)
            except TypeError:
                del kwargs["forcefield_files"]
                executor = PointMutationExecutor(**kwargs)

            print("[PASS] Created with chain=", cid, "(protein-only)")
            return executor
        except Exception as e:
            last_error = e
            print("[INFO] chain=", cid, "failed:", str(e)[:80])
            continue

    raise RuntimeError("All options failed. Last: " + str(last_error))

# FIXED: Use protein-only PDB for complex (perses adds ligand from SDF)
print("[DEBUG] Starting complex executor creation...", flush=True)
import signal

def timeout_handler(signum, frame):
    print("\n[FAIL] Timeout - PointMutationExecutor took too long (>10 min)", flush=True)
    sys.exit(1)

# Set timeout for PART 3A (10 minutes)
try:
    signal.signal(signal.SIGALRM, timeout_handler)
    signal.alarm(600)  # 10 minute timeout
except AttributeError:
    print("[WARNING] Signal not available on this platform, no timeout", flush=True)

try:
    complex_executor = try_create_executor(
        "fep_study/inputs/WT_protein_for_complex.pdb",
        chain_options,
        ligand_sdf="fep_study/inputs/febuxostat.sdf"
    )
    # Cancel the alarm
    try:
        signal.alarm(0)
    except:
        pass
except KeyboardInterrupt:
    print("\n[FAIL] Interrupted by user", flush=True)
    sys.exit(1)
except MemoryError:
    print("[FAIL] Out of memory! System too large.", flush=True)
    sys.exit(1)
except Exception as e:
    print(f"\n[FAIL] Exception type: {type(e).__name__}", flush=True)
    print(f"[FAIL] Message: {e}", flush=True)
    import traceback
    traceback.print_exc()
    sys.exit(1)

# Save hybrid system components (NOT the executor - it contains unpicklable SwigPyObjects)
with open("fep_study/complex/hybrid_system.xml", "w") as f:
    f.write(XmlSerializer.serialize(complex_executor.hybrid_system))
with open("fep_study/complex/hybrid_topology.pkl", "wb") as f:
    pickle.dump(complex_executor.hybrid_topology, f)
np.save("fep_study/complex/hybrid_positions.npy", complex_executor.hybrid_positions)
# FIXED: Do NOT pickle executor (contains SwigPyObjects that can't be pickled)
print("[PASS] Complex saved (system, topology, positions)")

# Verify charge
verify_system_charge(complex_executor.hybrid_system, "COMPLEX")

print()
print("="*60)
print("PART 3B: PROTEIN HYBRID")
print("="*60)

try:
    protein_executor = try_create_protein_executor("fep_study/inputs/WT_protein_stripped.pdb",
                                                   chain_options)
except Exception as e:
    print("[FAIL]", e)
    import traceback
    traceback.print_exc()
    sys.exit(1)

with open("fep_study/protein/hybrid_system.xml", "w") as f:
    f.write(XmlSerializer.serialize(protein_executor.hybrid_system))
with open("fep_study/protein/hybrid_topology.pkl", "wb") as f:
    pickle.dump(protein_executor.hybrid_topology, f)
np.save("fep_study/protein/hybrid_positions.npy", protein_executor.hybrid_positions)
# FIXED: Do NOT pickle executor
print("[PASS] Protein saved (system, topology, positions)")

# Verify charge
verify_system_charge(protein_executor.hybrid_system, "PROTEIN")

def save_preview_pdb(executor, positions, path):
    try:
        topo = executor.hybrid_topology
        if hasattr(topo, "to_openmm"):
            topo = topo.to_openmm()
        if not isinstance(topo, Topology):
            print("[WARNING] Topo is", type(topo), ", skip preview")
            return False
        with open(path, "w") as f:
            PDBFile.writeFile(topo, positions, f)
        print("[PASS] Preview:", os.path.basename(path))
        return True
    except Exception as e:
        print("[WARNING] Preview failed:", e)
        return False

def run_stability_test(executor, name, output_dir):
    """Run stability test at lambda=0.0, 0.5, and 1.0 (endpoints + midpoint)."""
    print()
    print("--- Stability:", name, "---")

    topo_for_sim = executor.hybrid_topology
    if hasattr(topo_for_sim, "to_openmm"):
        topo_for_sim = topo_for_sim.to_openmm()

    integrator = LangevinMiddleIntegrator(310*unit.kelvin, 1/unit.picosecond, 0.002*unit.picoseconds)

    try:
        platform = Platform.getPlatformByName("CUDA")
        print("[INFO] CUDA")
    except:
        platform = Platform.getPlatformByName("CPU")
        print("[INFO] CPU")

    sim = Simulation(topo_for_sim, executor.hybrid_system, integrator, platform)
    sim.context.setPositions(executor.hybrid_positions)

    params = list(sim.context.getParameters())
    lam_params = [p for p in params if p.startswith("lambda_") or p == "lambda"]
    if not lam_params:
        lam_params = [p for p in params if "lambda" in p.lower()]

    if not lam_params:
        print("[FAIL] No lambda params!")
        return False

    print("  Lambda params found:", len(lam_params))

    # FIXED: Test all three critical lambda values (endpoints + midpoint)
    for lam_val in [0.0, 0.5, 1.0]:
        print()
        print(f"  [....] Testing lambda={lam_val}...")

        for p in lam_params:
            sim.context.setParameter(p, lam_val)

        # Minimize
        try:
            sim.minimizeEnergy(maxIterations=500)
            state = sim.context.getState(getEnergy=True, getPositions=True)
            energy = state.getPotentialEnergy()
            print(f"    Minimized energy (lambda={lam_val}):", energy)

            if "nan" in str(energy).lower():
                print(f"[FAIL] NaN at lambda={lam_val}!")
                return False

            # Check for unreasonable energies
            energy_val = energy.value_in_unit(unit.kilojoules_per_mole)
            if abs(energy_val) > 1e8:
                print(f"[FAIL] Unreasonable energy at lambda={lam_val}: {energy_val} kJ/mol")
                return False

        except Exception as e:
            print(f"[FAIL] Minimization failed at lambda={lam_val}:", e)
            return False

    # Run short dynamics at lambda=0.5
    print()
    print("  [....] 10ps dynamics at lambda=0.5...")
    for p in lam_params:
        sim.context.setParameter(p, 0.5)

    sim.minimizeEnergy(maxIterations=500)
    sim.step(5000)
    state = sim.context.getState(getEnergy=True, getPositions=True)
    energy = state.getPotentialEnergy()
    print("    Final energy:", energy)

    if "nan" in str(energy).lower():
        print("[FAIL] NaN after dynamics!")
        return False

    save_preview_pdb(executor, state.getPositions(), output_dir + "/preview_final.pdb")

    print("[PASS]", name, "- all lambda values stable")
    return True

print()
print("="*60)
print("PART 4: STABILITY TESTS")
print("="*60)

if not run_stability_test(complex_executor, "COMPLEX", "fep_study/complex"):
    sys.exit(1)

if not run_stability_test(protein_executor, "PROTEIN", "fep_study/protein"):
    sys.exit(1)

print()
print("="*60)
print("COMPLETE")
print("="*60)
print()
print("FEP setup complete!")
print()
print("Output files:")
print("  fep_study/complex/hybrid_system.xml")
print("  fep_study/complex/hybrid_topology.pkl")
print("  fep_study/complex/hybrid_positions.npy")
print("  fep_study/protein/hybrid_system.xml")
print("  fep_study/protein/hybrid_topology.pkl")
print("  fep_study/protein/hybrid_positions.npy")
print()
print("Next: Run FEP sampling with run_fep_sampling.py")
print("="*60)
