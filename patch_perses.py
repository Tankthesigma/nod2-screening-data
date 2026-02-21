"""
Patch perses to use RDKit/OpenFF instead of OpenEye.

This module must be imported BEFORE importing perses.
"""
import sys
import numpy as np
from simtk import unit

# Create mock oechem module to prevent import errors
class MockOEMol:
    """Placeholder class - should never be instantiated."""
    pass

class MockOEChem:
    """Mock oechem module - returns False for all license checks."""
    OEMol = MockOEMol
    OEGraphMol = MockOEMol

    @staticmethod
    def OEChemIsLicensed():
        return False

    @staticmethod
    def OEIsLicensed(*args, **kwargs):
        return False

class MockOEQuacpac:
    """Mock oequacpac module."""
    @staticmethod
    def OEQuacPacIsLicensed():
        return False

class MockOEOmega:
    """Mock oeomega module."""
    @staticmethod
    def OEOmegaIsLicensed():
        return False

class MockOEIUPAC:
    """Mock oeiupac module."""
    @staticmethod
    def OEIUPACIsLicensed():
        return False

class MockOEStructChecker:
    """Mock oestructchecker module."""
    @staticmethod
    def OEStructCheckerIsLicensed():
        return False

# Install mocks before any openeye imports
class MockOpenEye:
    oechem = MockOEChem()
    oequacpac = MockOEQuacpac()
    oeomega = MockOEOmega()
    oeiupac = MockOEIUPAC()
    oestructchecker = MockOEStructChecker()

sys.modules['openeye'] = MockOpenEye()
sys.modules['openeye.oechem'] = MockOEChem()
sys.modules['openeye.oequacpac'] = MockOEQuacpac()
sys.modules['openeye.oeomega'] = MockOEOmega()
sys.modules['openeye.oeiupac'] = MockOEIUPAC()
sys.modules['openeye.oestructchecker'] = MockOEStructChecker()

print("[PATCH] Installed mock OpenEye module")

# Now we can import RDKit and OpenFF safely
from rdkit import Chem
from openff.toolkit.topology import Molecule
from openmm import app

def createOEMolFromSDF_rdkit(filename, index=0, allow_undefined_stereo=True, add_hydrogens=False):
    """
    Create an RDKit molecule from an SDF file.

    This replaces the OpenEye-based function in perses.
    """
    # Handle PDB files as well (perses passes PDB files to this function for residue templates)
    if filename.endswith('.pdb'):
        mol = Chem.MolFromPDBFile(filename, removeHs=False)
        if mol is None:
            raise ValueError(f"Could not read PDB file: {filename}")
        if add_hydrogens:
            mol = Chem.AddHs(mol)
        return mol

    # Handle SDF files
    suppl = Chem.SDMolSupplier(filename, removeHs=False)
    mols = [m for m in suppl if m is not None]

    if not mols:
        raise ValueError(f"No molecules found in {filename}")

    if index >= len(mols):
        raise ValueError(f"Index {index} out of range (only {len(mols)} molecules)")

    mol = mols[index]
    if add_hydrogens:
        mol = Chem.AddHs(mol)

    return mol

def extractPositionsFromOEMol_rdkit(mol):
    """
    Extract positions from an RDKit molecule.

    This replaces the OpenEye-based function in perses.
    Returns positions as OpenMM Quantity with units of angstroms.
    """
    if not isinstance(mol, Chem.Mol):
        raise TypeError(f"Expected RDKit Mol, got {type(mol)}")

    conf = mol.GetConformer()
    positions = []
    for i in range(mol.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        positions.append([pos.x, pos.y, pos.z])

    return np.array(positions) * unit.angstrom

def generateTopologyFromOEMol_rdkit(mol):
    """
    Generate OpenMM topology from RDKit molecule.

    This replaces the OpenEye-based function.
    """
    # Use OpenFF to convert RDKit to topology
    off_mol = Molecule.from_rdkit(mol, allow_undefined_stereo=True)
    off_top = off_mol.to_topology()
    return off_top.to_openmm()

# Patch perses.utils.openeye functions
def patch_perses():
    """Apply patches to perses modules."""
    import perses.utils.openeye as oe_utils

    # Replace functions
    oe_utils.createOEMolFromSDF = createOEMolFromSDF_rdkit
    oe_utils.extractPositionsFromOEMol = extractPositionsFromOEMol_rdkit

    print("[PATCH] Replaced perses.utils.openeye functions with RDKit versions")

    # Also need to patch forcefield_generators if it's used
    try:
        from openmmforcefields import generators as fg
        fg.generateTopologyFromOEMol = generateTopologyFromOEMol_rdkit
        print("[PATCH] Replaced forcefield_generators.generateTopologyFromOEMol")
    except:
        pass

# Apply patches
try:
    patch_perses()
except ImportError as e:
    print(f"[PATCH WARNING] Could not apply patches: {e}")

print("[PATCH] Perses patching complete")
print()
