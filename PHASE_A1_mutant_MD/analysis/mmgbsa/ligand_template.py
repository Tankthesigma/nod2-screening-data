#!/usr/bin/env python
"""
Create OpenMM residue template for ligand from SDF file
Uses RDKit for molecular structure and generates AMBER-like parameters
"""

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
import numpy as np
import os

# GAFF-like atom type parameters
# Format: (sigma in nm, epsilon in kJ/mol, mass in amu)
GAFF_TYPES = {
    # Carbon types
    'c': (0.3400, 0.3598, 12.01),
    'c1': (0.3400, 0.3598, 12.01),
    'c2': (0.3400, 0.3598, 12.01),
    'c3': (0.3400, 0.4577, 12.01),
    'ca': (0.3400, 0.3598, 12.01),

    # Nitrogen types
    'n': (0.3250, 0.7113, 14.01),
    'n1': (0.3250, 0.7113, 14.01),
    'n2': (0.3250, 0.7113, 14.01),
    'n3': (0.3250, 0.7113, 14.01),
    'nb': (0.3250, 0.7113, 14.01),

    # Oxygen types
    'o': (0.2960, 0.8786, 16.00),
    'oh': (0.3066, 0.8803, 16.00),
    'os': (0.3000, 0.7113, 16.00),

    # Sulfur types
    's': (0.3560, 1.0460, 32.06),
    's4': (0.3560, 1.0460, 32.06),
    's6': (0.3560, 1.0460, 32.06),
    'ss': (0.3560, 1.0460, 32.06),

    # Hydrogen types
    'ha': (0.2599, 0.0628, 1.008),
    'hc': (0.2649, 0.0657, 1.008),
    'hn': (0.1069, 0.0657, 1.008),
    'ho': (0.0001, 0.0000, 1.008),

    # Halogens
    'f': (0.3118, 0.2552, 19.00),
    'cl': (0.3470, 1.1088, 35.45),
    'br': (0.3910, 1.7023, 79.90),
    'i': (0.4190, 1.6736, 126.9),

    # Phosphorus
    'p5': (0.3742, 0.8368, 30.97),
}


def assign_gaff_type(atom, mol):
    """Assign GAFF-like atom type based on atom properties"""
    elem = atom.GetSymbol().lower()
    hyb = atom.GetHybridization()
    aromatic = atom.GetIsAromatic()

    # Hydrogens
    if elem == 'h':
        neighbors = list(atom.GetNeighbors())
        if len(neighbors) > 0:
            neighbor = neighbors[0]
            n_elem = neighbor.GetSymbol().lower()
            if n_elem == 'o':
                return 'ho'
            elif n_elem == 'n':
                return 'hn'
            elif neighbor.GetIsAromatic():
                return 'ha'
        return 'hc'

    # Carbons
    if elem == 'c':
        if aromatic:
            return 'ca'
        if hyb == Chem.HybridizationType.SP3:
            return 'c3'
        if hyb == Chem.HybridizationType.SP2:
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'O':
                    bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                    if bond and bond.GetBondType() == Chem.BondType.DOUBLE:
                        return 'c'
            return 'c2'
        return 'c3'

    # Nitrogens
    if elem == 'n':
        if aromatic:
            return 'nb'
        if hyb == Chem.HybridizationType.SP3:
            return 'n3'
        if hyb == Chem.HybridizationType.SP2:
            return 'n2'
        return 'n3'

    # Oxygens
    if elem == 'o':
        for neighbor in atom.GetNeighbors():
            if neighbor.GetSymbol() == 'H':
                return 'oh'
        if hyb == Chem.HybridizationType.SP3:
            return 'os'
        return 'o'

    # Sulfurs
    if elem == 's':
        o_neighbors = [n for n in atom.GetNeighbors() if n.GetSymbol() == 'O']
        o_count = len(o_neighbors)
        if o_count >= 2:
            return 's6'
        if o_count == 1:
            return 's4'
        return 'ss'

    # Halogens
    if elem in ['f', 'cl', 'br', 'i']:
        return elem

    # Phosphorus
    if elem == 'p':
        return 'p5'

    return 'c3'


def estimate_partial_charges(mol):
    """Estimate partial charges using Gasteiger method"""
    AllChem.ComputeGasteigerCharges(mol)
    charges = []
    for atom in mol.GetAtoms():
        charge = atom.GetDoubleProp('_GasteigerCharge')
        if np.isnan(charge):
            charge = 0.0
        charges.append(charge)
    return charges


def create_ligand_template(sdf_path, resname='UNK'):
    """
    Create ligand template info from SDF file
    """
    mol = Chem.SDMolSupplier(sdf_path, removeHs=False)[0]
    if mol is None:
        raise ValueError(f"Could not load molecule from {sdf_path}")

    # Add hydrogens if needed
    mol = Chem.AddHs(mol, addCoords=True)

    # Compute partial charges
    charges = estimate_partial_charges(mol)

    # Assign atom types
    atoms = []
    for i, atom in enumerate(mol.GetAtoms()):
        gaff_type = assign_gaff_type(atom, mol)
        params = GAFF_TYPES.get(gaff_type, GAFF_TYPES['c3'])

        atoms.append({
            'name': f"{atom.GetSymbol()}{i+1}",
            'element': atom.GetSymbol(),
            'gaff_type': gaff_type,
            'charge': charges[i],
            'sigma': params[0],
            'epsilon': params[1],
            'mass': params[2]
        })

    # Get bonds
    bonds = []
    for bond in mol.GetBonds():
        bonds.append((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))

    return {
        'resname': resname,
        'atoms': atoms,
        'bonds': bonds,
        'n_atoms': len(atoms),
        'total_charge': sum(charges)
    }


if __name__ == '__main__':
    base = r'C:\Users\vasud\nod2-screening-data'
    sdf = os.path.join(base, r'PHASE_6\structures\febuxostat_docked.sdf')

    print("="*60)
    print("Ligand Template Generation Test")
    print("="*60)

    template = create_ligand_template(sdf)

    print(f"\nFebuxostat template:")
    print(f"  Atoms: {template['n_atoms']}")
    print(f"  Net charge: {template['total_charge']:.3f}")
    print(f"\nAtom details:")
    for i, atom in enumerate(template['atoms'][:15]):
        print(f"  {i+1:2d}. {atom['name']:4s} ({atom['gaff_type']:3s}) q={atom['charge']:+.3f} sigma={atom['sigma']:.4f}")
    if template['n_atoms'] > 15:
        print(f"  ... and {template['n_atoms'] - 15} more atoms")

    # Test other ligands
    print("\n" + "="*60)
    print("Testing all ligands:")
    print("="*60)

    ligands = [
        ('febuxostat_docked.sdf', 'Febuxostat'),
        ('natural_top_docked.sdf', 'Natural'),
        ('budesonide_docked.sdf', 'Budesonide'),
        ('ursodiol_docked.sdf', 'Ursodiol'),
        ('decoy_docked.sdf', 'Decoy'),
    ]

    for sdf_file, name in ligands:
        sdf_path = os.path.join(base, r'PHASE_6\structures', sdf_file)
        template = create_ligand_template(sdf_path)
        print(f"  {name:12s}: {template['n_atoms']:3d} atoms, charge={template['total_charge']:+.3f}")
