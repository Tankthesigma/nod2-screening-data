#!/usr/bin/env python3
"""
NOD2-Scout Step 2: Generate Molecular Features
- Morgan fingerprints (2048 bits, radius 2)
- Physicochemical descriptors
- Drug-likeness scores
- PAINS alerts
"""

import pandas as pd
import numpy as np
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, QED, Lipinski
from rdkit.Chem import rdMolDescriptors, FilterCatalog
from rdkit import RDLogger
import pickle

# Suppress RDKit warnings
RDLogger.DisableLog('rdApp.*')

print("=" * 70)
print("NOD2-SCOUT: FEATURE GENERATION")
print("=" * 70)

# Paths
BASE_DIR = Path(__file__).resolve().parents[2]
DATA_DIR = BASE_DIR / 'nod2-scout' / 'data'

# Load valid compounds
df = pd.read_csv(DATA_DIR / 'valid_compounds.csv')
print(f"Loaded {len(df)} compounds with SMILES")

# Initialize PAINS filter
pains_params = FilterCatalog.FilterCatalogParams()
pains_params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
pains_catalog = FilterCatalog.FilterCatalog(pains_params)

def compute_features(smiles):
    """Compute all features for a SMILES string"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        features = {}

        # Morgan fingerprint (2048 bits, radius 2)
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
        fp_array = np.zeros(2048, dtype=np.int8)
        for i in fp.GetOnBits():
            fp_array[i] = 1
        features['fingerprint'] = fp_array

        # Physicochemical descriptors
        features['mw'] = Descriptors.MolWt(mol)
        features['logp'] = Descriptors.MolLogP(mol)
        features['tpsa'] = Descriptors.TPSA(mol)
        features['hbd'] = Descriptors.NumHDonors(mol)
        features['hba'] = Descriptors.NumHAcceptors(mol)
        features['rotatable_bonds'] = Descriptors.NumRotatableBonds(mol)
        features['aromatic_rings'] = rdMolDescriptors.CalcNumAromaticRings(mol)
        features['fraction_csp3'] = rdMolDescriptors.CalcFractionCSP3(mol)
        features['num_heavy_atoms'] = mol.GetNumHeavyAtoms()
        features['num_rings'] = rdMolDescriptors.CalcNumRings(mol)

        # Drug-likeness
        features['qed'] = QED.qed(mol)

        # Lipinski violations
        lipinski_violations = 0
        if features['mw'] > 500: lipinski_violations += 1
        if features['logp'] > 5: lipinski_violations += 1
        if features['hbd'] > 5: lipinski_violations += 1
        if features['hba'] > 10: lipinski_violations += 1
        features['lipinski_violations'] = lipinski_violations

        # Veber violations (oral bioavailability)
        veber_violations = 0
        if features['rotatable_bonds'] > 10: veber_violations += 1
        if features['tpsa'] > 140: veber_violations += 1
        features['veber_violations'] = veber_violations

        # PAINS alerts
        pains_matches = pains_catalog.GetMatches(mol)
        features['pains_alerts'] = len(pains_matches)

        # Intracellular access predictors (simplified)
        # High TPSA = poor membrane permeability
        # P-gp efflux more likely with MW > 400 and multiple aromatic rings
        features['permeability_score'] = 1 - min(features['tpsa'] / 140, 1)
        features['pgp_risk'] = int(features['mw'] > 400 and features['aromatic_rings'] >= 2)

        return features

    except Exception as e:
        return None

# Compute features for all compounds
print("\nComputing features...")
all_features = []
all_fingerprints = []
failed = 0

for idx, row in df.iterrows():
    features = compute_features(row['smiles'])

    if features is None:
        failed += 1
        continue

    # Store fingerprint separately
    fp = features.pop('fingerprint')
    all_fingerprints.append(fp)

    # Add metadata
    features['compound_id'] = row['compound_id']
    features['smiles'] = row['smiles']
    features['cnn_affinity'] = row['cnn_affinity']
    features['type'] = row['type']
    features['drug_name'] = row['drug_name'] if pd.notna(row['drug_name']) else None

    all_features.append(features)

    if (idx + 1) % 500 == 0:
        print(f"  Processed {idx + 1}/{len(df)} compounds...")

print(f"\nSuccessfully processed: {len(all_features)} compounds")
print(f"Failed: {failed} compounds")

# Create DataFrames
features_df = pd.DataFrame(all_features)
fingerprints_array = np.array(all_fingerprints)

print(f"\nFeature matrix shape: {features_df.shape}")
print(f"Fingerprint matrix shape: {fingerprints_array.shape}")

# Save
features_df.to_csv(DATA_DIR / 'features.csv', index=False)
np.save(DATA_DIR / 'fingerprints.npy', fingerprints_array)

print(f"\nSaved: {DATA_DIR / 'features.csv'}")
print(f"Saved: {DATA_DIR / 'fingerprints.npy'}")

# Summary statistics
print("\n--- FEATURE STATISTICS ---")
print(f"\nPhysicochemical Properties:")
for col in ['mw', 'logp', 'tpsa', 'hbd', 'hba', 'qed']:
    print(f"  {col}: {features_df[col].mean():.2f} ± {features_df[col].std():.2f}")

print(f"\nDrug-likeness:")
print(f"  Lipinski compliant (0 violations): {(features_df['lipinski_violations'] == 0).sum()}")
print(f"  Veber compliant (0 violations): {(features_df['veber_violations'] == 0).sum()}")
print(f"  PAINS-free: {(features_df['pains_alerts'] == 0).sum()}")

print(f"\nCNN Affinity:")
print(f"  Range: {features_df['cnn_affinity'].min():.2f} - {features_df['cnn_affinity'].max():.2f}")
print(f"  Mean: {features_df['cnn_affinity'].mean():.2f}")
print(f"  90th percentile: {features_df['cnn_affinity'].quantile(0.9):.2f}")
print(f"  10th percentile: {features_df['cnn_affinity'].quantile(0.1):.2f}")

# Distribution by type
print(f"\nBy Type:")
print(features_df.groupby('type').agg({
    'cnn_affinity': ['mean', 'max'],
    'qed': 'mean',
    'lipinski_violations': lambda x: (x == 0).mean()
}).round(3))
