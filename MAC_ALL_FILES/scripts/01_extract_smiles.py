#!/usr/bin/env python3
"""
NOD2-Scout: Extract SMILES from e-Drug3D SDF and map to screening results
"""

import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit import RDLogger
from pathlib import Path

RDLogger.DisableLog('rdApp.*')

print("=" * 70)
print("NOD2-SCOUT: EXTRACTING SMILES FROM e-Drug3D")
print("=" * 70)

# Paths
BASE_DIR = Path(__file__).resolve().parents[2]
DATA_DIR = BASE_DIR / 'nod2-screening-data'
OUTPUT_DIR = BASE_DIR / 'nod2-scout' / 'data'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Load SDF file
print("\n1. Loading e-Drug3D SDF file...")
suppl = Chem.SDMolSupplier(str(DATA_DIR / 'e-Drug3D_2162.sdf'))

# Extract SMILES and create ID mapping
sdf_data = {}
for mol in suppl:
    if mol is None:
        continue

    try:
        drug_id = int(mol.GetProp('ID'))
        name = mol.GetProp('name').strip() if mol.HasProp('name') else None
        smiles = Chem.MolToSmiles(mol)
        sdf_data[drug_id] = {'name': name, 'smiles': smiles}
    except:
        continue

print(f"   Extracted {len(sdf_data)} FDA drugs from SDF")

# Load screening results
print("\n2. Loading screening results...")
screening = pd.read_csv(DATA_DIR / 'screening_results.csv')
print(f"   Total compounds: {len(screening)}")

# Separate FDA and Natural
fda_df = screening[screening['type'] == 'FDA'].copy()
nat_df = screening[screening['type'] == 'Natural'].copy()

print(f"   FDA drugs: {len(fda_df)}")
print(f"   Natural products: {len(nat_df)}")

# Map FDA drugs to SMILES
print("\n3. Mapping FDA drugs to SMILES...")
fda_df['drug_id'] = fda_df['compound_id'].str.extract(r'fda_(\d+)')[0].astype(int)

def get_fda_smiles(drug_id):
    if drug_id in sdf_data:
        return sdf_data[drug_id]['smiles']
    return None

def get_fda_name(drug_id):
    if drug_id in sdf_data:
        return sdf_data[drug_id]['name']
    return None

fda_df['smiles_new'] = fda_df['drug_id'].apply(get_fda_smiles)
fda_df['drug_name'] = fda_df['drug_id'].apply(get_fda_name)

fda_matched = fda_df['smiles_new'].notna().sum()
print(f"   FDA drugs with SMILES: {fda_matched}/{len(fda_df)}")

# For Natural products, the compound_id might be PubChem CIDs
# We'll try to fetch these later, for now mark as needing fetch
print("\n4. Processing Natural products...")
nat_df['smiles_new'] = None
nat_df['drug_name'] = None

# Check if any natural products already have valid SMILES in the original column
from rdkit import Chem
valid_nat_smiles = 0
for idx, row in nat_df.iterrows():
    if pd.notna(row['smiles']):
        mol = Chem.MolFromSmiles(str(row['smiles']))
        if mol is not None:
            nat_df.at[idx, 'smiles_new'] = str(row['smiles'])
            valid_nat_smiles += 1

print(f"   Natural products with existing valid SMILES: {valid_nat_smiles}")

# Combine datasets
print("\n5. Combining datasets...")
fda_df['final_smiles'] = fda_df['smiles_new']
nat_df['final_smiles'] = nat_df['smiles_new']

combined = pd.concat([fda_df, nat_df], ignore_index=True)

# Keep only compounds with SMILES for ML
valid_df = combined[combined['final_smiles'].notna()].copy()
valid_df = valid_df[['compound_id', 'final_smiles', 'cnn_affinity', 'type', 'drug_name']].copy()
valid_df.columns = ['compound_id', 'smiles', 'cnn_affinity', 'type', 'drug_name']

print(f"\n   Total compounds with SMILES: {len(valid_df)}")
print(f"   - FDA: {len(valid_df[valid_df['type']=='FDA'])}")
print(f"   - Natural: {len(valid_df[valid_df['type']=='Natural'])}")

# Save
print("\n6. Saving results...")
valid_df.to_csv(OUTPUT_DIR / 'compounds_with_smiles.csv', index=False)
print(f"   Saved: {OUTPUT_DIR / 'compounds_with_smiles.csv'}")

# Summary statistics
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)

print(f"\nCNN Affinity by Type:")
print(valid_df.groupby('type')['cnn_affinity'].describe().round(2))

print(f"\nTop 20 FDA Drugs by CNN Affinity:")
top_fda = valid_df[valid_df['type'] == 'FDA'].nlargest(20, 'cnn_affinity')
for i, row in top_fda.iterrows():
    print(f"  {row['cnn_affinity']:.3f}: {row['drug_name']}")

# Verify SMILES are valid
print(f"\nValidating SMILES...")
invalid_count = 0
for idx, row in valid_df.iterrows():
    mol = Chem.MolFromSmiles(row['smiles'])
    if mol is None:
        invalid_count += 1

print(f"  Invalid SMILES: {invalid_count}/{len(valid_df)}")
print(f"  Valid SMILES: {len(valid_df) - invalid_count}/{len(valid_df)}")
