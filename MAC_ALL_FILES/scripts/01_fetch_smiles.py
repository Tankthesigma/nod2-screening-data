#!/usr/bin/env python3
"""
NOD2-Scout Step 1: Fetch SMILES for all compounds
- FDA drugs: Extract drug names, query PubChem
- Natural products: Use compound_id as PubChem CID
"""

import pandas as pd
import numpy as np
import urllib.request
import json
import time
from pathlib import Path

print("=" * 70)
print("NOD2-SCOUT: FETCHING SMILES")
print("=" * 70)

# Load data
df = pd.read_csv('../nod2-screening-data/screening_results.csv')
print(f"Loaded {len(df)} compounds")
print(f"  FDA drugs: {len(df[df['type']=='FDA'])}")
print(f"  Natural products: {len(df[df['type']=='Natural'])}")

# Separate FDA and Natural
fda_df = df[df['type'] == 'FDA'].copy()
nat_df = df[df['type'] == 'Natural'].copy()

def fetch_smiles_by_name(name, timeout=10):
    """Fetch SMILES from PubChem by drug name"""
    try:
        # Clean name
        name = str(name).strip().upper()
        if name.startswith('_C') or name == 'NAN' or len(name) < 3:
            return None

        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/CanonicalSMILES/JSON"
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req, timeout=timeout) as response:
            data = json.loads(response.read().decode())
            return data['PropertyTable']['Properties'][0]['CanonicalSMILES']
    except:
        return None

def fetch_smiles_by_cid(cid, timeout=10):
    """Fetch SMILES from PubChem by CID"""
    try:
        cid = int(cid)
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/CanonicalSMILES/JSON"
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req, timeout=timeout) as response:
            data = json.loads(response.read().decode())
            return data['PropertyTable']['Properties'][0]['CanonicalSMILES']
    except:
        return None

# Process FDA drugs
print("\n--- Processing FDA Drugs ---")
fda_smiles = []
fda_names = []

for idx, row in fda_df.iterrows():
    # The 'smiles' column contains either drug name or placeholder
    name_field = str(row['smiles'])

    # Skip placeholders
    if name_field.startswith('_c') or name_field == 'nan':
        fda_smiles.append(None)
        fda_names.append(None)
        continue

    fda_names.append(name_field)

    # Try to fetch
    smiles = fetch_smiles_by_name(name_field)
    fda_smiles.append(smiles)

    if len(fda_smiles) % 100 == 0:
        success = sum(1 for s in fda_smiles if s is not None)
        print(f"  Processed {len(fda_smiles)}/{len(fda_df)} FDA drugs ({success} with SMILES)")

    time.sleep(0.1)  # Rate limiting

fda_df['drug_name'] = fda_names
fda_df['smiles_fetched'] = fda_smiles

fda_success = sum(1 for s in fda_smiles if s is not None)
print(f"FDA SMILES fetched: {fda_success}/{len(fda_df)}")

# Process Natural Products (use compound_id as CID)
print("\n--- Processing Natural Products ---")
nat_smiles = []

# Batch fetch for efficiency
batch_size = 100
nat_ids = nat_df['compound_id'].tolist()

for i in range(0, len(nat_ids), batch_size):
    batch = nat_ids[i:i+batch_size]

    # Try batch fetch
    try:
        cids = ','.join(str(int(cid)) for cid in batch if str(cid).isdigit())
        if cids:
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cids}/property/CanonicalSMILES/JSON"
            req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
            with urllib.request.urlopen(req, timeout=30) as response:
                data = json.loads(response.read().decode())

                # Create CID->SMILES map
                cid_smiles = {str(p['CID']): p['CanonicalSMILES']
                             for p in data['PropertyTable']['Properties']}

                for cid in batch:
                    nat_smiles.append(cid_smiles.get(str(int(cid)), None))
        else:
            nat_smiles.extend([None] * len(batch))
    except Exception as e:
        # Fall back to individual fetch
        for cid in batch:
            nat_smiles.append(fetch_smiles_by_cid(cid))
            time.sleep(0.05)

    if len(nat_smiles) % 500 == 0:
        success = sum(1 for s in nat_smiles if s is not None)
        print(f"  Processed {len(nat_smiles)}/{len(nat_df)} Natural products ({success} with SMILES)")

    time.sleep(0.2)

nat_df['smiles_fetched'] = nat_smiles

nat_success = sum(1 for s in nat_smiles if s is not None)
print(f"Natural SMILES fetched: {nat_success}/{len(nat_df)}")

# Combine and save
print("\n--- Combining Results ---")

# Use fetched SMILES
fda_df['final_smiles'] = fda_df['smiles_fetched']
nat_df['final_smiles'] = nat_df['smiles_fetched']

# Combine
combined_df = pd.concat([fda_df, nat_df], ignore_index=True)

# Summary
total_with_smiles = combined_df['final_smiles'].notna().sum()
print(f"\nTotal compounds with SMILES: {total_with_smiles}/{len(combined_df)}")

# Save
output_path = Path('../data/compounds_with_smiles.csv')
output_path.parent.mkdir(parents=True, exist_ok=True)
combined_df.to_csv(output_path, index=False)
print(f"\nSaved to: {output_path}")

# Also save just the compounds with SMILES for next steps
valid_df = combined_df[combined_df['final_smiles'].notna()].copy()
valid_df.to_csv('../data/valid_compounds.csv', index=False)
print(f"Valid compounds saved: {len(valid_df)}")
