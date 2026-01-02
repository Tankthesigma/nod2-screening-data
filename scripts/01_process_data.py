#!/usr/bin/env python3
"""
NOD2-Scout Step 1: Process screening data and fetch missing SMILES
"""

import pandas as pd
import numpy as np
import urllib.request
import json
import time
from pathlib import Path
import re

print("=" * 70)
print("NOD2-SCOUT: DATA PROCESSING")
print("=" * 70)

# Paths
BASE_DIR = Path(__file__).resolve().parents[2]
DATA_DIR = BASE_DIR / 'nod2-screening-data'
OUTPUT_DIR = BASE_DIR / 'nod2-scout' / 'data'

# Load data
df = pd.read_csv(DATA_DIR / 'screening_results.csv')
print(f"Loaded {len(df)} compounds")

def is_valid_smiles(s):
    """Check if string looks like a SMILES"""
    if pd.isna(s):
        return False
    s = str(s).strip()
    # SMILES typically has these characteristics
    if len(s) < 3:
        return False
    if s.startswith('_') or s.startswith('M '):
        return False
    if s.upper() == s and not any(c in s for c in '()=#[]'):
        # All uppercase with no SMILES characters = probably drug name
        return False
    # Contains typical SMILES characters
    return bool(re.search(r'[CcNnOoSsPpFfClBrI]', s))

def extract_drug_name(smiles_field):
    """Extract drug name from SMILES field"""
    if pd.isna(smiles_field):
        return None
    s = str(smiles_field).strip()
    if s.startswith('_c') or s == 'Molecule-1':
        return None
    if s.startswith('M '):
        return s[2:].strip()
    # If all uppercase and looks like a name
    if s.upper() == s and len(s) > 2:
        return s
    return None

def fetch_smiles_by_name(name, timeout=10):
    """Fetch SMILES from PubChem by drug name"""
    try:
        name = str(name).strip()
        # Clean name - remove common suffixes
        name = re.sub(r'\s+(HYDROCHLORIDE|HCL|SODIUM|POTASSIUM|ACETATE|SULFATE|MALEATE|TARTRATE|PHOSPHATE|MESYLATE|BESYLATE)$', '', name, flags=re.I)
        name = re.sub(r'_c\d+$', '', name)  # Remove _c001 etc

        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{urllib.request.quote(name)}/property/CanonicalSMILES/JSON"
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req, timeout=timeout) as response:
            data = json.loads(response.read().decode())
            return data['PropertyTable']['Properties'][0]['CanonicalSMILES']
    except:
        return None

def fetch_smiles_batch_cid(cids, timeout=30):
    """Batch fetch SMILES from PubChem by CIDs"""
    try:
        cid_str = ','.join(str(c) for c in cids)
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid_str}/property/CanonicalSMILES/JSON"
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req, timeout=timeout) as response:
            data = json.loads(response.read().decode())
            return {str(p['CID']): p['CanonicalSMILES'] for p in data['PropertyTable']['Properties']}
    except:
        return {}

# Separate by type
fda_df = df[df['type'] == 'FDA'].copy()
nat_df = df[df['type'] == 'Natural'].copy()

# Process FDA drugs
print("\n--- Processing FDA Drugs ---")
fda_df['has_smiles'] = fda_df['smiles'].apply(is_valid_smiles)
fda_df['drug_name'] = fda_df['smiles'].apply(extract_drug_name)
fda_df['final_smiles'] = fda_df.apply(
    lambda row: row['smiles'] if row['has_smiles'] else None, axis=1
)

# Fetch missing FDA SMILES
fda_missing = fda_df[(fda_df['final_smiles'].isna()) & (fda_df['drug_name'].notna())]
print(f"FDA with valid SMILES: {fda_df['final_smiles'].notna().sum()}")
print(f"FDA needing fetch: {len(fda_missing)}")

fetched_count = 0
for idx, row in fda_missing.iterrows():
    smiles = fetch_smiles_by_name(row['drug_name'])
    if smiles:
        fda_df.at[idx, 'final_smiles'] = smiles
        fetched_count += 1

    if fetched_count % 50 == 0 and fetched_count > 0:
        print(f"  Fetched {fetched_count} FDA SMILES...")
    time.sleep(0.15)

print(f"FDA SMILES fetched: {fetched_count}")
print(f"FDA total with SMILES: {fda_df['final_smiles'].notna().sum()}")

# Process Natural Products
print("\n--- Processing Natural Products ---")
nat_df['has_smiles'] = nat_df['smiles'].apply(is_valid_smiles)
nat_df['final_smiles'] = nat_df.apply(
    lambda row: row['smiles'] if row['has_smiles'] else None, axis=1
)

print(f"Natural with valid SMILES: {nat_df['final_smiles'].notna().sum()}")

# For missing, try compound_id as PubChem CID
nat_missing = nat_df[nat_df['final_smiles'].isna()].copy()
nat_missing_ids = nat_missing['compound_id'].tolist()

print(f"Natural needing fetch: {len(nat_missing_ids)}")

# Batch fetch
batch_size = 100
fetched_nat = 0
for i in range(0, min(len(nat_missing_ids), 2000), batch_size):  # Limit to 2000 fetches
    batch_ids = nat_missing_ids[i:i+batch_size]

    # Filter to numeric IDs only
    valid_ids = [cid for cid in batch_ids if str(cid).replace('_', '').isdigit()]

    if valid_ids:
        cid_smiles = fetch_smiles_batch_cid(valid_ids)

        for cid in valid_ids:
            smiles = cid_smiles.get(str(cid))
            if smiles:
                nat_df.loc[nat_df['compound_id'] == cid, 'final_smiles'] = smiles
                fetched_nat += 1

    if (i + batch_size) % 500 == 0:
        print(f"  Processed {i + batch_size} natural products ({fetched_nat} fetched)...")

    time.sleep(0.25)

print(f"Natural SMILES fetched: {fetched_nat}")
print(f"Natural total with SMILES: {nat_df['final_smiles'].notna().sum()}")

# Combine
print("\n--- Combining Results ---")
combined_df = pd.concat([fda_df, nat_df], ignore_index=True)

# Summary
total = len(combined_df)
with_smiles = combined_df['final_smiles'].notna().sum()
print(f"\nTotal compounds: {total}")
print(f"With SMILES: {with_smiles} ({100*with_smiles/total:.1f}%)")

# Save all data
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
combined_df.to_csv(OUTPUT_DIR / 'all_compounds.csv', index=False)
print(f"\nSaved: {OUTPUT_DIR / 'all_compounds.csv'}")

# Save only valid compounds for ML
valid_df = combined_df[combined_df['final_smiles'].notna()].copy()
valid_df = valid_df[['compound_id', 'final_smiles', 'cnn_affinity', 'type', 'drug_name']].copy()
valid_df.columns = ['compound_id', 'smiles', 'cnn_affinity', 'type', 'drug_name']
valid_df.to_csv(OUTPUT_DIR / 'valid_compounds.csv', index=False)
print(f"Saved: {OUTPUT_DIR / 'valid_compounds.csv'} ({len(valid_df)} compounds)")

# Show top hits
print("\n--- TOP 20 HITS ---")
top20 = valid_df.nlargest(20, 'cnn_affinity')
for i, row in top20.iterrows():
    name = row['drug_name'] if pd.notna(row['drug_name']) else row['compound_id']
    print(f"  {row['cnn_affinity']:.3f}: {name} ({row['type']})")
