#!/usr/bin/env python3
"""
Generate ROC curve data and individual fold AUCs for ISEF poster.
Uses the CORRECTED size-matched model (not leaked).

Output:
- Individual fold AUCs (5 values)
- ROC curve points (FPR, TPR) for each fold
- Shuffled control ROC data
"""

import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.model_selection import GroupKFold, KFold
import xgboost as xgb
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
import pickle
import json
import warnings
warnings.filterwarnings('ignore')

print("=" * 70)
print("GENERATING ROC CURVE DATA FOR ISEF POSTER")
print("=" * 70)

# Paths - try multiple locations
BASE_DIR = Path(__file__).resolve().parent.parent
DATA_PATHS = [
    BASE_DIR / 'MAC_ALL_FILES' / 'data',
    BASE_DIR / 'nod2-scout' / 'data',
    BASE_DIR / 'data'
]

# Find data directory
DATA_DIR = None
for p in DATA_PATHS:
    if (p / 'features.csv').exists():
        DATA_DIR = p
        break

if DATA_DIR is None:
    print("ERROR: Could not find features.csv")
    print("Searched:", DATA_PATHS)
    exit(1)

print(f"Using data from: {DATA_DIR}")

# Load data
df = pd.read_csv(DATA_DIR / 'features.csv')
fps = np.load(DATA_DIR / 'fingerprints.npy')
print(f"Loaded {len(df)} compounds")

# =============================================================================
# STEP 1: Create Size-Matched Labels (CORRECTED - not leaked)
# =============================================================================
print("\n" + "=" * 70)
print("STEP 1: SIZE-MATCHED PSEUDO-LABELS (CORRECTED METHOD)")
print("=" * 70)

# Create MW bins (quintiles)
df['mw_bin'] = pd.qcut(df['mw'], q=5, labels=['XS', 'S', 'M', 'L', 'XL'])

# Within each bin, label top 20% and bottom 20% by cnn_affinity
df['label'] = -1

for bin_name in ['XS', 'S', 'M', 'L', 'XL']:
    mask = df['mw_bin'] == bin_name
    bin_data = df[mask]

    top_thresh = bin_data['cnn_affinity'].quantile(0.80)
    bottom_thresh = bin_data['cnn_affinity'].quantile(0.20)

    df.loc[mask & (df['cnn_affinity'] >= top_thresh), 'label'] = 1
    df.loc[mask & (df['cnn_affinity'] <= bottom_thresh), 'label'] = 0

labeled_df = df[df['label'] >= 0].copy()
labeled_fps = fps[df['label'] >= 0]

print(f"Labeled compounds: {len(labeled_df)}")
print(f"  Positive (top 20% per bin): {(labeled_df['label']==1).sum()}")
print(f"  Negative (bottom 20% per bin): {(labeled_df['label']==0).sum()}")

# =============================================================================
# STEP 2: Features (size-independent to prevent leakage)
# =============================================================================
descriptor_cols = [
    'logp', 'tpsa', 'hbd', 'hba',
    'aromatic_rings', 'fraction_csp3',
    'qed', 'pains_alerts',
    'permeability_score', 'pgp_risk'
]

print(f"\nFeatures: {len(descriptor_cols)} descriptors + 2048 fingerprint bits")

# =============================================================================
# STEP 3: Scaffold extraction
# =============================================================================
print("\n" + "=" * 70)
print("STEP 2: SCAFFOLD EXTRACTION")
print("=" * 70)

def get_scaffold(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return "INVALID"
        return Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(mol))
    except:
        return "INVALID"

labeled_df['scaffold'] = labeled_df['smiles'].apply(get_scaffold)
scaffold_groups = labeled_df.groupby('scaffold').indices

np.random.seed(42)
scaffolds = list(scaffold_groups.keys())
np.random.shuffle(scaffolds)

scaffold_to_group = {s: i for i, s in enumerate(scaffolds)}
groups = np.array([scaffold_to_group[labeled_df.iloc[i]['scaffold']] for i in range(len(labeled_df))])

print(f"Unique scaffolds: {len(scaffolds)}")

# Prepare feature matrix
X_desc = labeled_df[descriptor_cols].values
X = np.hstack([X_desc, labeled_fps])
y = labeled_df['label'].values

print(f"Feature matrix shape: {X.shape}")

# =============================================================================
# STEP 4: 5-Fold Scaffold CV with ROC curves
# =============================================================================
print("\n" + "=" * 70)
print("STEP 3: 5-FOLD SCAFFOLD CV WITH ROC DATA")
print("=" * 70)

# Use SAME params as original 08_corrected_model.py CV loop
params = {
    'n_estimators': 100,
    'max_depth': 4,
    'random_state': 42,
    'n_jobs': -1
}

gkf = GroupKFold(n_splits=5)

fold_aucs = []
fold_roc_data = []

print("\nRunning 5-fold scaffold-split CV...")
for fold_idx, (train_idx, test_idx) in enumerate(gkf.split(X, y, groups)):
    # Train model
    model = xgb.XGBClassifier(**params)
    model.fit(X[train_idx], y[train_idx], verbose=False)

    # Predict
    y_pred = model.predict_proba(X[test_idx])[:, 1]

    # Calculate AUC
    auc = roc_auc_score(y[test_idx], y_pred)
    fold_aucs.append(auc)

    # Calculate ROC curve points
    fpr, tpr, thresholds = roc_curve(y[test_idx], y_pred)

    fold_roc_data.append({
        'fold': fold_idx + 1,
        'auc': float(auc),
        'fpr': fpr.tolist(),
        'tpr': tpr.tolist(),
        'n_test': len(test_idx)
    })

    print(f"  Fold {fold_idx + 1}: AUC = {auc:.4f} (n={len(test_idx)})")

print(f"\n  Mean AUC: {np.mean(fold_aucs):.4f} ± {np.std(fold_aucs):.4f}")
print(f"  Range: {min(fold_aucs):.2f} - {max(fold_aucs):.2f}")

# =============================================================================
# STEP 5: Shuffled control
# =============================================================================
print("\n" + "=" * 70)
print("STEP 4: SHUFFLED LABEL CONTROL")
print("=" * 70)

# Use simple train/test split for shuffle control
train_size = int(0.8 * len(scaffolds))
train_scaffolds = set(scaffolds[:train_size])

train_idx = []
test_idx = []
for s in train_scaffolds:
    train_idx.extend(scaffold_groups[s])
for s in set(scaffolds) - train_scaffolds:
    test_idx.extend(scaffold_groups[s])

X_train, X_test = X[train_idx], X[test_idx]
y_train, y_test = y[train_idx], y[test_idx]

shuffle_aucs = []
shuffle_roc_data = []

print("Running 10 shuffle iterations...")
for i in range(10):
    y_shuf = np.random.permutation(y_train)
    model = xgb.XGBClassifier(n_estimators=50, max_depth=3, random_state=i, n_jobs=-1)
    model.fit(X_train, y_shuf, verbose=False)
    y_pred = model.predict_proba(X_test)[:, 1]

    auc = roc_auc_score(y_test, y_pred)
    shuffle_aucs.append(auc)

    if i == 0:  # Save first shuffle ROC for plotting
        fpr, tpr, _ = roc_curve(y_test, y_pred)
        shuffle_roc_data = {'fpr': fpr.tolist(), 'tpr': tpr.tolist(), 'auc': float(auc)}

print(f"  Shuffled AUC: {np.mean(shuffle_aucs):.4f} ± {np.std(shuffle_aucs):.4f}")

# =============================================================================
# STEP 6: Save results
# =============================================================================
print("\n" + "=" * 70)
print("STEP 5: SAVING RESULTS")
print("=" * 70)

output_dir = BASE_DIR / 'isef_figures'
output_dir.mkdir(exist_ok=True)

results = {
    'model_type': 'XGBoost (corrected, size-matched labels)',
    'n_labeled': len(labeled_df),
    'n_features': X.shape[1],
    'fold_aucs': fold_aucs,
    'mean_auc': float(np.mean(fold_aucs)),
    'std_auc': float(np.std(fold_aucs)),
    'min_auc': float(min(fold_aucs)),
    'max_auc': float(max(fold_aucs)),
    'fold_roc_data': fold_roc_data,
    'shuffle_mean_auc': float(np.mean(shuffle_aucs)),
    'shuffle_std_auc': float(np.std(shuffle_aucs)),
    'shuffle_roc_data': shuffle_roc_data
}

# Save as JSON
output_file = output_dir / 'roc_curve_data.json'
with open(output_file, 'w') as f:
    json.dump(results, f, indent=2)
print(f"Saved: {output_file}")

# Also save as pickle for easier Python loading
output_pkl = output_dir / 'roc_curve_data.pkl'
with open(output_pkl, 'wb') as f:
    pickle.dump(results, f)
print(f"Saved: {output_pkl}")

# =============================================================================
# SUMMARY
# =============================================================================
print("\n" + "=" * 70)
print("SUMMARY - DATA FOR ISEF POSTER")
print("=" * 70)

print(f"""
INDIVIDUAL FOLD AUCs (5-fold scaffold-split CV):
  Fold 1: {fold_aucs[0]:.4f}
  Fold 2: {fold_aucs[1]:.4f}
  Fold 3: {fold_aucs[2]:.4f}
  Fold 4: {fold_aucs[3]:.4f}
  Fold 5: {fold_aucs[4]:.4f}

  Mean: {np.mean(fold_aucs):.4f} ± {np.std(fold_aucs):.4f}
  Range: {min(fold_aucs):.2f} - {max(fold_aucs):.2f}

SHUFFLED CONTROL:
  Mean AUC: {np.mean(shuffle_aucs):.4f} ± {np.std(shuffle_aucs):.4f}

ROC CURVE DATA:
  Saved to: {output_file}
  Contains FPR/TPR points for all 5 folds + shuffled control
""")

print("Done!")
