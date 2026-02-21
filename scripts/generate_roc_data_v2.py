#!/usr/bin/env python3
"""
Generate ROC curve data using PRE-COMPUTED labels from features_corrected.csv
This ensures we match the original model exactly.
"""

import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.model_selection import GroupKFold
import xgboost as xgb
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
import json
import pickle
import warnings
warnings.filterwarnings('ignore')

print("=" * 70)
print("GENERATING ROC DATA (Using pre-computed labels)")
print("=" * 70)

# Paths
BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / 'MAC_ALL_FILES' / 'data'

# Load pre-computed corrected data
df = pd.read_csv(DATA_DIR / 'features_corrected.csv')
fps = np.load(DATA_DIR / 'fingerprints.npy')

print(f"Loaded {len(df)} compounds from features_corrected.csv")

# Use PRE-COMPUTED labels (already size-matched)
labeled_df = df[df['label'] >= 0].copy()
labeled_fps = fps[df['label'] >= 0]

print(f"Labeled compounds: {len(labeled_df)}")
print(f"  Positive: {(labeled_df['label']==1).sum()}")
print(f"  Negative: {(labeled_df['label']==0).sum()}")

# Features (same as original - no size-related features)
descriptor_cols = [
    'logp', 'tpsa', 'hbd', 'hba',
    'aromatic_rings', 'fraction_csp3',
    'qed', 'pains_alerts',
    'permeability_score', 'pgp_risk'
]

# Scaffold extraction
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

# Prepare features
X_desc = labeled_df[descriptor_cols].values
X = np.hstack([X_desc, labeled_fps])
y = labeled_df['label'].values

print(f"Feature matrix: {X.shape}")

# 5-Fold Scaffold CV (SAME params as original)
print("\n" + "=" * 70)
print("5-FOLD SCAFFOLD CV")
print("=" * 70)

gkf = GroupKFold(n_splits=5)

fold_aucs = []
fold_roc_data = []

for fold_idx, (train_idx, test_idx) in enumerate(gkf.split(X, y, groups)):
    # SAME params as original 08_corrected_model.py
    model = xgb.XGBClassifier(n_estimators=100, max_depth=4, random_state=42, n_jobs=-1)
    model.fit(X[train_idx], y[train_idx], verbose=False)

    y_pred = model.predict_proba(X[test_idx])[:, 1]
    auc = roc_auc_score(y[test_idx], y_pred)
    fold_aucs.append(auc)

    fpr, tpr, _ = roc_curve(y[test_idx], y_pred)
    fold_roc_data.append({
        'fold': fold_idx + 1,
        'auc': float(auc),
        'fpr': fpr.tolist(),
        'tpr': tpr.tolist()
    })

    print(f"  Fold {fold_idx + 1}: AUC = {auc:.4f}")

print(f"\n  Mean: {np.mean(fold_aucs):.4f} ± {np.std(fold_aucs):.4f}")
print(f"  Range: {min(fold_aucs):.2f} - {max(fold_aucs):.2f}")

# Shuffled control
print("\n" + "=" * 70)
print("SHUFFLED CONTROL")
print("=" * 70)

train_size = int(0.8 * len(scaffolds))
train_scaffolds = set(scaffolds[:train_size])
train_idx = [i for s in train_scaffolds for i in scaffold_groups[s]]
test_idx = [i for s in set(scaffolds) - train_scaffolds for i in scaffold_groups[s]]

X_train, X_test = X[train_idx], X[test_idx]
y_train, y_test = y[train_idx], y[test_idx]

shuffle_aucs = []
for i in range(10):
    y_shuf = np.random.permutation(y_train)
    model = xgb.XGBClassifier(n_estimators=50, max_depth=3, random_state=i, n_jobs=-1)
    model.fit(X_train, y_shuf, verbose=False)
    shuffle_aucs.append(roc_auc_score(y_test, model.predict_proba(X_test)[:, 1]))

print(f"  Shuffled AUC: {np.mean(shuffle_aucs):.4f} ± {np.std(shuffle_aucs):.4f}")

# Save
output_dir = BASE_DIR / 'isef_figures'
results = {
    'fold_aucs': fold_aucs,
    'mean_auc': float(np.mean(fold_aucs)),
    'std_auc': float(np.std(fold_aucs)),
    'fold_roc_data': fold_roc_data,
    'shuffle_auc': float(np.mean(shuffle_aucs))
}

with open(output_dir / 'roc_curve_data.json', 'w') as f:
    json.dump(results, f, indent=2)
with open(output_dir / 'roc_curve_data.pkl', 'wb') as f:
    pickle.dump(results, f)

print(f"\nSaved to: {output_dir / 'roc_curve_data.json'}")

print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print(f"""
FOLD AUCs:
  Fold 1: {fold_aucs[0]:.4f}
  Fold 2: {fold_aucs[1]:.4f}
  Fold 3: {fold_aucs[2]:.4f}
  Fold 4: {fold_aucs[3]:.4f}
  Fold 5: {fold_aucs[4]:.4f}

  Mean: {np.mean(fold_aucs):.4f} ± {np.std(fold_aucs):.4f}
  Range: {min(fold_aucs):.2f} - {max(fold_aucs):.2f}

SHUFFLED: {np.mean(shuffle_aucs):.2f}
""")
