#!/usr/bin/env python3
"""
NOD2-Scout Step 3: Train XGBoost Classifier
- Pseudo-labels from top/bottom 10% by CNN affinity
- Scaffold split for train/test
- Negative control (label shuffle)
"""

import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.metrics import roc_auc_score, precision_recall_curve, auc
from sklearn.metrics import confusion_matrix, classification_report
import xgboost as xgb
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
import pickle
import warnings
warnings.filterwarnings('ignore')

print("=" * 70)
print("NOD2-SCOUT: MODEL TRAINING")
print("=" * 70)

# Paths
BASE_DIR = Path(__file__).resolve().parents[2]
DATA_DIR = BASE_DIR / 'nod2-scout' / 'data'
MODELS_DIR = BASE_DIR / 'nod2-scout' / 'models'
MODELS_DIR.mkdir(parents=True, exist_ok=True)

# Load data
features_df = pd.read_csv(DATA_DIR / 'features.csv')
fingerprints = np.load(DATA_DIR / 'fingerprints.npy')

print(f"Loaded {len(features_df)} compounds")
print(f"Fingerprint shape: {fingerprints.shape}")

# =============================================================================
# STEP 1: Create Pseudo-Labels
# =============================================================================
print("\n--- CREATING PSEUDO-LABELS ---")

# Top 10% = high binding propensity (label 1)
# Bottom 10% = low binding propensity (label 0)
top_threshold = features_df['cnn_affinity'].quantile(0.90)
bottom_threshold = features_df['cnn_affinity'].quantile(0.10)

print(f"Top 10% threshold: {top_threshold:.3f}")
print(f"Bottom 10% threshold: {bottom_threshold:.3f}")

# Create labels
features_df['label'] = -1  # -1 = middle (excluded)
features_df.loc[features_df['cnn_affinity'] >= top_threshold, 'label'] = 1
features_df.loc[features_df['cnn_affinity'] <= bottom_threshold, 'label'] = 0

# Filter to labeled compounds only
labeled_mask = features_df['label'] >= 0
labeled_df = features_df[labeled_mask].copy()
labeled_fps = fingerprints[labeled_mask]

print(f"High binding (label=1): {(labeled_df['label'] == 1).sum()}")
print(f"Low binding (label=0): {(labeled_df['label'] == 0).sum()}")
print(f"Total labeled: {len(labeled_df)}")

# =============================================================================
# STEP 2: Scaffold Split
# =============================================================================
print("\n--- SCAFFOLD SPLIT ---")

def get_scaffold(smiles):
    """Get Murcko scaffold for a molecule"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
        return Chem.MolToSmiles(scaffold)
    except:
        return None

# Generate scaffolds
labeled_df['scaffold'] = labeled_df['smiles'].apply(get_scaffold)

# Group by scaffold
scaffold_groups = labeled_df.groupby('scaffold').indices

# Split scaffolds into train/test (80/20)
np.random.seed(42)
scaffolds = list(scaffold_groups.keys())
np.random.shuffle(scaffolds)

train_size = int(0.8 * len(scaffolds))
train_scaffolds = set(scaffolds[:train_size])
test_scaffolds = set(scaffolds[train_size:])

# Assign compounds to train/test
train_idx = []
test_idx = []

for scaffold in train_scaffolds:
    train_idx.extend(scaffold_groups[scaffold])
for scaffold in test_scaffolds:
    test_idx.extend(scaffold_groups[scaffold])

print(f"Unique scaffolds: {len(scaffolds)}")
print(f"Train scaffolds: {len(train_scaffolds)}")
print(f"Test scaffolds: {len(test_scaffolds)}")

# Get features and labels
descriptor_cols = ['mw', 'logp', 'tpsa', 'hbd', 'hba', 'rotatable_bonds',
                   'aromatic_rings', 'fraction_csp3', 'num_heavy_atoms', 'num_rings',
                   'qed', 'lipinski_violations', 'veber_violations', 'pains_alerts',
                   'permeability_score', 'pgp_risk']

X_desc = labeled_df[descriptor_cols].values
X_fp = labeled_fps
X = np.hstack([X_desc, X_fp])  # Combine descriptors and fingerprints
y = labeled_df['label'].values

X_train = X[train_idx]
X_test = X[test_idx]
y_train = y[train_idx]
y_test = y[test_idx]

print(f"\nTrain set: {len(X_train)} compounds ({sum(y_train)} positive)")
print(f"Test set: {len(X_test)} compounds ({sum(y_test)} positive)")

# =============================================================================
# STEP 3: Train XGBoost Model
# =============================================================================
print("\n--- TRAINING XGBOOST MODEL ---")

# Model parameters
params = {
    'n_estimators': 500,
    'max_depth': 6,
    'learning_rate': 0.05,
    'subsample': 0.8,
    'colsample_bytree': 0.8,
    'min_child_weight': 3,
    'gamma': 0.1,
    'reg_alpha': 0.1,
    'reg_lambda': 1.0,
    'objective': 'binary:logistic',
    'eval_metric': 'auc',
    'random_state': 42,
    'n_jobs': -1
}

model = xgb.XGBClassifier(**params)
model.fit(
    X_train, y_train,
    eval_set=[(X_test, y_test)],
    verbose=False
)

# Predictions
y_pred_proba = model.predict_proba(X_test)[:, 1]
y_pred = model.predict(X_test)

# Metrics
auroc = roc_auc_score(y_test, y_pred_proba)
precision, recall, _ = precision_recall_curve(y_test, y_pred_proba)
auprc = auc(recall, precision)

print(f"\n--- MODEL PERFORMANCE ---")
print(f"AUROC: {auroc:.4f}")
print(f"AUPRC: {auprc:.4f}")
print(f"\nClassification Report:")
print(classification_report(y_test, y_pred, target_names=['Low Binding', 'High Binding']))

# Confusion matrix
cm = confusion_matrix(y_test, y_pred)
print(f"\nConfusion Matrix:")
print(cm)

# =============================================================================
# STEP 4: Negative Control (Label Shuffle)
# =============================================================================
print("\n--- NEGATIVE CONTROL (LABEL SHUFFLE) ---")

n_shuffles = 10
shuffle_aurocs = []

for i in range(n_shuffles):
    # Shuffle labels
    y_train_shuffled = np.random.permutation(y_train)

    # Train model
    model_shuffle = xgb.XGBClassifier(**params)
    model_shuffle.fit(X_train, y_train_shuffled, verbose=False)

    # Evaluate
    y_pred_shuffle = model_shuffle.predict_proba(X_test)[:, 1]
    auroc_shuffle = roc_auc_score(y_test, y_pred_shuffle)
    shuffle_aurocs.append(auroc_shuffle)

shuffle_mean = np.mean(shuffle_aurocs)
shuffle_std = np.std(shuffle_aurocs)

print(f"Shuffled AUROC: {shuffle_mean:.4f} ± {shuffle_std:.4f}")
print(f"Expected (random): ~0.50")
print(f"\nReal model lift over random: {auroc - shuffle_mean:.4f}")

# =============================================================================
# STEP 5: Generate ML Scores for All Compounds
# =============================================================================
print("\n--- GENERATING ML SCORES FOR ALL COMPOUNDS ---")

# Get features for all compounds
X_all_desc = features_df[descriptor_cols].values
X_all = np.hstack([X_all_desc, fingerprints])

# Predict
ml_scores = model.predict_proba(X_all)[:, 1]
features_df['ml_score'] = ml_scores

# Normalize CNN affinity to 0-1 scale
cnn_min = features_df['cnn_affinity'].min()
cnn_max = features_df['cnn_affinity'].max()
features_df['cnn_norm'] = (features_df['cnn_affinity'] - cnn_min) / (cnn_max - cnn_min)

print(f"ML scores generated for {len(features_df)} compounds")
print(f"ML score range: {ml_scores.min():.4f} - {ml_scores.max():.4f}")

# =============================================================================
# STEP 6: Save Results
# =============================================================================
print("\n--- SAVING RESULTS ---")

# Save model
model.save_model(str(MODELS_DIR / 'xgboost_model.json'))
print(f"Model saved: {MODELS_DIR / 'xgboost_model.json'}")

# Save feature names
feature_names = descriptor_cols + [f'fp_{i}' for i in range(2048)]
with open(MODELS_DIR / 'feature_names.pkl', 'wb') as f:
    pickle.dump(feature_names, f)
print(f"Feature names saved: {MODELS_DIR / 'feature_names.pkl'}")

# Save predictions
features_df.to_csv(DATA_DIR / 'features_with_scores.csv', index=False)
print(f"Features with scores saved: {DATA_DIR / 'features_with_scores.csv'}")

# Save training metadata
metadata = {
    'n_compounds': len(features_df),
    'n_labeled': len(labeled_df),
    'n_train': len(X_train),
    'n_test': len(X_test),
    'top_threshold': float(top_threshold),
    'bottom_threshold': float(bottom_threshold),
    'auroc': float(auroc),
    'auprc': float(auprc),
    'shuffle_auroc_mean': float(shuffle_mean),
    'shuffle_auroc_std': float(shuffle_std),
    'n_scaffolds': len(scaffolds),
    'train_scaffolds': len(train_scaffolds),
    'test_scaffolds': len(test_scaffolds)
}

with open(MODELS_DIR / 'training_metadata.pkl', 'wb') as f:
    pickle.dump(metadata, f)
print(f"Metadata saved: {MODELS_DIR / 'training_metadata.pkl'}")

# =============================================================================
# Summary
# =============================================================================
print("\n" + "=" * 70)
print("TRAINING SUMMARY")
print("=" * 70)
print(f"""
Dataset:
  Total compounds: {len(features_df)}
  Labeled (top/bottom 10%): {len(labeled_df)}
  Train set: {len(X_train)} ({sum(y_train)} positive)
  Test set: {len(X_test)} ({sum(y_test)} positive)

Scaffold Split:
  Unique scaffolds: {len(scaffolds)}
  Train scaffolds: {len(train_scaffolds)} (80%)
  Test scaffolds: {len(test_scaffolds)} (20%)

Model Performance:
  AUROC: {auroc:.4f} (target ≥ 0.70)
  AUPRC: {auprc:.4f}

Negative Control:
  Shuffled AUROC: {shuffle_mean:.4f} ± {shuffle_std:.4f}
  Model lift: {auroc - shuffle_mean:.4f}

Status: {'PASSED' if auroc >= 0.70 else 'NEEDS IMPROVEMENT'} (AUROC {'≥' if auroc >= 0.70 else '<'} 0.70)
""")
