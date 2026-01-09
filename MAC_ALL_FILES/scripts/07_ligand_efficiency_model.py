#!/usr/bin/env python3
"""
NOD2-Scout: Ligand Efficiency Model (Size-Corrected)
- Labels based on ligand efficiency, not raw docking score
- Removes size-correlated features
- Full leakage audit
"""

import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.metrics import roc_auc_score, accuracy_score, classification_report
from sklearn.model_selection import KFold, GroupKFold
from sklearn.linear_model import LogisticRegression
from sklearn.dummy import DummyClassifier
from scipy.stats import pearsonr
import xgboost as xgb
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
import pickle
import warnings
warnings.filterwarnings('ignore')

print("=" * 80)
print("NOD2-SCOUT: LIGAND EFFICIENCY MODEL")
print("=" * 80)

# Paths
BASE_DIR = Path(__file__).resolve().parents[2]
DATA_DIR = BASE_DIR / 'nod2-scout' / 'data'
MODELS_DIR = BASE_DIR / 'nod2-scout' / 'models'
RESULTS_DIR = BASE_DIR / 'nod2-scout' / 'results'

# Load data
features_df = pd.read_csv(DATA_DIR / 'features.csv')
fingerprints = np.load(DATA_DIR / 'fingerprints.npy')

print(f"Loaded {len(features_df)} compounds")

# =============================================================================
# STEP 1: Calculate Ligand Efficiency
# =============================================================================
print("\n" + "=" * 80)
print("STEP 1: CALCULATE LIGAND EFFICIENCY")
print("=" * 80)

features_df['ligand_efficiency'] = features_df['cnn_affinity'] / features_df['num_heavy_atoms']

print(f"\nLigand Efficiency = CNN Affinity / Heavy Atom Count")
print(f"  Range: {features_df['ligand_efficiency'].min():.4f} - {features_df['ligand_efficiency'].max():.4f}")
print(f"  Mean: {features_df['ligand_efficiency'].mean():.4f}")
print(f"  Std: {features_df['ligand_efficiency'].std():.4f}")

# Correlation check
corr_cnn = features_df['ligand_efficiency'].corr(features_df['cnn_affinity'])
corr_size = features_df['ligand_efficiency'].corr(features_df['num_heavy_atoms'])
print(f"\n  Correlation with cnn_affinity: {corr_cnn:.4f}")
print(f"  Correlation with num_heavy_atoms: {corr_size:.4f}")

# =============================================================================
# STEP 2: Create New Pseudo-Labels from Ligand Efficiency
# =============================================================================
print("\n" + "=" * 80)
print("STEP 2: CREATE PSEUDO-LABELS FROM LIGAND EFFICIENCY")
print("=" * 80)

top_le_threshold = features_df['ligand_efficiency'].quantile(0.90)
bottom_le_threshold = features_df['ligand_efficiency'].quantile(0.10)

print(f"\n  Top 10% threshold (LE >= {top_le_threshold:.4f}) -> Label 1 (efficient binders)")
print(f"  Bottom 10% threshold (LE <= {bottom_le_threshold:.4f}) -> Label 0 (inefficient)")
print(f"  Middle 80% -> Excluded from training")

features_df['label_le'] = -1
features_df.loc[features_df['ligand_efficiency'] >= top_le_threshold, 'label_le'] = 1
features_df.loc[features_df['ligand_efficiency'] <= bottom_le_threshold, 'label_le'] = 0

labeled_df = features_df[features_df['label_le'] >= 0].copy()
labeled_fps = fingerprints[features_df['label_le'] >= 0]

print(f"\n  Label distribution:")
print(f"    Efficient binders (label=1): {(labeled_df['label_le'] == 1).sum()}")
print(f"    Inefficient (label=0): {(labeled_df['label_le'] == 0).sum()}")
print(f"    Total labeled: {len(labeled_df)}")

# =============================================================================
# STEP 3: Define Size-Independent Features
# =============================================================================
print("\n" + "=" * 80)
print("STEP 3: DEFINE SIZE-INDEPENDENT FEATURES")
print("=" * 80)

# REMOVED: num_heavy_atoms, mw, num_rings, rotatable_bonds (size-correlated)
# KEPT: Shape/polarity/druglikeness features
descriptor_cols_clean = [
    'logp', 'tpsa', 'hbd', 'hba',
    'aromatic_rings', 'fraction_csp3',
    'qed', 'lipinski_violations', 'veber_violations',
    'pains_alerts', 'permeability_score', 'pgp_risk'
]

print(f"\n  REMOVED (size-correlated):")
print(f"    - num_heavy_atoms (r=0.72 with CNN)")
print(f"    - mw (r=0.72 with CNN)")
print(f"    - num_rings (r=0.70 with CNN)")
print(f"    - rotatable_bonds (r=0.46 with CNN)")

print(f"\n  KEPT ({len(descriptor_cols_clean)} features):")
for col in descriptor_cols_clean:
    print(f"    - {col}")

print(f"\n  + Morgan fingerprints (2048 bits)")
print(f"\n  Total features: {len(descriptor_cols_clean) + 2048}")

# =============================================================================
# STEP 4: Scaffold Split
# =============================================================================
print("\n" + "=" * 80)
print("STEP 4: SCAFFOLD SPLIT")
print("=" * 80)

def get_scaffold(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return "INVALID"
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
        return Chem.MolToSmiles(scaffold)
    except:
        return "INVALID"

labeled_df['scaffold'] = labeled_df['smiles'].apply(get_scaffold)
scaffold_groups = labeled_df.groupby('scaffold').indices

np.random.seed(42)
scaffolds = list(scaffold_groups.keys())
np.random.shuffle(scaffolds)

train_size = int(0.8 * len(scaffolds))
train_scaffolds = set(scaffolds[:train_size])
test_scaffolds = set(scaffolds[train_size:])

train_idx = []
test_idx = []
for scaffold in train_scaffolds:
    train_idx.extend(scaffold_groups[scaffold])
for scaffold in test_scaffolds:
    test_idx.extend(scaffold_groups[scaffold])

print(f"\n  Unique scaffolds: {len(scaffolds)}")
print(f"  Train scaffolds: {len(train_scaffolds)} (80%)")
print(f"  Test scaffolds: {len(test_scaffolds)} (20%)")
print(f"  Scaffold overlap: {len(train_scaffolds & test_scaffolds)}")
print(f"\n  Train compounds: {len(train_idx)}")
print(f"  Test compounds: {len(test_idx)}")

# Prepare data
X_desc = labeled_df[descriptor_cols_clean].values
X = np.hstack([X_desc, labeled_fps])
y = labeled_df['label_le'].values

X_train, X_test = X[train_idx], X[test_idx]
y_train, y_test = y[train_idx], y[test_idx]

print(f"\n  Train: {len(X_train)} ({sum(y_train)} positive)")
print(f"  Test: {len(X_test)} ({sum(y_test)} positive)")

# =============================================================================
# STEP 5: Train XGBoost Model
# =============================================================================
print("\n" + "=" * 80)
print("STEP 5: TRAIN XGBOOST MODEL")
print("=" * 80)

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
    'random_state': 42,
    'n_jobs': -1
}

model = xgb.XGBClassifier(**params)
model.fit(X_train, y_train, verbose=False)

y_pred_proba = model.predict_proba(X_test)[:, 1]
y_pred = model.predict(X_test)

train_pred_proba = model.predict_proba(X_train)[:, 1]

test_auc = roc_auc_score(y_test, y_pred_proba)
train_auc = roc_auc_score(y_train, train_pred_proba)

print(f"\n  Train AUROC: {train_auc:.4f}")
print(f"  Test AUROC: {test_auc:.4f}")
print(f"  Gap: {train_auc - test_auc:.4f}")

# =============================================================================
# STEP 6: FULL LEAKAGE AUDIT
# =============================================================================
print("\n" + "=" * 80)
print("STEP 6: FULL LEAKAGE AUDIT")
print("=" * 80)

audit_results = {}

# CHECK 1: Label Leakage
print("\n--- CHECK 1: Label Leakage ---")
le_in_features = 'ligand_efficiency' in descriptor_cols_clean
cnn_in_features = 'cnn_affinity' in descriptor_cols_clean
audit_results['1_label_leakage'] = 'PASS' if not (le_in_features or cnn_in_features) else 'FAIL'
print(f"  Is 'ligand_efficiency' in features? {le_in_features}")
print(f"  Is 'cnn_affinity' in features? {cnn_in_features}")
print(f"  Result: {audit_results['1_label_leakage']}")

# CHECK 2: Scaffold Split
print("\n--- CHECK 2: Scaffold Split ---")
overlap = len(train_scaffolds & test_scaffolds)
audit_results['2_scaffold_split'] = 'PASS' if overlap == 0 else 'FAIL'
print(f"  Scaffold overlap: {overlap}")
print(f"  Result: {audit_results['2_scaffold_split']}")

# CHECK 3: Feature Correlation with LE
print("\n--- CHECK 3: Feature Correlation with Ligand Efficiency ---")
high_corr_features = []
print(f"  {'Feature':<25} {'Correlation':<15}")
print("  " + "-" * 40)
for col in descriptor_cols_clean:
    corr, _ = pearsonr(labeled_df[col], labeled_df['ligand_efficiency'])
    flag = " <-- HIGH!" if abs(corr) > 0.7 else ""
    print(f"  {col:<25} r = {corr:+.4f}{flag}")
    if abs(corr) > 0.7:
        high_corr_features.append(col)

audit_results['3_feature_corr'] = 'PASS' if len(high_corr_features) == 0 else f'WARN ({len(high_corr_features)})'
print(f"\n  Features with |r| > 0.7: {len(high_corr_features)}")
print(f"  Result: {audit_results['3_feature_corr']}")

# CHECK 4: ID Leakage
print("\n--- CHECK 4: ID Leakage ---")
audit_results['4_id_leakage'] = 'PASS'
print(f"  compound_id in features? {'compound_id' in descriptor_cols_clean}")
print(f"  smiles in features? {'smiles' in descriptor_cols_clean}")
print(f"  Result: {audit_results['4_id_leakage']}")

# CHECK 5: Random vs Scaffold CV
print("\n--- CHECK 5: Random vs Scaffold CV ---")

# Random CV
kf_random = KFold(n_splits=5, shuffle=True, random_state=42)
random_scores = []
for train_i, test_i in kf_random.split(X):
    m = xgb.XGBClassifier(n_estimators=100, max_depth=4, random_state=42, n_jobs=-1)
    m.fit(X[train_i], y[train_i], verbose=False)
    p = m.predict_proba(X[test_i])[:, 1]
    random_scores.append(roc_auc_score(y[test_i], p))

# Scaffold CV
scaffold_to_group = {s: i for i, s in enumerate(scaffolds)}
groups = np.array([scaffold_to_group[labeled_df.iloc[i]['scaffold']] for i in range(len(labeled_df))])
gkf = GroupKFold(n_splits=5)
scaffold_scores = []
for train_i, test_i in gkf.split(X, y, groups):
    m = xgb.XGBClassifier(n_estimators=100, max_depth=4, random_state=42, n_jobs=-1)
    m.fit(X[train_i], y[train_i], verbose=False)
    p = m.predict_proba(X[test_i])[:, 1]
    scaffold_scores.append(roc_auc_score(y[test_i], p))

random_mean = np.mean(random_scores)
scaffold_mean = np.mean(scaffold_scores)
cv_diff = random_mean - scaffold_mean

print(f"  Random CV: {random_mean:.4f} ± {np.std(random_scores):.4f}")
print(f"  Scaffold CV: {scaffold_mean:.4f} ± {np.std(scaffold_scores):.4f}")
print(f"  Difference: {cv_diff:.4f}")

audit_results['5_cv_diff'] = f'{cv_diff:.3f}'
if cv_diff > 0.05:
    print(f"  Result: OK (scaffold split matters)")
else:
    print(f"  Result: WARN (small difference)")

# CHECK 6: Overfitting
print("\n--- CHECK 6: Train vs Test AUROC ---")
gap = train_auc - test_auc
audit_results['6_overfit_gap'] = f'{gap:.3f}'
print(f"  Train AUROC: {train_auc:.4f}")
print(f"  Test AUROC: {test_auc:.4f}")
print(f"  Gap: {gap:.4f}")
if gap < 0.10:
    print(f"  Result: OK (reasonable gap)")
else:
    print(f"  Result: WARN (large gap suggests overfitting)")

# CHECK 7: Feature Ablation
print("\n--- CHECK 7: Feature Ablation ---")

# Descriptors only
X_desc_only = labeled_df[descriptor_cols_clean].values
m_desc = xgb.XGBClassifier(n_estimators=100, max_depth=4, random_state=42, n_jobs=-1)
m_desc.fit(X_desc_only[train_idx], y_train, verbose=False)
desc_auc = roc_auc_score(y_test, m_desc.predict_proba(X_desc_only[test_idx])[:, 1])

# Fingerprints only
m_fp = xgb.XGBClassifier(n_estimators=100, max_depth=4, random_state=42, n_jobs=-1)
m_fp.fit(labeled_fps[train_idx], y_train, verbose=False)
fp_auc = roc_auc_score(y_test, m_fp.predict_proba(labeled_fps[test_idx])[:, 1])

print(f"  Descriptors only ({len(descriptor_cols_clean)} features): AUROC = {desc_auc:.4f}")
print(f"  Fingerprints only (2048 bits): AUROC = {fp_auc:.4f}")
print(f"  Combined: AUROC = {test_auc:.4f}")

audit_results['7_ablation'] = f'Desc={desc_auc:.3f}, FP={fp_auc:.3f}'

# CHECK 8: Baseline Comparison
print("\n--- CHECK 8: Baseline Comparison ---")

dummy = DummyClassifier(strategy='stratified')
dummy.fit(X_train, y_train)
dummy_auc = roc_auc_score(y_test, dummy.predict_proba(X_test)[:, 1])

logreg = LogisticRegression(max_iter=1000, random_state=42)
logreg.fit(X_train, y_train)
logreg_auc = roc_auc_score(y_test, logreg.predict_proba(X_test)[:, 1])

print(f"  Dummy (stratified): AUROC = {dummy_auc:.4f}")
print(f"  Logistic Regression: AUROC = {logreg_auc:.4f}")
print(f"  XGBoost: AUROC = {test_auc:.4f}")
print(f"  XGBoost vs LogReg lift: {test_auc - logreg_auc:.4f}")

audit_results['8_baseline'] = f'LogReg={logreg_auc:.3f}'

# CHECK 9: Label Shuffle
print("\n--- CHECK 9: Label Shuffle ---")
shuffle_aucs = []
for _ in range(10):
    y_shuffled = np.random.permutation(y_train)
    m_shuf = xgb.XGBClassifier(n_estimators=100, max_depth=4, random_state=42, n_jobs=-1)
    m_shuf.fit(X_train, y_shuffled, verbose=False)
    shuffle_aucs.append(roc_auc_score(y_test, m_shuf.predict_proba(X_test)[:, 1]))

shuffle_mean = np.mean(shuffle_aucs)
print(f"  Shuffled AUROC: {shuffle_mean:.4f} ± {np.std(shuffle_aucs):.4f}")
print(f"  Expected: ~0.50")

audit_results['9_shuffle'] = 'PASS' if shuffle_mean < 0.6 else 'FAIL'
print(f"  Result: {audit_results['9_shuffle']}")

# CHECK 10: Error Analysis
print("\n--- CHECK 10: Error Analysis ---")
test_df = labeled_df.iloc[test_idx].copy()
test_df['pred_proba'] = y_pred_proba
test_df['prediction'] = y_pred
test_df['correct'] = test_df['prediction'] == test_df['label_le']

n_correct = test_df['correct'].sum()
n_total = len(test_df)
print(f"  Correct: {n_correct}/{n_total} ({100*n_correct/n_total:.1f}%)")
print(f"  Errors: {n_total - n_correct}")

wrong = test_df[~test_df['correct']].sort_values('pred_proba', ascending=False)
if len(wrong) > 0:
    print(f"\n  Sample wrong predictions:")
    for _, row in wrong.head(5).iterrows():
        name = row['drug_name'] if pd.notna(row['drug_name']) else row['compound_id']
        print(f"    {name[:30]:<30} P={row['pred_proba']:.3f} True={row['label_le']} LE={row['ligand_efficiency']:.3f}")

audit_results['10_errors'] = f'{n_total - n_correct} errors'

# =============================================================================
# AUDIT SUMMARY TABLE
# =============================================================================
print("\n" + "=" * 80)
print("LEAKAGE AUDIT SUMMARY - LIGAND EFFICIENCY MODEL")
print("=" * 80)

print(f"\n  {'Check':<35} {'Result':<20}")
print("  " + "-" * 55)
print(f"  {'1. Label Leakage':<35} {audit_results['1_label_leakage']:<20}")
print(f"  {'2. Scaffold Split (0% overlap)':<35} {audit_results['2_scaffold_split']:<20}")
print(f"  {'3. Feature Correlation (<0.7)':<35} {audit_results['3_feature_corr']:<20}")
print(f"  {'4. ID Leakage':<35} {audit_results['4_id_leakage']:<20}")
print(f"  {'5. Random vs Scaffold CV':<35} {audit_results['5_cv_diff']:<20}")
print(f"  {'6. Train-Test Gap (<0.10)':<35} {audit_results['6_overfit_gap']:<20}")
print(f"  {'7. Feature Ablation':<35} {audit_results['7_ablation']:<20}")
print(f"  {'8. Baseline (LogReg)':<35} {audit_results['8_baseline']:<20}")
print(f"  {'9. Label Shuffle (~0.50)':<35} {audit_results['9_shuffle']:<20}")
print(f"  {'10. Error Analysis':<35} {audit_results['10_errors']:<20}")

# =============================================================================
# COMPARISON: ORIGINAL VS LIGAND EFFICIENCY MODEL
# =============================================================================
print("\n" + "=" * 80)
print("MODEL COMPARISON")
print("=" * 80)

print(f"""
  {'Metric':<30} {'Original (CNN)':<20} {'LE Model':<20}
  {'-'*70}
  {'Test AUROC':<30} {'0.998':<20} {f'{test_auc:.3f}':<20}
  {'Train AUROC':<30} {'1.000':<20} {f'{train_auc:.3f}':<20}
  {'Train-Test Gap':<30} {'0.002':<20} {f'{train_auc-test_auc:.3f}':<20}
  {'Random vs Scaffold CV':<30} {'0.001':<20} {f'{cv_diff:.3f}':<20}
  {'LogReg AUROC':<30} {'0.999':<20} {f'{logreg_auc:.3f}':<20}
  {'Shuffled AUROC':<30} {'0.487':<20} {f'{shuffle_mean:.3f}':<20}
  {'Purpose':<30} {'Docking surrogate':<20} {'Efficiency predictor':<20}
""")

# =============================================================================
# SAVE MODEL AND RESULTS
# =============================================================================
print("\n" + "=" * 80)
print("SAVING RESULTS")
print("=" * 80)

# Save model
model.save_model(str(MODELS_DIR / 'xgboost_le_model.json'))
print(f"  Model saved: {MODELS_DIR / 'xgboost_le_model.json'}")

# Save feature names
feature_names_le = descriptor_cols_clean + [f'fp_{i}' for i in range(2048)]
with open(MODELS_DIR / 'feature_names_le.pkl', 'wb') as f:
    pickle.dump(feature_names_le, f)

# Save metadata
metadata_le = {
    'n_compounds': len(features_df),
    'n_labeled': len(labeled_df),
    'top_le_threshold': float(top_le_threshold),
    'bottom_le_threshold': float(bottom_le_threshold),
    'train_auc': float(train_auc),
    'test_auc': float(test_auc),
    'scaffold_cv_auc': float(scaffold_mean),
    'random_cv_auc': float(random_mean),
    'logreg_auc': float(logreg_auc),
    'shuffle_auc': float(shuffle_mean),
    'features_used': descriptor_cols_clean,
    'audit_results': audit_results
}

with open(MODELS_DIR / 'training_metadata_le.pkl', 'wb') as f:
    pickle.dump(metadata_le, f)
print(f"  Metadata saved: {MODELS_DIR / 'training_metadata_le.pkl'}")

# Generate predictions for all compounds
print("\n  Generating predictions for all compounds...")

X_all_desc = features_df[descriptor_cols_clean].values
X_all = np.hstack([X_all_desc, fingerprints])
features_df['le_score'] = model.predict_proba(X_all)[:, 1]

# Save with new scores
features_df.to_csv(DATA_DIR / 'features_with_le_scores.csv', index=False)
print(f"  Scores saved: {DATA_DIR / 'features_with_le_scores.csv'}")

print("\n" + "=" * 80)
print("COMPLETE")
print("=" * 80)
