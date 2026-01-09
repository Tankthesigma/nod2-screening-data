#!/usr/bin/env python3
"""
NOD2-Scout Step 4: Composite Drugability Score & Ablation Study
- Build composite score with weighted components
- Run ablation study (Docking-only vs ML-only vs NOD2-Scout)
- SHAP feature importance analysis
- Validate with known IBD drugs
"""

import pandas as pd
import numpy as np
from pathlib import Path
import xgboost as xgb
import shap
import pickle
import warnings
warnings.filterwarnings('ignore')

print("=" * 70)
print("NOD2-SCOUT: COMPOSITE SCORE & ABLATION STUDY")
print("=" * 70)

# Paths
BASE_DIR = Path(__file__).resolve().parents[2]
DATA_DIR = BASE_DIR / 'nod2-scout' / 'data'
MODELS_DIR = BASE_DIR / 'nod2-scout' / 'models'
RESULTS_DIR = BASE_DIR / 'nod2-scout' / 'results'
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# Load data
df = pd.read_csv(DATA_DIR / 'features_with_scores.csv')
fingerprints = np.load(DATA_DIR / 'fingerprints.npy')

print(f"Loaded {len(df)} compounds")

# =============================================================================
# STEP 1: Build Composite Drugability Score
# =============================================================================
print("\n--- BUILDING COMPOSITE SCORE ---")

# Component weights (tunable)
weights = {
    'ml_score': 0.35,          # ML prediction
    'cnn_affinity': 0.30,      # Docking score
    'qed': 0.15,               # Drug-likeness
    'permeability': 0.10,      # Membrane permeability
    'safety': 0.10             # Safety (PAINS-free, Lipinski-compliant)
}

# Normalize components to 0-1 scale
df['ml_norm'] = df['ml_score']  # Already 0-1

cnn_min = df['cnn_affinity'].min()
cnn_max = df['cnn_affinity'].max()
df['cnn_norm'] = (df['cnn_affinity'] - cnn_min) / (cnn_max - cnn_min)

df['qed_norm'] = df['qed']  # Already 0-1

df['perm_norm'] = df['permeability_score']  # Already 0-1

# Safety score: PAINS-free (1 if no alerts) and Lipinski-compliant (1 if ≤1 violation)
df['pains_score'] = (df['pains_alerts'] == 0).astype(float)
df['lipinski_score'] = (df['lipinski_violations'] <= 1).astype(float)
df['safety_norm'] = 0.5 * df['pains_score'] + 0.5 * df['lipinski_score']

# Calculate composite score
df['nod2_scout_score'] = (
    weights['ml_score'] * df['ml_norm'] +
    weights['cnn_affinity'] * df['cnn_norm'] +
    weights['qed'] * df['qed_norm'] +
    weights['permeability'] * df['perm_norm'] +
    weights['safety'] * df['safety_norm']
)

print(f"Composite score range: {df['nod2_scout_score'].min():.4f} - {df['nod2_scout_score'].max():.4f}")
print(f"Mean: {df['nod2_scout_score'].mean():.4f}")

# =============================================================================
# STEP 2: Ablation Study
# =============================================================================
print("\n--- ABLATION STUDY ---")

# Method 1: Docking Only (rank by CNN affinity)
df['rank_docking'] = df['cnn_affinity'].rank(ascending=False)

# Method 2: ML Only (rank by ML score)
df['rank_ml'] = df['ml_score'].rank(ascending=False)

# Method 3: NOD2-Scout Composite (rank by composite score)
df['rank_composite'] = df['nod2_scout_score'].rank(ascending=False)

# Analyze overlap in top candidates
top_n_values = [10, 25, 50, 100]

print("\nTop-N Overlap Analysis:")
print("-" * 60)
print(f"{'Top-N':<10} {'Dock∩ML':<12} {'Dock∩Scout':<12} {'ML∩Scout':<12}")
print("-" * 60)

for top_n in top_n_values:
    top_docking = set(df.nsmallest(top_n, 'rank_docking')['compound_id'])
    top_ml = set(df.nsmallest(top_n, 'rank_ml')['compound_id'])
    top_composite = set(df.nsmallest(top_n, 'rank_composite')['compound_id'])

    dock_ml = len(top_docking & top_ml)
    dock_scout = len(top_docking & top_composite)
    ml_scout = len(top_ml & top_composite)

    print(f"{top_n:<10} {dock_ml:<12} {dock_scout:<12} {ml_scout:<12}")

# Compare rankings
print("\n\nRank Correlation Analysis:")
from scipy.stats import spearmanr

r_dock_ml, _ = spearmanr(df['rank_docking'], df['rank_ml'])
r_dock_scout, _ = spearmanr(df['rank_docking'], df['rank_composite'])
r_ml_scout, _ = spearmanr(df['rank_ml'], df['rank_composite'])

print(f"  Docking vs ML: Spearman r = {r_dock_ml:.4f}")
print(f"  Docking vs NOD2-Scout: Spearman r = {r_dock_scout:.4f}")
print(f"  ML vs NOD2-Scout: Spearman r = {r_ml_scout:.4f}")

# =============================================================================
# STEP 3: SHAP Feature Importance
# =============================================================================
print("\n--- SHAP FEATURE IMPORTANCE ---")

# Load model
model = xgb.XGBClassifier()
model.load_model(str(MODELS_DIR / 'xgboost_model.json'))

# Load feature names
with open(MODELS_DIR / 'feature_names.pkl', 'rb') as f:
    feature_names = pickle.load(f)

# Prepare feature matrix
descriptor_cols = ['mw', 'logp', 'tpsa', 'hbd', 'hba', 'rotatable_bonds',
                   'aromatic_rings', 'fraction_csp3', 'num_heavy_atoms', 'num_rings',
                   'qed', 'lipinski_violations', 'veber_violations', 'pains_alerts',
                   'permeability_score', 'pgp_risk']

X_desc = df[descriptor_cols].values
X = np.hstack([X_desc, fingerprints])

# Calculate SHAP values (sample for speed)
np.random.seed(42)
sample_idx = np.random.choice(len(X), min(500, len(X)), replace=False)
X_sample = X[sample_idx]

explainer = shap.TreeExplainer(model)
shap_values = explainer.shap_values(X_sample)

# Get feature importance from SHAP
shap_importance = np.abs(shap_values).mean(axis=0)

# Top 20 most important features
top_indices = np.argsort(shap_importance)[-20:][::-1]
print("\nTop 20 Most Important Features:")
print("-" * 40)
for i, idx in enumerate(top_indices):
    print(f"  {i+1:2d}. {feature_names[idx]:<25} {shap_importance[idx]:.4f}")

# Save SHAP importance
shap_df = pd.DataFrame({
    'feature': feature_names,
    'importance': shap_importance
}).sort_values('importance', ascending=False)
shap_df.to_csv(RESULTS_DIR / 'shap_importance.csv', index=False)

# Save SHAP values for plotting
np.save(RESULTS_DIR / 'shap_values.npy', shap_values)
np.save(RESULTS_DIR / 'shap_sample_idx.npy', sample_idx)

print(f"\nSHAP importance saved: {RESULTS_DIR / 'shap_importance.csv'}")

# =============================================================================
# STEP 4: Validate with Known IBD Drugs
# =============================================================================
print("\n--- VALIDATION WITH KNOWN IBD DRUGS ---")

# Known IBD-related drugs to look for
known_ibd_drugs = [
    'CYCLOSPORINE', 'CYCLOSPORIN', 'CICLOSPORIN',
    'BETAMETHASONE',
    'BUDESONIDE',
    'MESALAZINE', 'MESALAMINE',
    'SULFASALAZINE',
    'AZATHIOPRINE',
    'METHOTREXATE',
    'INFLIXIMAB',
    'ADALIMUMAB',
    'VEDOLIZUMAB',
    'PREDNISONE', 'PREDNISOLONE',
    'HYDROCORTISONE',
    'DEXAMETHASONE',
    'TACROLIMUS',
    'THALIDOMIDE',
    'METRONIDAZOLE',
    'CIPROFLOXACIN'
]

# Find matches
df['drug_name_upper'] = df['drug_name'].fillna('').str.upper()

found_drugs = []
for drug in known_ibd_drugs:
    matches = df[df['drug_name_upper'].str.contains(drug, na=False)]
    if len(matches) > 0:
        for _, row in matches.iterrows():
            found_drugs.append({
                'drug_name': row['drug_name'],
                'compound_id': row['compound_id'],
                'cnn_affinity': row['cnn_affinity'],
                'ml_score': row['ml_score'],
                'nod2_scout_score': row['nod2_scout_score'],
                'rank_docking': int(row['rank_docking']),
                'rank_ml': int(row['rank_ml']),
                'rank_composite': int(row['rank_composite'])
            })

if found_drugs:
    ibd_df = pd.DataFrame(found_drugs)
    ibd_df = ibd_df.sort_values('nod2_scout_score', ascending=False)

    print("\nKnown IBD Drugs Found in Dataset:")
    print("-" * 90)
    print(f"{'Drug':<25} {'CNN Aff':<10} {'ML Score':<10} {'Scout':<10} {'Rank':<10}")
    print("-" * 90)

    for _, row in ibd_df.iterrows():
        print(f"{row['drug_name'][:24]:<25} {row['cnn_affinity']:.3f}     {row['ml_score']:.3f}     {row['nod2_scout_score']:.3f}     {int(row['rank_composite']):<10}")

    # Save IBD validation results
    ibd_df.to_csv(RESULTS_DIR / 'ibd_drug_validation.csv', index=False)
    print(f"\nIBD validation saved: {RESULTS_DIR / 'ibd_drug_validation.csv'}")

    # Calculate percentile rankings
    print("\nIBD Drug Percentile Rankings:")
    total_compounds = len(df)
    for _, row in ibd_df.iterrows():
        percentile = 100 * (1 - row['rank_composite'] / total_compounds)
        print(f"  {row['drug_name']}: Top {percentile:.1f}%")
else:
    print("No known IBD drugs found in dataset")

# =============================================================================
# STEP 5: Final Rankings
# =============================================================================
print("\n--- FINAL CANDIDATE RANKINGS ---")

# Get top 50 by NOD2-Scout composite score
top50 = df.nsmallest(50, 'rank_composite')[
    ['drug_name', 'compound_id', 'cnn_affinity', 'ml_score', 'qed',
     'nod2_scout_score', 'rank_composite', 'pains_alerts', 'lipinski_violations']
].copy()

print("\nTop 50 Drug Candidates by NOD2-Scout Score:")
print("-" * 100)
print(f"{'Rank':<6} {'Drug':<30} {'CNN Aff':<10} {'ML':<8} {'QED':<8} {'Score':<10} {'PAINS':<8} {'Lipinski':<8}")
print("-" * 100)

for idx, row in top50.iterrows():
    name = str(row['drug_name'])[:29] if pd.notna(row['drug_name']) else row['compound_id'][:29]
    print(f"{int(row['rank_composite']):<6} {name:<30} {row['cnn_affinity']:.3f}     {row['ml_score']:.3f}   {row['qed']:.3f}   {row['nod2_scout_score']:.4f}    {int(row['pains_alerts']):<8} {int(row['lipinski_violations']):<8}")

# Save final rankings
df_ranked = df.sort_values('rank_composite')
df_ranked.to_csv(RESULTS_DIR / 'final_rankings.csv', index=False)
print(f"\nFinal rankings saved: {RESULTS_DIR / 'final_rankings.csv'}")

# =============================================================================
# STEP 6: Ablation Summary Statistics
# =============================================================================
print("\n--- ABLATION STUDY SUMMARY ---")

# For each method, calculate mean score metrics for top-N
methods = {
    'Docking Only': 'rank_docking',
    'ML Only': 'rank_ml',
    'NOD2-Scout': 'rank_composite'
}

print("\nMean QED and Safety for Top 100 by Each Method:")
print("-" * 60)
print(f"{'Method':<20} {'Mean QED':<12} {'Mean Safety':<12} {'PAINS-free %':<12}")
print("-" * 60)

for method_name, rank_col in methods.items():
    top100 = df.nsmallest(100, rank_col)
    mean_qed = top100['qed'].mean()
    mean_safety = top100['safety_norm'].mean()
    pains_free_pct = (top100['pains_alerts'] == 0).mean() * 100
    print(f"{method_name:<20} {mean_qed:.4f}       {mean_safety:.4f}       {pains_free_pct:.1f}%")

# =============================================================================
# Save Ablation Results
# =============================================================================
ablation_results = {
    'weights': weights,
    'rank_correlations': {
        'docking_vs_ml': float(r_dock_ml),
        'docking_vs_scout': float(r_dock_scout),
        'ml_vs_scout': float(r_ml_scout)
    },
    'top_n_overlap': {}
}

for top_n in top_n_values:
    top_docking = set(df.nsmallest(top_n, 'rank_docking')['compound_id'])
    top_ml = set(df.nsmallest(top_n, 'rank_ml')['compound_id'])
    top_composite = set(df.nsmallest(top_n, 'rank_composite')['compound_id'])

    ablation_results['top_n_overlap'][top_n] = {
        'dock_ml': len(top_docking & top_ml),
        'dock_scout': len(top_docking & top_composite),
        'ml_scout': len(top_ml & top_composite)
    }

with open(RESULTS_DIR / 'ablation_results.pkl', 'wb') as f:
    pickle.dump(ablation_results, f)

print(f"\nAblation results saved: {RESULTS_DIR / 'ablation_results.pkl'}")

print("\n" + "=" * 70)
print("COMPOSITE SCORE & ABLATION COMPLETE")
print("=" * 70)
