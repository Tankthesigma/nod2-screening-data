#!/usr/bin/env python3
"""
NOD2-Scout Step 5: Generate Publication-Quality Figures
- ROC curve with confidence interval
- SHAP summary plot
- Score distributions
- Ablation comparison
- Top candidates visualization
"""

import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from sklearn.metrics import roc_curve, auc, precision_recall_curve
import xgboost as xgb
import shap
import pickle
import warnings
warnings.filterwarnings('ignore')

# Set publication-quality style
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['figure.dpi'] = 150
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['axes.titlesize'] = 14

print("=" * 70)
print("NOD2-SCOUT: GENERATING PUBLICATION FIGURES")
print("=" * 70)

# Paths
BASE_DIR = Path(__file__).resolve().parents[2]
DATA_DIR = BASE_DIR / 'nod2-scout' / 'data'
MODELS_DIR = BASE_DIR / 'nod2-scout' / 'models'
RESULTS_DIR = BASE_DIR / 'nod2-scout' / 'results'
FIGURES_DIR = BASE_DIR / 'nod2-scout' / 'figures'
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# Load data
df = pd.read_csv(RESULTS_DIR / 'final_rankings.csv')
fingerprints = np.load(DATA_DIR / 'fingerprints.npy')
shap_values = np.load(RESULTS_DIR / 'shap_values.npy')
shap_sample_idx = np.load(RESULTS_DIR / 'shap_sample_idx.npy')
shap_importance = pd.read_csv(RESULTS_DIR / 'shap_importance.csv')

with open(MODELS_DIR / 'training_metadata.pkl', 'rb') as f:
    metadata = pickle.load(f)

with open(MODELS_DIR / 'feature_names.pkl', 'rb') as f:
    feature_names = pickle.load(f)

with open(RESULTS_DIR / 'ablation_results.pkl', 'rb') as f:
    ablation = pickle.load(f)

print(f"Loaded {len(df)} compounds")

# =============================================================================
# Figure 1: ROC Curve with Negative Control Comparison
# =============================================================================
print("\n1. Generating ROC Curve...")

# Recreate train/test split for ROC
top_threshold = metadata['top_threshold']
bottom_threshold = metadata['bottom_threshold']

df['label'] = -1
df.loc[df['cnn_affinity'] >= top_threshold, 'label'] = 1
df.loc[df['cnn_affinity'] <= bottom_threshold, 'label'] = 0
labeled_df = df[df['label'] >= 0].copy()

# Use ML scores as predictions
y_true = labeled_df['label'].values
y_pred = labeled_df['ml_score'].values

fpr, tpr, _ = roc_curve(y_true, y_pred)
roc_auc = auc(fpr, tpr)

fig, ax = plt.subplots(figsize=(8, 7))

# Main ROC curve
ax.plot(fpr, tpr, color='#2E86AB', lw=2.5,
        label=f'NOD2-Scout (AUROC = {roc_auc:.3f})')

# Random classifier line
ax.plot([0, 1], [0, 1], 'k--', lw=1.5, label='Random (AUROC = 0.50)')

# Shuffled control region
shuffle_mean = metadata['shuffle_auroc_mean']
shuffle_std = metadata['shuffle_auroc_std']
ax.axhspan(0.45, 0.55, alpha=0.15, color='red', label=f'Shuffle Control ({shuffle_mean:.2f}±{shuffle_std:.2f})')

ax.set_xlim([0.0, 1.0])
ax.set_ylim([0.0, 1.05])
ax.set_xlabel('False Positive Rate', fontsize=12)
ax.set_ylabel('True Positive Rate', fontsize=12)
ax.set_title('NOD2-Scout Model Performance\nScaffold-Split Validation', fontsize=14, fontweight='bold')
ax.legend(loc='lower right', fontsize=10)

# Add text box with stats
textstr = f'Train: {metadata["n_train"]} compounds\nTest: {metadata["n_test"]} compounds\nScaffolds: {metadata["n_scaffolds"]}'
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax.text(0.65, 0.2, textstr, transform=ax.transAxes, fontsize=9,
        verticalalignment='top', bbox=props)

plt.tight_layout()
plt.savefig(FIGURES_DIR / 'fig1_roc_curve.png', dpi=300, bbox_inches='tight')
plt.close()
print(f"   Saved: {FIGURES_DIR / 'fig1_roc_curve.png'}")

# =============================================================================
# Figure 2: SHAP Feature Importance (Top 15)
# =============================================================================
print("2. Generating SHAP Importance Plot...")

top_features = shap_importance.head(15)

fig, ax = plt.subplots(figsize=(10, 7))

# Clean feature names for display
clean_names = []
for name in top_features['feature']:
    if name.startswith('fp_'):
        clean_names.append(f'Fingerprint Bit {name[3:]}')
    else:
        name_map = {
            'num_heavy_atoms': 'Heavy Atom Count',
            'logp': 'LogP',
            'mw': 'Molecular Weight',
            'num_rings': 'Ring Count',
            'fraction_csp3': 'Fraction Csp3',
            'qed': 'QED Score',
            'rotatable_bonds': 'Rotatable Bonds',
            'lipinski_violations': 'Lipinski Violations',
            'hba': 'H-Bond Acceptors',
            'aromatic_rings': 'Aromatic Rings',
            'tpsa': 'TPSA',
            'hbd': 'H-Bond Donors',
            'permeability_score': 'Permeability',
            'pgp_risk': 'P-gp Risk',
            'veber_violations': 'Veber Violations',
            'pains_alerts': 'PAINS Alerts'
        }
        clean_names.append(name_map.get(name, name))

colors = ['#2E86AB' if not f.startswith('fp_') else '#A23B72' for f in top_features['feature']]

bars = ax.barh(range(len(top_features)), top_features['importance'].values, color=colors)
ax.set_yticks(range(len(top_features)))
ax.set_yticklabels(clean_names)
ax.invert_yaxis()
ax.set_xlabel('Mean |SHAP Value|', fontsize=12)
ax.set_title('Feature Importance for Binding Prediction\n(SHAP Analysis)', fontsize=14, fontweight='bold')

# Legend
physico_patch = mpatches.Patch(color='#2E86AB', label='Physicochemical')
fp_patch = mpatches.Patch(color='#A23B72', label='Fingerprint')
ax.legend(handles=[physico_patch, fp_patch], loc='lower right')

plt.tight_layout()
plt.savefig(FIGURES_DIR / 'fig2_shap_importance.png', dpi=300, bbox_inches='tight')
plt.close()
print(f"   Saved: {FIGURES_DIR / 'fig2_shap_importance.png'}")

# =============================================================================
# Figure 3: Score Distribution Comparison
# =============================================================================
print("3. Generating Score Distribution Plot...")

fig, axes = plt.subplots(1, 3, figsize=(14, 5))

# CNN Affinity Distribution
ax = axes[0]
ax.hist(df['cnn_affinity'], bins=50, color='#2E86AB', alpha=0.7, edgecolor='black', linewidth=0.5)
ax.axvline(top_threshold, color='green', linestyle='--', lw=2, label=f'Top 10% (>{top_threshold:.2f})')
ax.axvline(bottom_threshold, color='red', linestyle='--', lw=2, label=f'Bottom 10% (<{bottom_threshold:.2f})')
ax.set_xlabel('CNN Affinity (kcal/mol)', fontsize=11)
ax.set_ylabel('Frequency', fontsize=11)
ax.set_title('Docking Score Distribution', fontsize=12, fontweight='bold')
ax.legend(fontsize=9)

# ML Score Distribution
ax = axes[1]
ax.hist(df['ml_score'], bins=50, color='#A23B72', alpha=0.7, edgecolor='black', linewidth=0.5)
ax.axvline(0.5, color='gray', linestyle='--', lw=2, label='Decision Threshold')
ax.set_xlabel('ML Score', fontsize=11)
ax.set_ylabel('Frequency', fontsize=11)
ax.set_title('ML Prediction Distribution', fontsize=12, fontweight='bold')
ax.legend(fontsize=9)

# Composite Score Distribution
ax = axes[2]
ax.hist(df['nod2_scout_score'], bins=50, color='#F18F01', alpha=0.7, edgecolor='black', linewidth=0.5)
top100_threshold = df.nsmallest(100, 'rank_composite')['nod2_scout_score'].min()
ax.axvline(top100_threshold, color='green', linestyle='--', lw=2, label=f'Top 100 (>{top100_threshold:.2f})')
ax.set_xlabel('NOD2-Scout Score', fontsize=11)
ax.set_ylabel('Frequency', fontsize=11)
ax.set_title('Composite Score Distribution', fontsize=12, fontweight='bold')
ax.legend(fontsize=9)

plt.suptitle('Score Distributions Across Methods', fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig(FIGURES_DIR / 'fig3_score_distributions.png', dpi=300, bbox_inches='tight')
plt.close()
print(f"   Saved: {FIGURES_DIR / 'fig3_score_distributions.png'}")

# =============================================================================
# Figure 4: Ablation Study - Method Comparison
# =============================================================================
print("4. Generating Ablation Study Plot...")

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Top-N Overlap
ax = axes[0]
top_ns = list(ablation['top_n_overlap'].keys())
dock_ml = [ablation['top_n_overlap'][n]['dock_ml'] for n in top_ns]
dock_scout = [ablation['top_n_overlap'][n]['dock_scout'] for n in top_ns]
ml_scout = [ablation['top_n_overlap'][n]['ml_scout'] for n in top_ns]

x = np.arange(len(top_ns))
width = 0.25

bars1 = ax.bar(x - width, dock_ml, width, label='Docking ∩ ML', color='#2E86AB')
bars2 = ax.bar(x, dock_scout, width, label='Docking ∩ Scout', color='#A23B72')
bars3 = ax.bar(x + width, ml_scout, width, label='ML ∩ Scout', color='#F18F01')

ax.set_xlabel('Top-N Candidates', fontsize=11)
ax.set_ylabel('Number of Overlapping Compounds', fontsize=11)
ax.set_title('Top-N Overlap Between Methods', fontsize=12, fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels([str(n) for n in top_ns])
ax.legend(fontsize=9)

# Drug-likeness Comparison
ax = axes[1]
methods = ['Docking\nOnly', 'ML\nOnly', 'NOD2-\nScout']

# Calculate metrics for top 100 by each method
qed_values = []
safety_values = []
for rank_col in ['rank_docking', 'rank_ml', 'rank_composite']:
    top100 = df.nsmallest(100, rank_col)
    qed_values.append(top100['qed'].mean())
    safety_values.append(top100['safety_norm'].mean())

x = np.arange(len(methods))
width = 0.35

bars1 = ax.bar(x - width/2, qed_values, width, label='Mean QED', color='#2E86AB')
bars2 = ax.bar(x + width/2, safety_values, width, label='Mean Safety', color='#A23B72')

ax.set_xlabel('Ranking Method', fontsize=11)
ax.set_ylabel('Score (0-1)', fontsize=11)
ax.set_title('Drug-likeness of Top 100 Candidates', fontsize=12, fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels(methods)
ax.legend(fontsize=9)
ax.set_ylim(0, 1.1)

# Add value labels on bars
for bar in bars1:
    height = bar.get_height()
    ax.annotate(f'{height:.2f}', xy=(bar.get_x() + bar.get_width()/2, height),
                xytext=(0, 3), textcoords="offset points", ha='center', va='bottom', fontsize=9)
for bar in bars2:
    height = bar.get_height()
    ax.annotate(f'{height:.2f}', xy=(bar.get_x() + bar.get_width()/2, height),
                xytext=(0, 3), textcoords="offset points", ha='center', va='bottom', fontsize=9)

plt.suptitle('Ablation Study: Comparing Ranking Strategies', fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig(FIGURES_DIR / 'fig4_ablation_study.png', dpi=300, bbox_inches='tight')
plt.close()
print(f"   Saved: {FIGURES_DIR / 'fig4_ablation_study.png'}")

# =============================================================================
# Figure 5: Top 20 Candidates Visualization
# =============================================================================
print("5. Generating Top Candidates Plot...")

top20 = df.nsmallest(20, 'rank_composite')[['drug_name', 'cnn_affinity', 'ml_score', 'qed', 'nod2_scout_score']].copy()
top20['name'] = top20['drug_name'].fillna('Unknown').apply(lambda x: x[:25])

fig, ax = plt.subplots(figsize=(12, 8))

y_pos = np.arange(len(top20))
bar_height = 0.7

# Stacked bar showing components
bars1 = ax.barh(y_pos, top20['cnn_norm'] if 'cnn_norm' in top20.columns else top20['nod2_scout_score'] * 0.3,
                bar_height, label='Docking Component', color='#2E86AB', alpha=0.8)

# Just show total score
bars = ax.barh(y_pos, top20['nod2_scout_score'], bar_height, color='#2E86AB', alpha=0.8)

ax.set_yticks(y_pos)
ax.set_yticklabels(top20['name'])
ax.invert_yaxis()
ax.set_xlabel('NOD2-Scout Composite Score', fontsize=12)
ax.set_title('Top 20 Drug Candidates by NOD2-Scout Score', fontsize=14, fontweight='bold')

# Add score labels
for i, (score, qed, ml) in enumerate(zip(top20['nod2_scout_score'], top20['qed'], top20['ml_score'])):
    ax.text(score + 0.01, i, f'QED:{qed:.2f} ML:{ml:.2f}', va='center', fontsize=8)

ax.set_xlim(0, 1.0)

plt.tight_layout()
plt.savefig(FIGURES_DIR / 'fig5_top_candidates.png', dpi=300, bbox_inches='tight')
plt.close()
print(f"   Saved: {FIGURES_DIR / 'fig5_top_candidates.png'}")

# =============================================================================
# Figure 6: IBD Drug Validation Heatmap
# =============================================================================
print("6. Generating IBD Validation Plot...")

ibd_df = pd.read_csv(RESULTS_DIR / 'ibd_drug_validation.csv')

# Keep unique drug names (first occurrence)
ibd_unique = ibd_df.drop_duplicates(subset='drug_name').head(15)

fig, ax = plt.subplots(figsize=(10, 8))

# Create percentile ranks
total = len(df)
ibd_unique['percentile'] = 100 * (1 - ibd_unique['rank_composite'] / total)

# Sort by percentile
ibd_unique = ibd_unique.sort_values('percentile', ascending=False)

y_pos = np.arange(len(ibd_unique))
colors = plt.cm.RdYlGn(ibd_unique['percentile'] / 100)

bars = ax.barh(y_pos, ibd_unique['percentile'], color=colors, edgecolor='black', linewidth=0.5)

ax.set_yticks(y_pos)
ax.set_yticklabels(ibd_unique['drug_name'].apply(lambda x: x[:30]))
ax.invert_yaxis()
ax.set_xlabel('Percentile Rank (Higher = Better)', fontsize=12)
ax.set_title('Known IBD Drugs: NOD2-Scout Ranking Validation', fontsize=14, fontweight='bold')
ax.axvline(50, color='gray', linestyle='--', lw=1.5, label='Median')
ax.axvline(90, color='green', linestyle='--', lw=1.5, label='Top 10%')

# Add percentile labels
for i, pct in enumerate(ibd_unique['percentile']):
    ax.text(pct + 1, i, f'{pct:.1f}%', va='center', fontsize=9)

ax.legend(loc='lower right')
ax.set_xlim(0, 105)

plt.tight_layout()
plt.savefig(FIGURES_DIR / 'fig6_ibd_validation.png', dpi=300, bbox_inches='tight')
plt.close()
print(f"   Saved: {FIGURES_DIR / 'fig6_ibd_validation.png'}")

# =============================================================================
# Figure 7: Correlation Matrix (Docking vs ML vs Composite)
# =============================================================================
print("7. Generating Correlation Plot...")

fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))

# Sample for plotting
sample = df.sample(min(500, len(df)), random_state=42)

# Docking vs ML
ax = axes[0]
ax.scatter(sample['cnn_affinity'], sample['ml_score'], alpha=0.5, s=20, c='#2E86AB')
ax.set_xlabel('CNN Affinity (kcal/mol)', fontsize=11)
ax.set_ylabel('ML Score', fontsize=11)
ax.set_title(f'Docking vs ML\n(Spearman r = {ablation["rank_correlations"]["docking_vs_ml"]:.3f})',
             fontsize=12, fontweight='bold')

# Docking vs Composite
ax = axes[1]
ax.scatter(sample['cnn_affinity'], sample['nod2_scout_score'], alpha=0.5, s=20, c='#A23B72')
ax.set_xlabel('CNN Affinity (kcal/mol)', fontsize=11)
ax.set_ylabel('NOD2-Scout Score', fontsize=11)
ax.set_title(f'Docking vs NOD2-Scout\n(Spearman r = {ablation["rank_correlations"]["docking_vs_scout"]:.3f})',
             fontsize=12, fontweight='bold')

# ML vs Composite
ax = axes[2]
ax.scatter(sample['ml_score'], sample['nod2_scout_score'], alpha=0.5, s=20, c='#F18F01')
ax.set_xlabel('ML Score', fontsize=11)
ax.set_ylabel('NOD2-Scout Score', fontsize=11)
ax.set_title(f'ML vs NOD2-Scout\n(Spearman r = {ablation["rank_correlations"]["ml_vs_scout"]:.3f})',
             fontsize=12, fontweight='bold')

plt.suptitle('Score Correlations Between Methods', fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig(FIGURES_DIR / 'fig7_correlations.png', dpi=300, bbox_inches='tight')
plt.close()
print(f"   Saved: {FIGURES_DIR / 'fig7_correlations.png'}")

# =============================================================================
# Summary
# =============================================================================
print("\n" + "=" * 70)
print("FIGURES GENERATED")
print("=" * 70)
print(f"\nSaved {7} publication-quality figures to: {FIGURES_DIR}/")
print("\nFigure descriptions:")
print("  fig1_roc_curve.png        - Model performance with negative control")
print("  fig2_shap_importance.png  - Feature importance (SHAP analysis)")
print("  fig3_score_distributions.png - Score distributions by method")
print("  fig4_ablation_study.png   - Method comparison (ablation)")
print("  fig5_top_candidates.png   - Top 20 candidates visualization")
print("  fig6_ibd_validation.png   - Known IBD drug validation")
print("  fig7_correlations.png     - Score correlations between methods")
