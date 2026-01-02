#!/usr/bin/env python3
"""
NOD2-Scout: Inference Script for New Compound Scoring

Usage:
    python predict.py --smiles "CCO"
    python predict.py --file compounds.csv --smiles_col smiles
    python predict.py --file compounds.sdf

Output includes:
    - ML binding probability score (0-1)
    - Composite NOD2-Scout drugability score
    - Drug-likeness metrics (QED, Lipinski, PAINS)
"""

import argparse
import pandas as pd
import numpy as np
from pathlib import Path
import sys
import warnings
warnings.filterwarnings('ignore')

# Check dependencies
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors, QED
    from rdkit.Chem import rdMolDescriptors, FilterCatalog
    from rdkit import RDLogger
    RDLogger.DisableLog('rdApp.*')
except ImportError:
    print("ERROR: RDKit is required. Install with: conda install -c conda-forge rdkit")
    sys.exit(1)

try:
    import xgboost as xgb
except ImportError:
    print("ERROR: XGBoost is required. Install with: pip install xgboost")
    sys.exit(1)

# Paths
SCRIPT_DIR = Path(__file__).resolve().parent
MODELS_DIR = SCRIPT_DIR / 'models'

def load_model():
    """Load trained XGBoost model"""
    model_path = MODELS_DIR / 'xgboost_model.json'
    if not model_path.exists():
        print(f"ERROR: Model not found at {model_path}")
        print("Run the training pipeline first: python scripts/03_train_model.py")
        sys.exit(1)

    model = xgb.XGBClassifier()
    model.load_model(str(model_path))
    return model

def compute_features(smiles, pains_catalog):
    """Compute molecular features for a single SMILES"""
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

        # Veber violations
        veber_violations = 0
        if features['rotatable_bonds'] > 10: veber_violations += 1
        if features['tpsa'] > 140: veber_violations += 1
        features['veber_violations'] = veber_violations

        # PAINS alerts
        pains_matches = pains_catalog.GetMatches(mol)
        features['pains_alerts'] = len(pains_matches)

        # Permeability and P-gp risk
        features['permeability_score'] = 1 - min(features['tpsa'] / 140, 1)
        features['pgp_risk'] = int(features['mw'] > 400 and features['aromatic_rings'] >= 2)

        return features

    except Exception as e:
        return None

def compute_composite_score(ml_score, features, cnn_affinity=None):
    """Compute NOD2-Scout composite score"""
    # Weights
    weights = {
        'ml_score': 0.35,
        'cnn_affinity': 0.30,
        'qed': 0.15,
        'permeability': 0.10,
        'safety': 0.10
    }

    # Normalize components
    ml_norm = ml_score

    # If no docking score provided, use ML score as proxy
    if cnn_affinity is not None:
        # Normalize to 0-1 (assuming typical range 2-8 kcal/mol)
        cnn_norm = min(max((cnn_affinity - 2) / 6, 0), 1)
    else:
        cnn_norm = ml_score  # Use ML as proxy

    qed_norm = features['qed']
    perm_norm = features['permeability_score']

    # Safety score
    pains_score = 1.0 if features['pains_alerts'] == 0 else 0.0
    lipinski_score = 1.0 if features['lipinski_violations'] <= 1 else 0.0
    safety_norm = 0.5 * pains_score + 0.5 * lipinski_score

    # Composite score
    composite = (
        weights['ml_score'] * ml_norm +
        weights['cnn_affinity'] * cnn_norm +
        weights['qed'] * qed_norm +
        weights['permeability'] * perm_norm +
        weights['safety'] * safety_norm
    )

    return composite

def predict_single(smiles, model, pains_catalog, cnn_affinity=None, verbose=True):
    """Predict scores for a single SMILES"""
    features = compute_features(smiles, pains_catalog)

    if features is None:
        if verbose:
            print(f"ERROR: Could not parse SMILES: {smiles}")
        return None

    # Prepare feature vector
    descriptor_cols = ['mw', 'logp', 'tpsa', 'hbd', 'hba', 'rotatable_bonds',
                       'aromatic_rings', 'fraction_csp3', 'num_heavy_atoms', 'num_rings',
                       'qed', 'lipinski_violations', 'veber_violations', 'pains_alerts',
                       'permeability_score', 'pgp_risk']

    X_desc = np.array([[features[col] for col in descriptor_cols]])
    X_fp = features['fingerprint'].reshape(1, -1)
    X = np.hstack([X_desc, X_fp])

    # Predict
    ml_score = model.predict_proba(X)[0, 1]

    # Compute composite score
    composite = compute_composite_score(ml_score, features, cnn_affinity)

    result = {
        'smiles': smiles,
        'ml_score': float(ml_score),
        'nod2_scout_score': float(composite),
        'qed': float(features['qed']),
        'mw': float(features['mw']),
        'logp': float(features['logp']),
        'tpsa': float(features['tpsa']),
        'lipinski_violations': int(features['lipinski_violations']),
        'pains_alerts': int(features['pains_alerts']),
        'permeability_score': float(features['permeability_score'])
    }

    if cnn_affinity is not None:
        result['cnn_affinity'] = cnn_affinity

    return result

def predict_batch(df, smiles_col, model, pains_catalog, cnn_col=None, verbose=True):
    """Predict scores for multiple compounds"""
    results = []

    for idx, row in df.iterrows():
        smiles = row[smiles_col]
        cnn_affinity = row[cnn_col] if cnn_col and cnn_col in df.columns else None

        result = predict_single(smiles, model, pains_catalog, cnn_affinity, verbose=False)

        if result is not None:
            if 'name' in df.columns:
                result['name'] = row['name']
            elif 'compound_id' in df.columns:
                result['name'] = row['compound_id']
            results.append(result)
        else:
            if verbose:
                print(f"  Skipping invalid SMILES at row {idx}")

        if verbose and (idx + 1) % 100 == 0:
            print(f"  Processed {idx + 1}/{len(df)} compounds...")

    return pd.DataFrame(results)

def load_sdf(sdf_path):
    """Load compounds from SDF file"""
    from rdkit import Chem

    suppl = Chem.SDMolSupplier(str(sdf_path))
    data = []

    for mol in suppl:
        if mol is None:
            continue

        smiles = Chem.MolToSmiles(mol)
        name = mol.GetProp('_Name') if mol.HasProp('_Name') else None

        data.append({'smiles': smiles, 'name': name})

    return pd.DataFrame(data)

def main():
    parser = argparse.ArgumentParser(
        description='NOD2-Scout: Score compounds for NOD2 binding potential',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python predict.py --smiles "CC(=O)Nc1ccc(O)cc1"
  python predict.py --file my_compounds.csv --smiles_col SMILES
  python predict.py --file library.sdf --output results.csv
        """
    )

    parser.add_argument('--smiles', type=str, help='Single SMILES string to score')
    parser.add_argument('--file', type=str, help='Input file (CSV or SDF)')
    parser.add_argument('--smiles_col', type=str, default='smiles',
                        help='Column name for SMILES in CSV (default: smiles)')
    parser.add_argument('--cnn_col', type=str,
                        help='Column name for CNN affinity scores (optional)')
    parser.add_argument('--output', type=str, help='Output CSV file path')
    parser.add_argument('--top_n', type=int, default=50,
                        help='Show top N results (default: 50)')

    args = parser.parse_args()

    if not args.smiles and not args.file:
        parser.print_help()
        print("\nERROR: Please provide --smiles or --file")
        sys.exit(1)

    # Load model
    print("Loading NOD2-Scout model...")
    model = load_model()

    # Initialize PAINS filter
    pains_params = FilterCatalog.FilterCatalogParams()
    pains_params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
    pains_catalog = FilterCatalog.FilterCatalog(pains_params)

    if args.smiles:
        # Single SMILES prediction
        print(f"\nScoring: {args.smiles}")
        print("=" * 60)

        result = predict_single(args.smiles, model, pains_catalog)

        if result:
            print(f"\n  ML Binding Score:      {result['ml_score']:.4f}")
            print(f"  NOD2-Scout Score:      {result['nod2_scout_score']:.4f}")
            print(f"\n  Drug-likeness (QED):   {result['qed']:.4f}")
            print(f"  Molecular Weight:      {result['mw']:.1f}")
            print(f"  LogP:                  {result['logp']:.2f}")
            print(f"  TPSA:                  {result['tpsa']:.1f}")
            print(f"  Lipinski Violations:   {result['lipinski_violations']}")
            print(f"  PAINS Alerts:          {result['pains_alerts']}")
            print(f"  Permeability Score:    {result['permeability_score']:.4f}")

            # Interpretation
            print("\n" + "-" * 60)
            if result['nod2_scout_score'] >= 0.75:
                print("  Interpretation: HIGH potential - Priority candidate")
            elif result['nod2_scout_score'] >= 0.50:
                print("  Interpretation: MODERATE potential - Worth investigating")
            else:
                print("  Interpretation: LOW potential - Deprioritize")

    else:
        # Batch prediction
        input_path = Path(args.file)

        if not input_path.exists():
            print(f"ERROR: File not found: {input_path}")
            sys.exit(1)

        print(f"\nLoading compounds from: {input_path}")

        if input_path.suffix.lower() == '.sdf':
            df = load_sdf(input_path)
            smiles_col = 'smiles'
        else:
            df = pd.read_csv(input_path)
            smiles_col = args.smiles_col

        print(f"Found {len(df)} compounds")
        print("\nProcessing...")

        results_df = predict_batch(df, smiles_col, model, pains_catalog,
                                   args.cnn_col, verbose=True)

        # Sort by NOD2-Scout score
        results_df = results_df.sort_values('nod2_scout_score', ascending=False)

        # Save results
        if args.output:
            output_path = Path(args.output)
            results_df.to_csv(output_path, index=False)
            print(f"\nResults saved to: {output_path}")

        # Show top results
        print(f"\n{'=' * 80}")
        print(f"TOP {min(args.top_n, len(results_df))} CANDIDATES BY NOD2-SCOUT SCORE")
        print(f"{'=' * 80}")
        print(f"{'Rank':<6} {'Name':<30} {'ML':<8} {'Scout':<8} {'QED':<8} {'PAINS':<6}")
        print("-" * 80)

        for i, (_, row) in enumerate(results_df.head(args.top_n).iterrows()):
            name = str(row.get('name', row['smiles'][:30]))[:29]
            print(f"{i+1:<6} {name:<30} {row['ml_score']:.4f}  {row['nod2_scout_score']:.4f}  {row['qed']:.4f}  {int(row['pains_alerts']):<6}")

        print(f"\nTotal compounds scored: {len(results_df)}")
        print(f"High potential (score >= 0.75): {(results_df['nod2_scout_score'] >= 0.75).sum()}")
        print(f"Moderate potential (0.50-0.75): {((results_df['nod2_scout_score'] >= 0.50) & (results_df['nod2_scout_score'] < 0.75)).sum()}")

if __name__ == '__main__':
    main()
